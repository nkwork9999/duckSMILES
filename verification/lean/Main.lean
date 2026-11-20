import Descriptors
open Lean Descriptors

abbrev ReadResult := Except String
abbrev field := Json.getObjVal?

def number? : Json → ReadResult (Option JsonNumber)
  | .null => pure none
  | .num n => pure (some n)
  | _ => throw "expected a number or null"

structure Column where
  id : String
  absolute : JsonNumber
  relative : JsonNumber

def catalog (j : Json) : ReadResult (Array Column) := do
  let rows ← (← field j "descriptors").getArr?
  let cols ← rows.mapM fun s => do
    let name ← (← field s "id").getStr?
    let a ← (← field s "tolerance").getNum?
    let r ← (← field s "relative_tolerance").getNum?
    if a.mantissa < 0 || r.mantissa < 0 then throw s!"negative tolerance for {name}"
    pure ({ id := name, absolute := a, relative := r } : Column)
  let names := (cols.map Column.id).toList
  if names.isEmpty || names.eraseDups.length != names.length then throw "invalid descriptor IDs"
  return cols

/-- Evaluate element counts from the graph, fragment cardinalities from
uniquified SMARTS matches, and the other fields from their named oracle. -/
def readList (obs : Json) (α : Type) [FromJson α] (key : String) : ReadResult (List α) := do
  fromJson? (← field obs key)

def reference (cols : Array Column) (obs : Json) : ReadResult (Array Json) := do
  let atoms ← fromJson? (α := List String) (← field obs "atoms")
  let embeddings ← field obs "matches"
  let oracle ← field obs "oracle"
  let charges ← readList obs Int "charges"
  let eh ← readList obs Nat "explicitH"
  let ih ← readList obs Nat "implicitH"
  let bonds ← readList obs String "bondOrders"
  let aromatic ← readList obs Bool "aromaticAtoms"
  let ringAtoms ← readList obs Bool "ringAtoms"
  let ringBonds ← readList obs Bool "ringBonds"
  let components ← readList obs Nat "componentSizes"
  let rings ← readList obs Nat "ringSizes"
  let totals : List (String × Json) := [
    ("HeavyAtomCount", toJson (atoms.length - atoms.count "H")),
    ("NumHeteroatoms", toJson ((atoms.filter fun a => a != "C" && a != "H").length)),
    ("NOCount", toJson (atoms.count "N" + atoms.count "O")),
    ("ds_BondCount", toJson bonds.length),
    ("ds_FormalCharge", toJson charges.sum),
    ("ds_FragmentCount", toJson components.length),
    ("ds_ExplicitHydrogens", toJson eh.sum),
    ("ds_ImplicitHydrogens", toJson ih.sum),
    ("ds_TotalHydrogens", toJson (eh.sum + ih.sum)),
    ("ds_HalogenCount", toJson ((atoms.filter fun a => ["F","Cl","Br","I"].contains a).length)),
    ("ds_SingleBonds", toJson (bonds.count "SINGLE")),
    ("ds_DoubleBonds", toJson (bonds.count "DOUBLE")),
    ("ds_TripleBonds", toJson (bonds.count "TRIPLE")),
    ("ds_AromaticBonds", toJson (bonds.count "AROMATIC")),
    ("ds_AromaticAtoms", toJson (aromatic.count true)),
    ("ds_RingAtoms", toJson (ringAtoms.count true)),
    ("ds_RingBonds", toJson (ringBonds.count true)),
    ("ds_LargestRingSize", toJson (rings.foldl max 0)),
    ("RingCount", toJson rings.length)]
  cols.mapM fun col => do
    if let some v := totals.lookup col.id then return v
    if col.id.startsWith "ds_Count" then
      return toJson (countAtoms atoms (col.id.drop 8))
    if col.id.startsWith "fr_" then
      return toJson (← (← field embeddings col.id).getArr?).size
    return ← field oracle col.id

def verify (cols : Array Column) (actual expected : Array Json) : ReadResult Unit := do
  if actual.size != cols.size || expected.size != cols.size then throw "profile length mismatch"
  let mut results := #[]
  for h : i in [:cols.size] do
    let c := cols[i]
    let ok := agrees c.absolute c.relative (← number? actual[i]!) (← number? expected[i]!)
    if !ok then throw s!"{c.id}: observed {actual[i]!}, expected {expected[i]!}"
    results := results.push ok
  unless checkAll cols.size actual.size expected.size (fun i => results[i]?.getD false) do
    throw "incomplete comparison"

def readJson (p : String) : IO Json := do
  IO.ofExcept <| Json.parse (← IO.FS.readFile p)

def run (args : List String) : IO Unit := do
  let mode :: cat :: ref :: rest := args | throw <| IO.userError "use reference CATALOG OBS or check CATALOG OBS ACTUAL"
  let cols ← IO.ofExcept <| catalog (← readJson cat)
  let observations ← IO.ofExcept (← readJson ref).getArr?
  if observations.isEmpty then throw <| IO.userError "empty input"
  let expected ← IO.ofExcept <| observations.mapM fun row => do
    let id ← (← field row "smiles").getStr?
    let obs ← field row "observations"
    let values ← if obs == Json.null then pure Json.null else Json.arr <$> reference cols obs
    pure (id, values)
  match mode, rest with
  | "reference", [] =>
    IO.println (Json.arr (expected.map fun (s, v) => Json.mkObj [("smiles", toJson s), ("values", v)])).compress
  | "check", [candidate] => do
    let rows ← IO.ofExcept (← readJson candidate).getArr?
    if rows.size != expected.size then throw <| IO.userError "row count mismatch"
    for (row, (s, v)) in rows.zip expected do
      IO.ofExcept <| do
        if (← (← field row "smiles").getStr?) != s then throw "input SMILES or order changed"
        let actual ← field row "values"
        if v == Json.null then
          if actual != Json.null then throw s!"{s}: expected unsupported molecule"
        else
          match verify cols (← actual.getArr?) (← v.getArr?) with
          | .error e => throw s!"{s}: {e}"
          | .ok _ => pure ()
    IO.println s!"Verified {expected.size} rows against {cols.size} descriptor columns."
  | _, _ => throw <| IO.userError "invalid command or argument count"

def main (args : List String) : IO Unit := run args
