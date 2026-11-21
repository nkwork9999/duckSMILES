import Descriptors
import Descriptors.DrugLikeness
import Descriptors.Burden
import Descriptors.Morgan
import Descriptors.Information
import Descriptors.Connectivity
import Descriptors.GraphIndices
import Descriptors.FragmentPatterns
open Lean Descriptors

abbrev ReadResult := Except String
abbrev field := Json.getObjVal?

def number? : Json → ReadResult (Option JsonNumber)
  | .null => pure none
  | .num n => pure (some n)
  | _ => throw "expected a number or null"

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

def readList (obs : Json) (α : Type) [FromJson α] (key : String) : ReadResult (List α) := do
  fromJson? (← field obs key)

structure PerceivedRing where
  elements : List String
  orders : List String
  deriving FromJson

def ringDescriptors (rings : List PerceivedRing) : List (String × Json) := Id.run do
  let aromatic := fun (r : PerceivedRing) => r.orders.all (· == "AROMATIC")
  let saturated := fun (r : PerceivedRing) => r.orders.all (· == "SINGLE")
  let hetero := fun (r : PerceivedRing) => r.elements.any (fun e => e != "C" && e != "H")
  let categories := [("Aromatic", aromatic), ("Aliphatic", fun r => !(aromatic r)),
                     ("Saturated", saturated)]
  let mut counts := [("NumHeterocycles", toJson (classCounts rings hetero).1)]
  for (name, predicate) in categories do
    let members := rings.filter predicate
    let (h, c) := classCounts members hetero
    counts := counts ++ [("Num" ++ name ++ "Rings", toJson members.length),
                         ("Num" ++ name ++ "Heterocycles", toJson h),
                         ("Num" ++ name ++ "Carbocycles", toJson c)]
  return counts

def reference (cols : Array Column) (obs : Json) : ReadResult (Array Json) := do
  let atoms ← fromJson? (α := List String) (← field obs "atoms")
  let vertices ← fromJson? (α := Array Fragment.Atom) (← field obs "vertices")
  let edges ← readList obs (List Nat) "edges"
  let graphBonds ← edges.mapM fun edge => do
    let [u, v, order] := edge | throw "invalid bond record"
    if u >= vertices.size || v >= vertices.size || u == v then throw "invalid bond endpoints"
    pure (u, v, order)
  let ringFlags ← readList obs Bool "ringBonds"
  unless ringFlags.length == graphBonds.length do throw "ring bond flag length mismatch"
  let ringEdges := (graphBonds.zip ringFlags).filterMap fun ((i,j,_), ring) =>
    if ring then some (min i j, max i j) else none
  let g := Fragment.Graph.ofBonds vertices graphBonds ringEdges
  let charges ← readList obs Int "charges"
  let eh ← readList obs Nat "explicitH"
  let ih ← readList obs Nat "implicitH"
  let bonds ← readList obs String "bondOrders"
  let aromatic ← readList obs Bool "aromaticAtoms"
  let ringAtoms ← readList obs Bool "ringAtoms"
  let components ← readList obs Nat "componentSizes"
  let rings ← readList obs Nat "ringSizes"
  let perceived ← readList obs PerceivedRing "perceivedRings"
  let fragmentOnly := cols.all fun c => c.id.startsWith "fr_"
  let atomic ← if fragmentOnly then pure [] else do
    let rows ← readList obs Atomic.Contribution "atomicContributions"
    unless rows.length == vertices.size do throw "atomic contribution length mismatch"
    let props ← fromJson? (α := Array Composition.AtomProperties) (← field obs "atomProperties")
    unless props.size == vertices.size do throw "atom property length mismatch"
    let initial ← fromJson? (α := Array Nat) (← field obs "morganInitial")
    unless initial.size == vertices.size do throw "connectivity invariant length mismatch"
    unless initial.all (· < 2^32) do throw "connectivity invariant exceeds 32 bits"
    let lp ← readList obs Nat "expandedLogP"
    let mr ← readList obs Nat "expandedMR"
    unless lp.length == mr.length do throw "expanded Crippen length mismatch"
    let q := ("qed", ← DrugLikeness.value (← readList obs Nat "qedProperties"))
    let burden ← Burden.values (← fromJson? (← field obs "burden"))
    let ipc ← Information.values g
    let ringsByBond ← readList obs (List Nat) "ringBondIndices"
    pure ([q] ++ burden ++ ipc ++ Morgan.values g graphBonds (initial.map Nat.toUInt32) ++
          Atomic.values rows (← fromJson? (← field obs "hydrogenArea")) ++
          Composition.values g props lp mr ++ Connectivity.values g graphBonds props ++
          Composition.ringJunctions ringsByBond graphBonds.toArray)
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
    ("ds_RingBonds", toJson (ringFlags.count true)),
    ("ds_LargestRingSize", toJson (rings.foldl max 0)),
    ("RingCount", toJson rings.length)] ++ ringDescriptors perceived ++ atomic ++
      (if fragmentOnly then [] else GraphIndices.values g)
  cols.mapM fun col => do
    if let some v := totals.lookup col.id then return v
    if ["C","N","O","S","P","F","Cl","Br","I"].any (fun e => col.id == "ds_Count" ++ e) then
      return toJson (countAtoms atoms (col.id.drop 8))
    if let some q := Fragment.pattern? g col.id then
      return toJson (Fragment.fragmentCount g q)
    throw s!"no Lean definition for catalog ID {col.id}"

def verify (cols : Array Column) (actual expected : Array Json) : ReadResult Unit := do
  if actual.size != cols.size || expected.size != cols.size then throw "profile length mismatch"
  let xs ← actual.mapM number?
  let ys ← expected.mapM number?
  if h : checkAll cols.size xs.size ys.size (atColumn cols xs ys) = true then
    have _certificate : ProfileMatches cols xs ys := (verify_profile_iff cols xs ys).mp h
    return ()
  for h : i in [:cols.size] do
    unless atColumn cols xs ys i do
      throw s!"{cols[i].id}: observed {actual[i]!}, expected {expected[i]!}"
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
