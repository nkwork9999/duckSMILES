import Descriptors
import Descriptors.EState
import Descriptors.QED
import Descriptors.Burden
import Descriptors.Morgan
import Descriptors.Information
import Descriptors.Connectivity
import Descriptors.GraphIndices
import Descriptors.ProfileCounts
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
  if cols.all (fun col => col.id == "Ipc" || col.id == "AvgIpc") then
    let values ← Information.values g
    return ← cols.mapM (fun col => match values.lookup col.id with | some val => pure val | none => throw "missing information column")
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
    let lp ← readList obs Nat "expandedLogP"
    let mr ← readList obs Nat "expandedMR"
    unless lp.length == mr.length do throw "expanded Crippen length mismatch"
    let q := ("qed", ← QED.value (← fromJson? (← field obs "qedGraph")))
    let burden ← Burden.values (← fromJson? (← field obs "burden"))
    let ipc ← Information.values g
    let ringsByBond ← readList obs (List Nat) "ringBondIndices"
    pure ([q] ++ burden ++ ipc ++ Morgan.values g graphBonds props ++
          Atomic.values rows (← Labute.areas g graphBonds) (← EState.indices g props) ++
          Composition.values g props lp mr ++ Connectivity.values g graphBonds props ++
          Composition.ringJunctions ringsByBond graphBonds.toArray)
  let counts : ProfileCounts.Input :=
    ⟨atoms, charges, eh, ih, bonds, aromatic, ringAtoms, ringFlags, components, rings⟩
  let extra ← if fragmentOnly then pure [] else GraphIndices.values g
  let totals : List (String × Json) := [
    ("HeavyAtomCount", toJson (atoms.length - atoms.count "H")),
    ("NumHeteroatoms", toJson ((atoms.filter fun a => a != "C" && a != "H").length)),
    ("NOCount", toJson (atoms.count "N" + atoms.count "O")),
    ("RingCount", toJson rings.length)] ++
    ProfileCounts.values counts ++ ringDescriptors perceived ++ atomic ++ extra
  cols.mapM fun col => do
    if let some v := totals.lookup col.id then return v
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
  | mode, [candidate] => do
    unless mode == "check" || mode == "audit" do throw <| IO.userError "invalid comparison mode"
    let rows ← IO.ofExcept (← readJson candidate).getArr?
    if rows.size != expected.size then throw <| IO.userError "row count mismatch"
    let mut differences := #[]
    for (row, (s, v)) in rows.zip expected do
      let actual ← IO.ofExcept (field row "values")
      let id ← IO.ofExcept ((field row "smiles") >>= Json.getStr?)
      if id != s then throw <| IO.userError "input SMILES or order changed"
      if v == Json.null then
        if actual != Json.null then throw <| IO.userError s!"{s}: expected unsupported molecule"
      else
        let xs ← IO.ofExcept actual.getArr?
        let ys ← IO.ofExcept v.getArr?
        if mode == "check" then
          if let .error e := verify cols xs ys then throw <| IO.userError s!"{s}: {e}"
        else
          unless xs.size == cols.size && ys.size == cols.size do throw <| IO.userError "profile length mismatch"
          let xn ← IO.ofExcept (xs.mapM number?)
          let yn ← IO.ofExcept (ys.mapM number?)
          for h : i in [:cols.size] do
            unless atColumn cols xn yn i do
              differences := differences.push <| Json.mkObj [
                ("smiles", toJson s), ("id", toJson cols[i].id), ("observed", xs[i]!), ("modeled", ys[i]!)]
    if mode == "audit" then IO.println (Json.arr differences).compress
    else IO.println s!"Verified {expected.size} rows against {cols.size} descriptor columns."
  | _, _ => throw <| IO.userError "invalid command or argument count"

def main (args : List String) : IO Unit := run args
