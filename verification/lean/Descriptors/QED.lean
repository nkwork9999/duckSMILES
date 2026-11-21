import Descriptors.DrugLikeness
import Descriptors.FragmentPatterns
import Descriptors.GraphIndices
import Descriptors.QEDPatterns

namespace Descriptors.QED
open Lean Fragment Atomic

structure AtomInput where
  mass : Nat
  attachedH : Nat
  tpsa : Nat
  deriving FromJson

structure Input where
  vertices : Array Atom
  edges : List (List Nat)
  ringFlags : List Bool
  atoms : List AtomInput
  expandedLogP : List Nat
  deriving FromJson

/-- Cycle rank of the graph left by RDKit's aliphatic-ring deletion query.
The ring rank excludes dative and hydrogen bonds, as does GetSSSR by default. -/
def aromaticRank (g : Graph) (edges : List (Nat × Nat × Nat)) : Except String Nat := do
  let deleted := fun i => (qedAliphatic g).any (rooted g i)
  let vs := (List.range g.atoms.size).filter (fun i => !deleted i)
  let keep := edges.filter fun (i,j,order) =>
    vs.contains i && vs.contains j && order != 17 && order != 14
  let n := g.atoms.size
  let reduced := Graph.ofBonds g.atoms keep
  let adj := fun i => (reduced.neighbors[i]?.getD []).map Prod.fst
  let d := GraphIndices.distances n vs adj
  let mut components := 0
  for i in vs do
    let row := fun j => d[i]![j]!
    unless GraphIndices.Distance.accepts n adj i row do
      throw "invalid QED ring distance row"
    if vs.all (fun j => j >= i || row j > n) then components := components + 1
  return keep.length + components - vs.length

def components (x : Input) : Except String (List Float) := do
  unless x.atoms.length == x.vertices.size && x.edges.length == x.ringFlags.length do
    throw "QED graph dimensions disagree with prepared properties"
  let edges ← x.edges.mapM fun e => do
    let [i,j,bond] := e | throw "invalid QED bond record"
    unless i < x.vertices.size && j < x.vertices.size && i != j do
      throw "invalid QED bond endpoints"
    pure (i,j,bond)
  let rings := (edges.zip x.ringFlags).filterMap fun ((i,j,_), ring) =>
    if ring then some (min i j, max i j) else none
  let g := Graph.ofBonds x.vertices edges rings
  let atomicMass := orderedSum (x.atoms.map (fun a => decode a.mass))
  let hydrogens := (x.atoms.map AtomInput.attachedH).sum
  let mw := atomicMass + hydrogens.toFloat * 1.008
  let lp := orderedSum (x.expandedLogP.map decode)
  let hba := ((qedAcceptors g).map (fragmentCount g)).sum
  let some donors := pattern? g "NumHDonors" | throw "missing registered donor pattern"
  let hbd := fragmentCount g donors
  let psa := orderedSum (x.atoms.map (fun a => decode a.tpsa))
  let some rotors := pattern? g "NumRotatableBonds" | throw "missing registered rotor pattern"
  let rotb := fragmentCount g rotors
  let arom ← aromaticRank g edges
  let alerts := ((qedAlerts g).filter (hasEmbedding g [])).length
  return [mw, lp, hba.toFloat, hbd.toFloat, psa, rotb.toFloat, arom.toFloat, alerts.toFloat]

def value (x : Input) : Except String Json := do
  DrugLikeness.value (← components x)

end Descriptors.QED
