import Descriptors.Composition

namespace Descriptors.Morgan
open Lean Fragment Atomic Composition

private def combine (seed value : UInt32) : UInt32 :=
  seed ^^^ (value + (0x9e3779b9 : UInt32) + (seed <<< 6) + (seed >>> 2))

private def encodeSigned (x : Int) : UInt32 :=
  (x % ((2 : Int)^32)).toNat.toUInt32

def connectivityInvariant (atom : Atom) (prop : AtomProperties) : UInt32 :=
  let isotopeDelta := ((decode prop.mass - decode prop.averageMass).toInt64).toInt
  let terms := [atom.z.toUInt32, atom.totalDegree.toUInt32, atom.hydrogens.toUInt32,
                encodeSigned atom.charge, encodeSigned isotopeDelta]
  (if atom.rings > 0 then terms ++ [1] else terms).foldl combine 0

/-- The bond-set mask is an unbounded natural, so ordering compares the most
significant differing bond first regardless of the number of machine words. -/
def values (g : Graph) (edges : List (Nat × Nat × Nat)) (props : Array AtomProperties) :
    List (String × Json) := Id.run do
  let n := g.atoms.size
  let heavy := (g.atoms.toList.filter (fun a => a.z > 1)).length
  let initial := g.atoms.mapIdx (fun i a => connectivityInvariant a props[i]!)
  let mut current := initial
  let mut previous := Array.replicate n (0 : Nat)
  let mut dead := Array.replicate n false
  let mut emitted : Std.HashSet Nat := {}
  let mut unique := Std.HashSet.ofList initial.toList
  let mut out := []
  let indexed := edges.zipIdx.foldl (init := Array.replicate n []) fun rows ((u,v,b),i) =>
    (rows.modify u ((v,b,i) :: ·)).modify v ((u,b,i) :: ·)
  for layer in [:3] do
    let mut next := Array.replicate n (0 : UInt32)
    let mut round := previous
    let mut candidates : List (Nat × UInt32 × Nat) := []
    for i in [:n] do
      let adj := indexed[i]!
      if !dead[i]! && !adj.isEmpty then
        let mut mask := previous[i]!
        for (v,_,index) in adj do mask := mask ||| (1 <<< index) ||| previous[v]!
        round := round.set! i mask
        let labels := ((adj.map fun (v,b,_) => (b.toUInt32,current[v]!)).mergeSort (fun a b => compare a b != .gt))
        let key := labels.foldl (fun seed (order,other) => combine seed (combine (combine 0 order) other))
            (combine layer.toUInt32 current[i]!)
        next := next.set! i key
        candidates := (mask,key,i) :: candidates
    for (mask,key,i) in (candidates.mergeSort (fun a b => compare a b != .gt)) do
      if emitted.contains mask then dead := dead.set! i true
      else
        emitted := emitted.insert mask
        unique := unique.insert key
    current := next
    previous := round
    let density := if heavy == 0 then 0.0 else unique.size.toFloat / heavy.toFloat
    out := (s!"FpDensityMorgan{layer+1}", floatJson density) :: out
  return out

end Descriptors.Morgan
