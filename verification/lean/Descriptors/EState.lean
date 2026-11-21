import Descriptors.GraphIndices
import Descriptors.Composition

namespace Descriptors.EState

/-- Hall–Kier intrinsic state before the distance-dependent perturbation. -/
def intrinsic (z degree outer attachedH : Nat) : Float :=
  if degree == 0 then 0 else
    let period : Nat :=
      if z ≤ 2 then 1 else if z ≤ 10 then 2 else if z ≤ 18 then 3 else
      if z ≤ 36 then 4 else if z ≤ 54 then 5 else if z ≤ 86 then 6 else 7
    let valence := if attachedH > outer then -(attachedH - outer).toFloat else (outer - attachedH).toFloat
    ((2.0 / period.toFloat) * (2.0 / period.toFloat) * valence + 1) / degree.toFloat

def indices (g : Fragment.Graph) (props : Array Composition.AtomProperties) :
    Except String (List Float) := do
  let n := g.atoms.size
  let vs := List.range n
  let adj := fun i => (g.neighbors[i]?.getD []).map Prod.fst
  let d := GraphIndices.distances n vs adj
  let base := Array.ofFn (n := n) fun i =>
    let a := g.atom i.val
    let p := props[i.val]!
    intrinsic a.z a.degree p.outer p.attachedH
  let mut shift := Array.replicate n (0 : Float)
  for i in vs do
    unless GraphIndices.Distance.accepts n adj i (fun j => d[i]![j]!) do
      throw "invalid EState distance row"
    for j in [i+1:n] do
      let length := d[i]![j]!
      if length ≤ n then
        let r := (length + 1).toFloat
        let delta := (base[i]! - base[j]!) / (r*r)
        shift := (shift.modify i (· + delta)).modify j (· - delta)
  return (base.toList.zipWith (· + ·) shift.toList)

end Descriptors.EState
