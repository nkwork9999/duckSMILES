import Descriptors.Fragment

namespace Descriptors.Labute
open Fragment

-- Rb0 radii, RDKit 2026.03.6. See RDKIT-DATA-LICENSE.txt.
private def radii : Array Float := #[
    0.33, 0.7, 1.23, 0.9, 0.82, 0.77, 0.7, 0.66, 0.611, 0.7, 1.54, 1.36, 1.18, 0.937, 0.89, 1.04,
    0.997, 1.74, 2.03, 1.74, 1.44, 1.32, 1.22, 1.18, 1.17, 1.17, 1.16, 1.15, 1.17, 1.25, 1.26,
    1.188, 1.2, 1.17, 1.167, 1.91, 2.16, 1.91, 1.62, 1.45, 1.34, 1.3, 1.27, 1.25, 1.25, 1.28, 1.34,
    1.48, 1.44, 1.385, 1.4, 1.378, 1.387, 1.98, 2.35, 1.98, 1.69, 1.83, 1.82, 1.81, 1.8, 1.8, 1.99,
    1.79, 1.76, 1.75, 1.74, 1.73, 1.72, 1.94, 1.72, 1.44, 1.34, 1.3, 1.28, 1.26, 1.27, 1.3, 1.34,
    1.49, 1.48, 1.48, 1.45, 1.46, 1.45, 2.4, 2.0, 1.9, 1.88, 1.79, 1.61, 1.58, 1.55, 1.53, 1.07,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0,
]

structure Areas where
  atoms : List Float
  hydrogen : Float

def overlap (a b correction : Float) : Float :=
  let r := min (a+b) (max (a-b).abs (a+b-correction))
  b*b - (a-r)*(a-r)/r

def areas (g : Graph) (edges : List (Nat × Nat × Nat)) : Except String Areas := do
  let n := g.atoms.size
  unless g.atoms.all (fun a => a.z > 0 && a.z <= radii.size) do
    throw "Labute radius unavailable for atom"
  let rs := g.atoms.map (fun a => radii[a.z-1]!)
  let mut sums := Array.replicate n (0 : Float)
  for (i,j,order) in edges do
    let correction := match order with
      | 2 => 0.2 | 3 => 0.3 | 12 => 0.1 | _ => 0
    sums := sums.modify i (· + overlap rs[i]! rs[j]! correction)
    sums := sums.modify j (· + overlap rs[j]! rs[i]! correction)
  let h := radii[0]!
  let pi : Float := 3.141592653589793
  let fromOverlap (r v : Float) := pi*r*(4*r-v)
  let hydrogen := rs.foldl (fun v r => v + overlap h r 0) 0.0
  let atomAreas := rs.mapIdx fun i r => fromOverlap r (sums[i]! + overlap r h 0)
  return ⟨atomAreas.toList, if hydrogen.abs > 1e-4 then fromOverlap h hydrogen else hydrogen⟩

end Descriptors.Labute
