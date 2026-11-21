import Descriptors.Atomic

namespace Descriptors.Burden
open Lean Atomic

structure Input where
  diagonal : Array (Array Nat)
  edges : List (List Nat)
  parameterized : Bool
  deriving FromJson

private def symmetricSet (a : Array Float) (n i j : Nat) (x : Float) :=
  (a.set! (i*n+j) x).set! (j*n+i) x

/-- Cyclic Jacobi rotations on a real symmetric Burden matrix. Termination
requires a small off-diagonal row sum, also bounding each Gershgorin radius. -/
def extrema (input : Array Float) (n : Nat) : Except String (Float × Float) := do
  if n == 0 then return (0,0)
  let mut a := input
  for _ in [:64] do
    let mut residual := 0.0
    for i in [:n] do
      let mut row := 0.0
      for j in [:n] do
        if i != j then row := row + a[i*n+j]!.abs
      residual := max residual row
    if residual <= 1.0e-12 then
      let diag := (List.range n).map fun i => a[i*n+i]!
      return (diag.foldl min a[0]!, diag.foldl max a[0]!)
    for p in [:n] do
      for q in [p+1:n] do
        let entry := a[p*n+q]!
        if entry.abs <= 1.0e-15 / n.toFloat then continue
        let app := a[p*n+p]!
        let aqq := a[q*n+q]!
        let tau := (aqq-app)/(2*entry)
        let t := (if tau < 0 then -1.0 else 1.0) / (tau.abs + (1 + tau*tau).sqrt)
        let c := 1.0 / (1+t*t).sqrt
        let s := t*c
        a := symmetricSet a n p q 0
        a := (a.set! (p*n+p) (app - t*entry)).set! (q*n+q) (aqq + t*entry)
        for k in [:n] do
          if k != p && k != q then
            let x := a[k*n+p]!
            let y := a[k*n+q]!
            a := symmetricSet (symmetricSet a n k p (c*x-s*y)) n k q (s*x+c*y)
  throw "Burden matrix eigensolver did not converge"

def values (obs : Input) : Except String (List (String × Json)) := do
  let n := obs.diagonal.size
  unless obs.diagonal.all (fun row => row.size == 4) do throw "invalid BCUT atomic properties"
  let mut bondsValid := true
  let mut base := Array.replicate (n*n) 0.001
  for edge in obs.edges do
    let [i,j,b] := edge | throw "invalid BCUT edge"
    unless i < n && j < n && i != j do throw "invalid BCUT endpoints"
    let val := match b with | 1 => 1.0 | 2 => 1.0/(2.0 : Float).sqrt | 3 => 1.0/(3.0 : Float).sqrt
                            | 12 => 1.0/(1.5 : Float).sqrt | _ => 0.0
    if val == 0 then bondsValid := false
    base := symmetricSet base n i j val
  let mut result := []
  for (key,col) in [("MW","LOW"),("CHG","LO"),("LOGP","LOW"),("MR","LOW")].zipIdx do
    let (stem, suffixLo) := key
    let diag := obs.diagonal.map fun row => decode row[col]!
    let (lo, hi) ← if obs.parameterized && bondsValid && diag.all Float.isFinite then do
        let mut matrix := base
        for i in [:n] do matrix := matrix.set! (i*n+i) diag[i]!
        let (lo, hi) ← extrema matrix n
        pure (floatJson lo, floatJson hi)
      else pure (Json.null, Json.null)
    result := ("BCUT2D_" ++ stem ++ suffixLo, lo) :: ("BCUT2D_" ++ stem ++ "HI", hi) :: result
  return result

end Descriptors.Burden
