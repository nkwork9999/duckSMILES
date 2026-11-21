import Descriptors.Numerical

namespace Descriptors.Atomic
open Lean

/-- Binary64 contribution bit patterns from the prepared molecular observations. -/
structure Contribution where
  area : Nat
  logp : Nat
  mr : Nat
  estate : Nat
  partialCharge : Nat
  deriving FromJson

def decode (bits : Nat) : Float := Float.ofBits bits.toUInt64

def bin (edges : List Float) (x : Float) : Nat :=
  (edges.takeWhile fun edge => !(x < edge)).length

/-- The first strictly larger boundary selects the bin. Values exactly at an
edge continue to the following bin. Unordered comparisons continue too. -/
theorem bin_at_first_larger (before after : List Float) (edge x : Float)
    (passed : ∀ v ∈ before, ¬ x < v) (stop : x < edge) :
    bin (before ++ edge :: after) x = before.length := by
  induction before with
  | nil => simp [bin, List.takeWhile, stop]
  | cons y ys ih =>
    have hy := passed y (by simp)
    have ht := ih (fun v hv => passed v (by simp [hv]))
    simp [bin, List.takeWhile, hy] at ht ⊢
    exact ht

def orderedSum (xs : List Float) := xs.foldl (· + ·) 0.0

def family (label : String) (edges : List Float) (pairs : List (Float × Float)) :
    List (String × Json) :=
  (List.range (edges.length + 1)).map fun k =>
    let weights := (pairs.filter fun p => bin edges p.1 == k).map Prod.snd
    (label ++ toString (k+1), floatJson (orderedSum weights))

private def pyMin (a b : Float) := if b < a then b else a
private def pyMax (a b : Float) := if b > a then b else a

def values (rows : List Contribution) (hydrogenArea : Nat) : List (String × Json) := Id.run do
  let area := rows.map fun a => decode a.area
  let estate := rows.map fun a => decode a.estate
  let q := rows.map fun a => decode a.partialCharge
  let logp := rows.map fun a => decode a.logp
  let mr := rows.map fun a => decode a.mr
  let minQ := q.foldl (fun acc x => pyMin x acc) 500.0
  let maxQ := q.foldl (fun acc x => pyMax x acc) (-500.0)
  let extrema := fun (xs : List Float) (f : Float → Float → Float) =>
    match xs with
    | [] => Json.null
    | x :: tail => floatJson (tail.foldl f x)
  let absEstate := estate.map Float.abs
  return [
    ("LabuteASA", floatJson (orderedSum area + decode hydrogenArea)),
    ("MinPartialCharge", floatJson minQ), ("MaxPartialCharge", floatJson maxQ),
    ("MinAbsPartialCharge", floatJson (pyMin minQ.abs maxQ.abs)),
    ("MaxAbsPartialCharge", floatJson (pyMax minQ.abs maxQ.abs)),
    ("MinEStateIndex", extrema estate pyMin), ("MaxEStateIndex", extrema estate pyMax),
    ("MinAbsEStateIndex", extrema absEstate pyMin), ("MaxAbsEStateIndex", extrema absEstate pyMax)] ++
    family "SMR_VSA" [1.29,1.82,2.24,2.45,2.75,3.05,3.63,3.8,4.0] (mr.zip area) ++
    family "SlogP_VSA" [-0.4,-0.2,0.0,0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.6] (logp.zip area) ++
    family "EState_VSA" [-0.39,0.29,0.717,1.165,1.54,1.807,2.05,4.69,9.17,15.0] (estate.zip area) ++
    family "VSA_EState" [4.78,5.0,5.41,5.74,6.0,6.07,6.45,7.0,11.0] (area.zip estate) ++
    family "PEOE_VSA" [-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,0.0,0.05,0.1,0.15,0.2,0.25,0.3] (q.zip area)

end Descriptors.Atomic
