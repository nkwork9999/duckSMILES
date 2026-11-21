import Descriptors.Fragment
import Descriptors.Numerical

namespace Descriptors.Information
open Lean

/-- Check each exact trace division in the Faddeev-LeVerrier recurrence. -/
def divideExact (t : Int) (k : Nat) : Option Int :=
  let q := t / (k : Int)
  if k > 0 && q * k == t then some q else none

theorem divideExact_correct {t q : Int} {k : Nat} (h : divideExact t k = some q) :
    k > 0 ∧ q * k = t := by
  simp only [divideExact] at h
  split at h
  next hp => cases h; simpa using hp
  next => contradiction

def coefficients (g : Fragment.Graph) : Except String (List Int) := do
  let n := g.atoms.size
  let mut previous := Array.ofFn (n := n) fun i => Array.ofFn (n := n) fun j =>
    if i == j then (1 : Int) else 0
  let mut cs : List Int := [1]
  for k in [1:n+1] do
    let mut product := Array.ofFn (n := n) fun i => Array.ofFn (n := n) fun j =>
      (g.neighbors[i.val]!).foldl (fun t (v, _) => t + previous[v]![j.val]!) (0 : Int)
    let trace := (List.range n).foldl (fun t i => t + product[i]![i]!) (0 : Int)
    let some q := divideExact (-trace) k | throw "inexact characteristic-polynomial coefficient"
    cs := q :: cs
    for i in [:n] do product := product.modify i (fun row => row.modify i (· + q))
    previous := product
  unless previous.all (fun row => row.all (· == 0)) do throw "nonzero final characteristic-polynomial residual"
  return cs.reverse

def values (g : Fragment.Graph) : Except String (List (String × Json)) := do
  let weights := ((← coefficients g).map Int.natAbs).filter (· != 0)
  let total := weights.sum
  let avg := -weights.foldl (init := 0.0) fun acc value =>
    let p := value.toFloat / total.toFloat
    acc + p * Float.log2 p
  return [("Ipc", floatJson (total.toFloat * avg)), ("AvgIpc", floatJson avg)]

end Descriptors.Information
