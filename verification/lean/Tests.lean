import Descriptors
import Descriptors.FragmentCount
import Descriptors.FragmentPatterns
import Descriptors.Atomic
import Descriptors.GraphIndices
import Descriptors.Information
open Descriptors Lean

example : agrees 0 0 none (some 0) = false := by decide
example : agrees 0 0 (some 0) none = false := by decide
example : agrees 0 0 none none = true := by decide
example : agrees 0 0 (some ⟨10, 1⟩) (some 1) = true := by decide
example : agrees 0 0 (some (-1)) (some 1) = false := by decide
example : agrees ⟨1, 2⟩ 0 (some ⟨101, 2⟩) (some 1) = true := by decide
example : agrees ⟨1, 2⟩ 0 (some ⟨1011, 3⟩) (some 1) = false := by decide
example : agrees 0 ⟨1, 2⟩ (some (-101)) (some (-100)) = true := by decide
example : agrees 0 ⟨1, 2⟩ (some (-101)) (some (-99)) = false := by decide

open Descriptors.Fragment

-- The two orientations of C-C use the same two atoms.
private def ethane : Graph := Graph.ofBonds #[{ z := 6 }, { z := 6 }] [(0, 1, 1)]
private def carbonLink : List Step := [
  ⟨fun i => (ethane.atom i).z == 6, []⟩,
  ⟨fun i => (ethane.atom i).z == 6, [(0, fun b => b == 1)]⟩]
example : fragmentCount ethane carbonLink = 1 := by native_decide
example : fragmentCount (Graph.ofBonds #[{z := 6}] []) carbonLink = 0 := by native_decide

-- A repeated orientation has one image. A second disconnected bond has another.
private def twoBonds : Graph := Graph.ofBonds (Array.replicate 4 {z := 6}) [(0,1,1),(2,3,1)]
-- With a query built from the actual target both disconnected components match.
private def bondIn (g : Graph) : List Step :=
  [⟨fun i => (g.atom i).z == 6, []⟩, ⟨fun i => (g.atom i).z == 6, [(0, fun b => b == 1)]⟩]
example : fragmentCount twoBonds (bondIn twoBonds) = 2 := by native_decide

example : Atomic.bin [0.0, 1.0] 0.0 = 1 := by native_decide
example : Atomic.bin [0.0, 1.0] (0.0 / 0.0) = 2 := by native_decide
example : Information.divideExact 7 2 = none := by decide
example : Information.divideExact (-6) 2 = some (-3) := by decide
example : GraphIndices.walkRow 2 (fun i => if i == 0 then [1] else [0]) 0 3 = #[0,1] := by native_decide

-- The order-80 constant term has the opposite behavior to the rounded
-- RDKit matrix recurrence on this control. Compute the matrix coefficients
-- for the actual molecular graph as well as the analytic path recurrence.
private def carbonPath (n : Nat) : Graph :=
  Graph.ofBonds (Array.replicate n {z := 6}) ((List.range (n-1)).map fun i => (i,i+1,1))

example : Information.pathCoefficient 80 0 = 1 := by
  rw [show 80 = 2*40 by decide, Information.path_even_constant]
  decide

#guard match Information.coefficients (carbonPath 80) with
  | .ok coefficients => coefficients.getLast? == some 1
  | .error _ => false

-- 0–1–2 and one isolated vertex. Reject a false shortcut, false
-- unreachability, and a finite value for the isolated vertex.
private def distanceChain : GraphIndices.Adj
  | 0 => [1]
  | 1 => [0,2]
  | 2 => [1]
  | _ => []

private def verifyDistance (row : List Nat) :=
  GraphIndices.Distance.accepts 4 distanceChain 0 (fun v => row[v]!)

example : verifyDistance [0,1,2,5] = true := by decide
example : verifyDistance [0,1,1,5] = false := by decide
example : verifyDistance [0,1,5,5] = false := by decide
example : verifyDistance [0,1,2,3] = false := by decide
