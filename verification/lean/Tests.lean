import Descriptors
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
