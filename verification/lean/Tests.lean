import Descriptors
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
