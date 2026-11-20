import Lean

namespace Descriptors
open Lean

/-- Bring every decimal in the error inequality to the same integer scale.
The relative term is relative to the reference, including for negative values. -/
def within (absolute relative actual expected : JsonNumber) : Prop :=
  let e := max (max actual.exponent expected.exponent)
               (max absolute.exponent (relative.exponent + expected.exponent))
  let scaled := fun (n : JsonNumber) => n.mantissa * (10 : Int) ^ (e - n.exponent)
  let product : JsonNumber :=
    ⟨relative.mantissa * expected.mantissa.natAbs, relative.exponent + expected.exponent⟩
  (scaled actual - scaled expected).natAbs ≤ max (scaled absolute).toNat (scaled product).toNat

instance (a r x y : JsonNumber) : Decidable (within a r x y) := by
  unfold within
  infer_instance

def agrees (absolute relative : JsonNumber) : Option JsonNumber → Option JsonNumber → Bool
  | none, none => true
  | some actual, some expected => decide (within absolute relative actual expected)
  | _, _ => false

theorem agrees_defined (a r x y : JsonNumber) :
    agrees a r (some x) (some y) = true ↔ within a r x y := by simp [agrees]

theorem agrees_presence (a r : JsonNumber) (x y : Option JsonNumber)
    (h : agrees a r x y = true) : x.isSome = y.isSome := by
  cases x <;> cases y <;> simp_all [agrees]

/-- The complete catalog must match. A short prefix, long suffix, or any
failing coordinate makes the profile fail. -/
def checkAll (n : Nat) (actualSize expectedSize : Nat) (check : Nat → Bool) : Bool :=
  actualSize == n && expectedSize == n && (List.range n).all check

theorem checkAll_spec {n a e : Nat} {check : Nat → Bool} :
    checkAll n a e check = true ↔ a = n ∧ e = n ∧ ∀ i < n, check i = true := by
  simp [checkAll, and_assoc]

/-- Counts in the supplied atom-label representation. Parsing and chemical
perception are outside the theorem. -/
def countAtoms (labels : List String) (symbol : String) : Nat := labels.count symbol

theorem countAtoms_perm {xs ys : List String} (p : xs.Perm ys) (s : String) :
    countAtoms xs s = countAtoms ys s := p.count_eq s

#print axioms agrees_defined
#print axioms agrees_presence
#print axioms checkAll_spec
#print axioms countAtoms_perm
end Descriptors
