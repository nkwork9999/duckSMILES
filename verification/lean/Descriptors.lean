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

/-- Per-cell numerical contract, including the two undefined/defined cases. -/
def matchesValue (a r : JsonNumber) : Option JsonNumber → Option JsonNumber → Prop
  | none, none => True
  | some x, some y => within a r x y
  | _, _ => False

theorem agrees_spec (a r : JsonNumber) (x y : Option JsonNumber) :
    agrees a r x y = true ↔ matchesValue a r x y := by
  cases x <;> cases y <;> simp [matchesValue, agrees]

structure Column where
  id : String
  absolute : JsonNumber
  relative : JsonNumber
  deriving Inhabited

abbrev Values := Array (Option JsonNumber)

def atColumn (cols : Array Column) (xs ys : Values) (i : Nat) : Bool :=
  let c := cols[i]?.getD default
  agrees c.absolute c.relative (xs[i]?.join) (ys[i]?.join)

/-- A successful check certifies both lengths and the numerical contract at
every valid catalog index. The caller supplies the observed and modeled arrays. -/
def ProfileMatches (cols : Array Column) (xs ys : Values) : Prop :=
  xs.size = cols.size ∧ ys.size = cols.size ∧ ∀ (i : Fin cols.size),
    matchesValue cols[i].absolute cols[i].relative xs[i.val]?.join ys[i.val]?.join

theorem verify_profile_iff (cols : Array Column) (xs ys : Values) :
    checkAll cols.size xs.size ys.size (atColumn cols xs ys) = true ↔ ProfileMatches cols xs ys := by
  rw [checkAll_spec]
  unfold ProfileMatches
  constructor
  · rintro ⟨hx, hy, h⟩
    refine ⟨hx, hy, fun i => ?_⟩
    have hi := h i.val i.isLt
    simpa [atColumn, i.isLt, agrees_spec] using hi
  · rintro ⟨hx, hy, h⟩
    refine ⟨hx, hy, fun i hi => ?_⟩
    simpa [atColumn, hi, agrees_spec] using h ⟨i, hi⟩

/-- Counts in the supplied atom-label representation. Parsing and chemical
perception are outside the theorem. -/
def countAtoms (labels : List String) (symbol : String) : Nat := labels.count symbol

theorem countAtoms_perm {xs ys : List String} (p : xs.Perm ys) (s : String) :
    countAtoms xs s = countAtoms ys s := p.count_eq s

#print axioms agrees_defined
#print axioms agrees_presence
#print axioms checkAll_spec
#print axioms countAtoms_perm
/-- Partition an already perceived ring basis by a Boolean predicate. -/
def classCounts (rings : List α) (belongs : α → Bool) : Nat × Nat :=
  ((rings.filter belongs).length, (rings.filter (fun r => !(belongs r))).length)

theorem classCounts_partition (rings : List α) (p : α → Bool) :
    (classCounts rings p).1 + (classCounts rings p).2 = rings.length := by
  induction rings with
  | nil => rfl
  | cons r rs ih =>
    cases h : p r <;> simp [classCounts, List.filter, h] at ih ⊢ <;> omega

end Descriptors
