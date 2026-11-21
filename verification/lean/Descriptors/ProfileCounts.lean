import Descriptors

namespace Descriptors.ProfileCounts
open Lean

structure Input where
  symbols : List String
  charges : List Int
  explicitH : List Nat
  implicitH : List Nat
  bondOrders : List String
  aromaticAtoms : List Bool
  ringAtoms : List Bool
  ringBonds : List Bool
  componentSizes : List Nat
  ringSizes : List Nat

def countedElements : List String := ["C", "N", "O", "S", "P", "F", "Cl", "Br", "I"]

def maximum (xs : List Nat) : Nat := xs.foldl max 0

theorem maximum_le (xs : List Nat) (n : Nat) :
    maximum xs ≤ n ↔ ∀ k ∈ xs, k ≤ n := by
  have fold (xs : List Nat) (initial : Nat) :
      xs.foldl max initial ≤ n ↔ initial ≤ n ∧ ∀ k ∈ xs, k ≤ n := by
    induction xs generalizing initial with
    | nil => simp
    | cons x xs ih => simp [List.foldl, ih, Nat.max_le, and_assoc]
  simp [maximum, fold]

/-- Complete evaluator for the 24 native composition and perceived-graph columns. -/
def values (x : Input) : List (String × Json) :=
  [("ds_BondCount", toJson x.bondOrders.length),
   ("ds_FormalCharge", toJson x.charges.sum),
   ("ds_FragmentCount", toJson x.componentSizes.length),
   ("ds_ExplicitHydrogens", toJson x.explicitH.sum),
   ("ds_ImplicitHydrogens", toJson x.implicitH.sum),
   ("ds_TotalHydrogens", toJson (x.explicitH.sum + x.implicitH.sum)),
   ("ds_HalogenCount", toJson ((x.symbols.filter fun e => ["F", "Cl", "Br", "I"].contains e).length)),
   ("ds_SingleBonds", toJson (x.bondOrders.count "SINGLE")),
   ("ds_DoubleBonds", toJson (x.bondOrders.count "DOUBLE")),
   ("ds_TripleBonds", toJson (x.bondOrders.count "TRIPLE")),
   ("ds_AromaticBonds", toJson (x.bondOrders.count "AROMATIC")),
   ("ds_AromaticAtoms", toJson (x.aromaticAtoms.count true)),
   ("ds_RingAtoms", toJson (x.ringAtoms.count true)),
   ("ds_RingBonds", toJson (x.ringBonds.count true)),
   ("ds_LargestRingSize", toJson (maximum x.ringSizes))] ++
  countedElements.map fun e => ("ds_Count" ++ e, toJson (countAtoms x.symbols e))

/-- The input properties can be enumerated in different orders. This
relation concerns already perceived graphs and their attached attributes. -/
structure Reordered (x y : Input) : Prop where
  symbols : x.symbols.Perm y.symbols
  charges : x.charges.Perm y.charges
  explicitH : x.explicitH.Perm y.explicitH
  implicitH : x.implicitH.Perm y.implicitH
  bondOrders : x.bondOrders.Perm y.bondOrders
  aromaticAtoms : x.aromaticAtoms.Perm y.aromaticAtoms
  ringAtoms : x.ringAtoms.Perm y.ringAtoms
  ringBonds : x.ringBonds.Perm y.ringBonds
  componentSizes : x.componentSizes.Perm y.componentSizes
  ringSizes : x.ringSizes.Perm y.ringSizes

private theorem natSum {xs ys : List Nat} (h : xs.Perm ys) : xs.sum = ys.sum := by
  induction h with
  | nil => rfl
  | cons _ _ ih => simp [ih]
  | swap => simp [Nat.add_left_comm]
  | trans _ _ ha hb => exact ha.trans hb

private theorem intSum {xs ys : List Int} (h : xs.Perm ys) : xs.sum = ys.sum := by
  induction h with
  | nil => rfl
  | cons _ _ ih => simp [ih]
  | swap => simp [Int.add_left_comm]
  | trans _ _ ha hb => exact ha.trans hb

theorem values_reordered {x y : Input} (h : Reordered x y) : values x = values y := by
  have ringMax : maximum x.ringSizes = maximum y.ringSizes := by
    apply Nat.le_antisymm
    · exact (maximum_le ..).mpr fun k hk =>
        (maximum_le ..).mp (Nat.le_refl _) k (h.ringSizes.mem_iff.mp hk)
    · exact (maximum_le ..).mpr fun k hk =>
        (maximum_le ..).mp (Nat.le_refl _) k (h.ringSizes.mem_iff.mpr hk)
  have elements : ∀ e, countAtoms x.symbols e = countAtoms y.symbols e :=
    countAtoms_perm h.symbols
  have halogen := (h.symbols.filter fun e => ["F", "Cl", "Br", "I"].contains e).length_eq
  unfold values
  rw [halogen, intSum h.charges, natSum h.explicitH, natSum h.implicitH, ringMax]
  simp [h.bondOrders.length_eq, h.componentSizes.length_eq, h.bondOrders.count_eq,
    h.aromaticAtoms.count_eq, h.ringAtoms.count_eq, h.ringBonds.count_eq, elements]

#print axioms values_reordered

end Descriptors.ProfileCounts
