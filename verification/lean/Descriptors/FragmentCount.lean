import Descriptors.Fragment

namespace Descriptors.Fragment

private theorem count_le_one [DecidableEq α] (xs : List α) (h : xs.Nodup) (v : α) :
    xs.count v ≤ 1 := by
  induction xs with
  | nil => simp
  | cons a as ih =>
    obtain ⟨fresh, rest⟩ := List.nodup_cons.mp h
    by_cases eq : a = v
    · subst a
      simp [List.count_cons, List.count_eq_zero_of_not_mem fresh]
    · simpa [List.count_cons, eq, beq_iff_eq] using ih rest

private theorem perm_of_same_members [DecidableEq α] (xs ys : List α)
    (hx : xs.Nodup) (hy : ys.Nodup) (same : ∀ x, x ∈ xs ↔ x ∈ ys) : xs.Perm ys := by
  apply List.perm_iff_count.mpr
  intro key
  by_cases h : key ∈ ys
  · have left := List.count_pos_iff.mpr ((same key).mpr h)
    have right := List.count_pos_iff.mpr h
    have lb := count_le_one _ hx key
    have rb := count_le_one _ hy key
    omega
  · rw [List.count_eq_zero_of_not_mem h, List.count_eq_zero_of_not_mem (fun h' => h ((same key).mp h'))]

/-- Any duplicate-free enumeration containing exactly the images of valid
SMARTS assignments yields the same capped descriptor, independently of traversal order. -/
theorem fragmentCount_spec (g : Graph) (query : List Step) (images : List (List Nat))
    (nodup : images.Nodup)
    (complete : ∀ key, key ∈ images ↔ ∃ vs, Fits g [] query vs ∧ vs.mergeSort = key) :
    fragmentCount g query = min 1000 images.length := by
  obtain ⟨nd, spec⟩ := uniqueImages_spec g query
  have perm := perm_of_same_members _ _ nd nodup (fun key => (spec key).trans (complete key).symm)
  simp [fragmentCount, perm.length_eq]

end Descriptors.Fragment
