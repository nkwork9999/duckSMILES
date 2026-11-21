import Lean

namespace Descriptors.GraphIndices.Distance

/-- A walk of exactly the given length, extended one edge at a time. -/
inductive Walk (adj : Nat → List Nat) (s : Nat) : Nat → Nat → Prop
  | zero : Walk adj s 0 s
  | next {k u v} : Walk adj s k u → v ∈ adj u → Walk adj s (k + 1) v

/-- Values above n encode unreachable vertices. Every finite vertex other
than the source must have an incoming edge from one level closer. -/
def Valid (n : Nat) (adj : Nat → List Nat) (s : Nat) (d : Nat → Nat) : Prop :=
  s < n ∧ d s = 0 ∧
  ∀ v < n, (d v = 0 → v = s) ∧
    (d v ≤ n → ∀ w ∈ adj v, w < n ∧ d w ≤ n ∧ d w ≤ d v + 1) ∧
    (0 < d v → d v ≤ n → ∃ u, u < n ∧ u ∈ adj v ∧ v ∈ adj u ∧ d u + 1 = d v)

def accepts (n : Nat) (adj : Nat → List Nat) (s : Nat) (d : Nat → Nat) : Bool :=
  decide (s < n ∧ d s = 0) &&
  (List.range n).all fun v =>
    decide (d v ≠ 0 ∨ v = s) &&
    (decide (n < d v) || (adj v).all fun w =>
      decide (w < n ∧ d w ≤ n ∧ d w ≤ d v + 1)) &&
    (decide (d v = 0 ∨ n < d v) || (adj v).any fun u =>
      decide (u < n ∧ v ∈ adj u ∧ d u + 1 = d v))

theorem accepts_iff (n : Nat) (adj : Nat → List Nat) (s : Nat) (d : Nat → Nat) :
    accepts n adj s d = true ↔ Valid n adj s d := by
  simp only [accepts, Valid, Bool.and_eq_true, Bool.or_eq_true, decide_eq_true_eq,
    List.all_eq_true, List.any_eq_true, List.mem_range]
  constructor
  · rintro ⟨⟨hs, hz⟩, h⟩
    refine ⟨hs,hz,fun v hv => ?_⟩
    obtain ⟨⟨zero, edge⟩, pred⟩ := h v hv
    refine ⟨?_, ?_, ?_⟩
    · intro eq
      exact zero.resolve_left (fun h => h eq)
    · intro hvn
      exact edge.resolve_left (by omega)
    · intro hp hn
      obtain ⟨u, mem, lt, reverse, eq⟩ := pred.resolve_left (by omega)
      exact ⟨u, lt, mem, reverse, eq⟩
  · rintro ⟨hs, hz, h⟩
    refine ⟨⟨hs,hz⟩,fun v hv => ?_⟩
    obtain ⟨zero,edge,pred⟩ := h v hv
    refine ⟨⟨?_, ?_⟩, ?_⟩
    · by_cases eq : d v = 0
      · exact Or.inr (zero eq)
      · exact Or.inl eq
    · by_cases hn : d v ≤ n
      · exact Or.inr (edge hn)
      · exact Or.inl (by omega)
    · by_cases hn : d v ≤ n
      · by_cases eq : d v = 0
        · exact Or.inl (Or.inl eq)
        · obtain ⟨u, lt, mem, reverse, step⟩ := pred (by omega) hn
          exact Or.inr ⟨u, mem, lt, reverse, step⟩
      · exact Or.inl (Or.inr (by omega))

theorem walk_lower_bound {n adj s d : _} (valid : Valid n adj s d)
    {k v : Nat} (path : Walk adj s k v) : v < n ∧ d v ≤ n ∧ d v ≤ k := by
  induction path with
  | zero =>
    obtain ⟨hs, hz, _⟩ := valid
    exact ⟨hs, by omega, by omega⟩
  | next _ edge ih =>
    obtain ⟨hv, finite, bound⟩ := ih
    obtain ⟨hw, hwfin, step⟩ := (valid.2.2 _ hv).2.1 finite _ edge
    exact ⟨hw, hwfin, by omega⟩

theorem finite_witness {n adj s d : _} (valid : Valid n adj s d)
    {v : Nat} (hv : v < n) (hf : d v ≤ n) :
    Walk adj s (d v) v := by
  generalize eq : d v = k
  induction k using Nat.strongRecOn generalizing v with
  | ind k ih =>
    by_cases hz : k = 0
    · subst k
      have zero : d v = 0 := by omega
      have same := (valid.2.2 v hv).1 zero
      subst v
      simpa [valid.2.1] using (Walk.zero (adj := adj) (s := s))
    · obtain ⟨u, hu, _, edge, step⟩ := (valid.2.2 v hv).2.2 (by omega) hf
      have path := ih (d u) (by omega) hu (by omega) rfl
      have extended := Walk.next path edge
      rwa [step, eq] at extended

/-- Successful validation is sufficient for exact distances and for
unreachability. The proof applies to a row computed by any implementation. -/
theorem accepts_shortest {n adj s d : _} (ok : accepts n adj s d = true)
    {v : Nat} (hv : v < n) :
    (d v ≤ n ↔ ∃ k, Walk adj s k v) ∧
    (d v ≤ n → Walk adj s (d v) v ∧ ∀ k, Walk adj s k v → d v ≤ k) := by
  have valid := (accepts_iff ..).mp ok
  refine ⟨⟨fun hf => ⟨d v, finite_witness valid hv hf⟩,
    fun ⟨_, p⟩ => (walk_lower_bound valid p).2.1⟩, ?_⟩
  exact fun hf => ⟨finite_witness valid hv hf, fun _ p => (walk_lower_bound valid p).2.2⟩

#print axioms accepts_shortest

end Descriptors.GraphIndices.Distance

namespace Descriptors.GraphIndices.Distance

/-- Any two accepted rows for the same graph and source agree on reachability
and on every finite distance. Unreachable sentinels may differ. -/
theorem verified_rows_agree {n adj s a b : _}
    (ha : accepts n adj s a = true) (hb : accepts n adj s b = true)
    {v : Nat} (hv : v < n) :
    (a v ≤ n ↔ b v ≤ n) ∧ (a v ≤ n → a v = b v) := by
  obtain ⟨ar, ashort⟩ := accepts_shortest ha hv
  obtain ⟨br, bshort⟩ := accepts_shortest hb hv
  have reach := ar.trans br.symm
  refine ⟨reach, fun finite => Nat.le_antisymm ?_ ?_⟩
  · exact (ashort finite).2 _ (bshort (reach.mp finite)).1
  · exact (bshort (reach.mp finite)).2 _ (ashort finite).1

#print axioms verified_rows_agree

end Descriptors.GraphIndices.Distance
