import Lean

namespace Descriptors.Fragment

/-- Molecular properties used by the pinned SMARTS vocabulary. -/
structure Atom where
  z : Nat
  aromatic : Bool := false
  hydrogens : Nat := 0
  charge : Int := 0
  totalDegree : Nat := 0
  degree : Nat := 0
  rings : Nat := 0
  valence : Nat := 0
  deriving Inhabited, Lean.FromJson

structure Graph where
  atoms : Array Atom
  neighbors : Array (List (Nat × Nat))
  ringEdges : List (Nat × Nat) := []

def Graph.ofBonds (atoms : Array Atom) (bonds : List (Nat × Nat × Nat)) (ringEdges : List (Nat × Nat) := []) : Graph :=
  ⟨atoms, bonds.foldl (init := Array.replicate atoms.size []) fun rows (i, j, order) =>
    (rows.modify i ((j, order) :: ·)).modify j ((i, order) :: ·), ringEdges⟩

def Graph.atom (g : Graph) (v : Nat) := g.atoms[v]?.getD default

structure Step where
  accepts : Nat → Bool
  /-- Offset in the reversed prefix. Bond predicates use the low five bits for
  the RDKit order and bit five for ring membership. -/
  links : List (Nat × (Nat → Bool))

def bondCode (g : Graph) (u v order : Nat) :=
  order + (if g.ringEdges.contains (min u v, max u v) then 32 else 0)

def allowed (g : Graph) (s : Step) (used : List Nat) (v : Nat) : Bool :=
  s.accepts v && !used.contains v && s.links.all fun (idx, p) =>
    (used[idx]?).any fun u => ((g.neighbors[v]?).getD []).any fun (w, order) => w == u && p (bondCode g v w order)

/-- Full suffixes in query vertex order. A recursive SMARTS starts a fresh
assignment with its root constrained to the enclosing atom. -/
def embeddings (g : Graph) (used : List Nat) : List Step → List (List Nat)
  | [] => [[]]
  | s :: ss => ((List.range g.atoms.size).filter (allowed g s used)).flatMap fun v =>
      (embeddings g (v :: used) ss).map (v :: ·)

/-- Declarative witness for an injection satisfying all atom and bond tests. -/
inductive Fits (g : Graph) : List Nat → List Step → List Nat → Prop
  | nil (used) : Fits g used [] []
  | cons {used s ss v vs} : v < g.atoms.size → s.accepts v = true → v ∉ used →
      (∀ link ∈ s.links, ∃ u, used[link.1]? = some u ∧
        ∃ edge ∈ ((g.neighbors[v]?).getD []), edge.1 = u ∧ link.2 (bondCode g v edge.1 edge.2) = true) →
      Fits g (v :: used) ss vs → Fits g used (s :: ss) (v :: vs)

private theorem allowed_spec :
    allowed g s used v = true ↔ s.accepts v = true ∧ v ∉ used ∧
      ∀ link ∈ s.links, ∃ u, used[link.1]? = some u ∧
        ∃ edge ∈ ((g.neighbors[v]?).getD []), edge.1 = u ∧ link.2 (bondCode g v edge.1 edge.2) = true := by
  simp [allowed, Option.any, and_assoc]
  intro _ _
  constructor
  · intro h i p hi
    specialize h i p hi
    cases hv : used[i]? <;> simp_all
  · intro h i p hi
    obtain ⟨u, hv, e⟩ := h i p hi
    simpa [hv] using e

theorem mem_embeddings_iff (g : Graph) (used : List Nat) (steps : List Step) (xs : List Nat) :
    xs ∈ embeddings g used steps ↔ Fits g used steps xs := by
  induction steps generalizing used xs with
  | nil =>
    constructor
    · intro h
      simp [embeddings] at h
      subst xs
      exact .nil _
    · intro h
      cases h
      simp [embeddings]
  | cons s ss ih =>
    simp only [embeddings, List.mem_flatMap, List.mem_filter, List.mem_range, List.mem_map,
               allowed_spec, ih]
    constructor
    · rintro ⟨v, ⟨lt, ha, fresh, links⟩, vs, tail, rfl⟩
      exact .cons lt ha fresh links tail
    · intro h
      cases h with
      | cons lt ha fresh links tail => exact ⟨_, ⟨lt, ha, fresh, links⟩, _, tail, rfl⟩

def rooted (g : Graph) (v : Nat) : List Step → Bool
  | [] => false
  | s :: ss => v < g.atoms.size && allowed g s [] v && !(embeddings g [v] ss).isEmpty

theorem rooted_iff (g : Graph) (v : Nat) (s : Step) (ss : List Step) :
    rooted g v (s :: ss) = true ↔ ∃ vs, Fits g [] (s :: ss) (v :: vs) := by
  simp only [rooted, Bool.and_eq_true, decide_eq_true_eq, Bool.not_eq_true', allowed_spec]
  constructor
  · rintro ⟨⟨lt, ha, fresh, links⟩, nonempty⟩
    cases h : embeddings g [v] ss with
    | nil => simp [h] at nonempty
    | cons x xs =>
      exact ⟨x, .cons lt ha fresh links ((mem_embeddings_iff ..).mp (by simp [h]))⟩
  · rintro ⟨vs, h⟩
    cases h with
    | cons lt ha fresh links tail =>
      refine ⟨⟨lt, ha, fresh, links⟩, ?_⟩
      have hm := (mem_embeddings_iff ..).mpr tail
      cases he : embeddings g [v] ss <;> simp_all

private def distinctImages : List (List Nat) → List (List Nat)
  | [] => []
  | x :: xs =>
    let rest := distinctImages xs
    if x ∈ rest then rest else x :: rest

private theorem mem_distinctImages (x : List Nat) (xs : List (List Nat)) :
    x ∈ distinctImages xs ↔ x ∈ xs := by
  induction xs with
  | nil => simp [distinctImages]
  | cons y ys ih =>
    simp only [distinctImages]
    split <;> simp_all [List.mem_cons]
    intro eq
    subst x
    exact ih.mp (by assumption)

private theorem nodup_distinctImages (xs : List (List Nat)) :
    (distinctImages xs).Nodup := by
  induction xs with
  | nil => simp [distinctImages]
  | cons x xs ih =>
    simp only [distinctImages]
    split <;> simp_all

def uniqueImages (g : Graph) (query : List Step) : List (List Nat) :=
  distinctImages ((embeddings g [] query).map (fun vs => vs.mergeSort))

theorem uniqueImages_spec (g : Graph) (query : List Step) :
    (uniqueImages g query).Nodup ∧
    ∀ key, key ∈ uniqueImages g query ↔ ∃ vs, Fits g [] query vs ∧ vs.mergeSort = key := by
  refine ⟨nodup_distinctImages _, ?_⟩
  intro key
  simp only [uniqueImages, mem_distinctImages, List.mem_map, mem_embeddings_iff]

def fragmentCount (g : Graph) (query : List Step) : Nat :=
  min 1000 (uniqueImages g query).length

end Descriptors.Fragment
