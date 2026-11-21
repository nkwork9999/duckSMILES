import Descriptors.Distance
import Descriptors.Fragment
import Descriptors.Numerical

namespace Descriptors.GraphIndices
open Lean

/-- Vertices of the hydrogen-depleted graph retain their original atom indices. -/
def heavyVertices (g : Fragment.Graph) : List Nat :=
  (List.range g.atoms.size).filter fun v => (g.atom v).z != 1

abbrev Adj := Nat → List Nat

def neighbors (g : Fragment.Graph) : Adj := fun i =>
  (((g.neighbors[i]?).getD []).map Prod.fst).filter fun v => (g.atom v).z != 1

def walkStep (adj : Adj) (previous : Nat → Nat) (i : Nat) : Nat :=
  ((adj i).map previous).sum

/-- Explicit endpoints of all length-k walks from i with multiplicity. -/
def endpoints (adj : Adj) : Nat → Nat → List Nat
  | 0, i => [i]
  | k+1, i => (adj i).flatMap (endpoints adj k)

private theorem count_flatMap (vs : List Nat) (f : Nat → List Nat) (j : Nat) :
    (vs.flatMap f).count j = (vs.map (fun v => (f v).count j)).sum := by
  induction vs with
  | nil => rfl
  | cons x xs ih => simp [ih]

def walkRow (n : Nat) (adj : Adj) (j : Nat) : Nat → Array Nat
  | 0 => Array.ofFn (n := n) fun i => if i.val = j then 1 else 0
  | k+1 =>
    let previous := walkRow n adj j k
    Array.ofFn (n := n) fun i => walkStep adj (fun v => previous[v]?.getD 0) i

theorem walkRow_counts (n : Nat) (adj : Adj)
    (closed : ∀ i < n, ∀ v ∈ adj i, v < n) (j k i : Nat) (hi : i < n) :
    (walkRow n adj j k)[i]?.getD 0 = (endpoints adj k i).count j := by
  induction k generalizing i with
  | zero =>
    by_cases h : i = j
    · subst j
      simp [walkRow, endpoints, hi, List.count_cons]
    · have h' := Ne.symm h
      simp [walkRow, endpoints, hi, List.count_cons, h, h', beq_iff_eq]
  | succ k ih =>
    simp [walkRow, hi, endpoints, count_flatMap, walkStep]
    congr 1
    exact List.map_congr_left fun v hv => ih v (closed i hi v hv)

/-- Floyd-Warshall on heavy atoms. The sentinel is greater than every
finite shortest simple path, including in disconnected inputs. -/
def distances (n : Nat) (vertices : List Nat) (adj : Adj) : Array (Array Nat) := Id.run do
  let inf := n + 1
  let mut d := Array.ofFn (n := n) fun i => Array.ofFn (n := n) fun j =>
    if i == j then 0 else if (adj i.val).contains j.val then 1 else inf
  for k in vertices do
    let via := d[k]!
    for i in vertices do
      let di := d[i]!
      let dik := di[k]!
      if dik < inf then
        d := d.set! i (di.mapIdx fun j direct => min direct (dik + via[j]!))
  return d

def mdeClasses : List (String × Nat × Nat × Nat) :=
  [("C", 6, 4), ("N", 7, 3), ("O", 8, 2)].flatMap fun (label, z, limit) =>
    (List.range' 1 limit).flatMap fun a => (List.range' a (limit + 1 - a)).map fun b =>
      (s!"ds_MDE{label}{a}{b}", z, a, b)

/-- Ordered pairs from a list containing each heavy vertex once.
The sentinel is excluded so disconnected components contribute separately. -/
def finitePairs (n : Nat) (vertices : List Nat) (d : Nat → Nat → Nat) :
    List (Nat × Nat × Nat) :=
  vertices.flatMap fun i => vertices.filterMap fun j =>
    if i < j && d i j ≤ n then some (i, j, d i j) else none

theorem mem_finitePairs (n : Nat) (vertices : List Nat) (d : Nat → Nat → Nat) (i j k : Nat) :
    (i,j,k) ∈ finitePairs n vertices d ↔
    i ∈ vertices ∧ j ∈ vertices ∧ i < j ∧ k = d i j ∧ k ≤ n := by
  simp only [finitePairs, List.mem_flatMap, List.mem_filterMap]
  constructor
  · rintro ⟨a, ha, b, hb, pair⟩
    split at pair
    · rename_i h
      simp only [Bool.and_eq_true] at h
      obtain ⟨lt, le⟩ := h
      simp only [decide_eq_true_eq] at lt le
      cases pair
      exact ⟨ha, hb, lt, rfl, le⟩
    · contradiction
  · rintro ⟨hi, hj, lt, rfl, le⟩
    exact ⟨i, hi, j, hj, by simp [lt, le]⟩

def allWalks (n : Nat) (vs : List Nat) (adj : Adj) (k : Nat) : Nat :=
  (vs.map (fun j => let row := walkRow n adj j k
                    (vs.map (fun i => row[i]?.getD 0)).sum)).sum

def closedWalks (n : Nat) (vs : List Nat) (adj : Adj) (k : Nat) : Nat :=
  (vs.map (fun j => (walkRow n adj j k)[j]?.getD 0)).sum

theorem allWalks_counts (n : Nat) (vs : List Nat) (adj : Adj) (k : Nat)
    (within : ∀ i ∈ vs, i < n)
    (closed : ∀ i < n, ∀ v ∈ adj i, v < n) :
    allWalks n vs adj k =
      (vs.map (fun j => (vs.map (fun i => (endpoints adj k i).count j)).sum)).sum := by
  unfold allWalks
  congr 1
  apply List.map_congr_left
  intro j _
  dsimp
  congr 1
  exact List.map_congr_left fun i hi => walkRow_counts n adj closed j k i (within i hi)

theorem closedWalks_counts (n : Nat) (vs : List Nat) (adj : Adj) (k : Nat)
    (within : ∀ i ∈ vs, i < n)
    (closed : ∀ i < n, ∀ v ∈ adj i, v < n) :
    closedWalks n vs adj k = (vs.map (fun i => (endpoints adj k i).count i)).sum := by
  unfold closedWalks
  congr 1
  exact List.map_congr_left fun i hi => walkRow_counts n adj closed i k i (within i hi)

private def walkColumns (g : Fragment.Graph) (vs : List Nat) (adj : Adj) : List (String × Json) :=
  (List.range' 3 8).flatMap fun k =>
    [(s!"ds_MWC{k}", floatJson (Float.log ((allWalks g.atoms.size vs adj k + 1).toFloat))),
     (s!"ds_SRW{k}", floatJson (Float.log ((closedWalks g.atoms.size vs adj k + 1).toFloat)))]

#print axioms mem_finitePairs
#print axioms allWalks_counts
#print axioms closedWalks_counts

def values (g : Fragment.Graph) : Except String (List (String × Json)) := do
  let vs := heavyVertices g
  let n := g.atoms.size
  let adjRows := Array.ofFn (n := n) fun i =>
    if (g.atom i.val).z == 1 then [] else neighbors g i.val
  let adj := fun v => adjRows[v]!
  let d := distances n vs adj
  for i in vs do
    unless Distance.accepts n adj i (fun j => d[i]![j]!) do
      throw s!"invalid heavy-atom distance row {i}"
  let pairs := finitePairs n vs (fun i j => d[i]![j]!)
  let degree := fun i => (adj i).length
  let mut result := walkColumns g vs adj
  let mut wiener := 0
  let mut degreeDistance := 0
  let mut gutman := 0
  let mut byDistance := Array.replicate (n+1) 0
  let mut ggi := Array.replicate 11 (0 : Float)
  let transfer := fun i j => (adj i).foldl (init := 0.0) fun acc v =>
    let dist := d[v]![j]!
    if dist == 0 || dist > n then acc else acc + 1.0 / (dist * dist).toFloat
  for (i,j,length) in pairs do
    wiener := wiener + length
    degreeDistance := degreeDistance + (degree i + degree j) * length
    gutman := gutman + degree i * degree j * length
    byDistance := byDistance.modify length (· + 1)
    if length <= 10 then
      ggi := ggi.modify length (· + Float.abs (transfer i j - transfer j i))
  let eccentric := vs.foldl (init := 0) fun acc i =>
    acc + degree i * ((d[i]!.toList.filter (· <= n)).foldl max 0)
  let harary := (List.range' 1 n).foldl (init := 0.0) fun acc k => acc + byDistance[k]!.toFloat / k.toFloat
  result := result ++ [("ds_Wiener", toJson wiener), ("ds_DegreeDistance", toJson degreeDistance),
                       ("ds_Gutman", toJson gutman), ("ds_EccentricConnectivity", toJson eccentric),
                       ("ds_Harary", floatJson harary)]
  for k in [1:11] do result := (s!"ds_GGI{k}", floatJson ggi[k]!) :: result
  for (id,z,a,b) in mdeClasses do
    let ds := pairs.filterMap fun (i,j,length) =>
      if (g.atom i).z == z && (g.atom j).z == z &&
         min (degree i) (degree j) == a && max (degree i) (degree j) == b
      then some length else none
    let count := ds.length.toFloat
    let logs := ds.foldl (init := 0.0) fun acc k => acc + Float.log k.toFloat
    result := (id, if ds.isEmpty then Json.null else floatJson (count * Float.exp (-logs/count))) :: result
  return if vs.isEmpty then result.map (fun kv => (kv.1, Json.null)) else result

end Descriptors.GraphIndices
