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

private def walkColumns (g : Fragment.Graph) (vs : List Nat) (adj : Adj) : List (String × Json) := Id.run do
  let mut total := Array.replicate 11 0
  let mut closed := Array.replicate 11 0
  for j in vs do
    for k in [3:11] do
      let row := walkRow g.atoms.size adj j k
      total := total.modify k (· + (vs.map (fun i => row[i]!)).sum)
      closed := closed.modify k (· + row[j]!)
  return (List.range' 3 8).flatMap fun k =>
    [(s!"ds_MWC{k}", floatJson (Float.log ((total[k]! + 1).toFloat))),
     (s!"ds_SRW{k}", floatJson (Float.log ((closed[k]! + 1).toFloat)))]

def values (g : Fragment.Graph) : List (String × Json) := Id.run do
  let vs := heavyVertices g
  let n := g.atoms.size
  let adjRows := Array.ofFn (n := n) fun i => neighbors g i.val
  let adj := fun v => adjRows[v]!
  let d := distances n vs adj
  let pairs := vs.flatMap fun i => (vs.filterMap fun j =>
    if i < j && d[i]![j]! <= n then some (i, j, d[i]![j]!) else none)
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
