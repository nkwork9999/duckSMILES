import Descriptors.Composition

namespace Descriptors.Connectivity
open Lean Fragment Composition Atomic

private def stepPaths (g : Graph) (visited : List Nat) (trail : List (Nat × Nat)) :
    Nat → List (List (Nat × Nat) × List Nat)
  | 0 => [((trail.mergeSort (fun a b => compare a b != .gt)), visited.reverse)]
  | k+1 => Id.run do
    let u := visited.head!
    let links := (((g.neighbors[u]?).getD []).mergeSort (fun a b => compare a b != .gt))
    return links.flatMap fun (v, _) =>
      let bond := (min u v, max u v)
      if (g.atom v).z == 1 || trail.contains bond || (k != 0 && visited.contains v) then []
      else stepPaths g (v :: visited) (bond :: trail) k

def paths (g : Graph) (k : Nat) : List (List Nat) := Id.run do
  let mut seen : Std.HashSet (List (Nat × Nat)) := {}
  let mut paths := []
  for i in [:g.atoms.size] do
    if (g.atom i).z != 1 then
      for (key, vs) in stepPaths g [i] [] k do
        if !seen.contains key then
          seen := seen.insert key
          paths := vs :: paths
  return paths.reverse

private def delta (x : Float) := if x == 0 then 0.0 else 1.0 / x.sqrt

private def divideSquare (a b c : Float) := if c == 0 then 0 else a*b*b/(c*c)

private def entropy (xs : List Float) :=
  let total := orderedSum xs
  if total == 0 then 0 else
    -(orderedSum (xs.filterMap (fun x =>
      let p := x/total
      if x == 0 then none else some (p * p.log)))) / Float.log 2

private def topology (g : Graph) (edges : List (Nat × Nat × Nat)) : List (String × Json) := Id.run do
  let n := g.atoms.size
  let inf := 600000000
  let mut d := Array.ofFn (n := n) fun i => Array.ofFn (n := n) fun j => if i == j then 0 else inf
  for (a,b,order) in edges do
    let weight := match order with | 2 => 3 | 3 => 2 | 12 => 4 | _ => 6
    d := (d.modify a (·.set! b weight)).modify b (·.set! a weight)
  for k in [:n] do
    let via := d[k]!
    for i in [:n] do
      let row := d[i]!
      let through := row[k]!
      if through < inf then
        d := d.set! i (row.mapIdx (fun j direct => min direct (through + via[j]!)))
  let sums := d.map fun row => row.toList.sum.toFloat / 6
  let q := edges.length.toFloat
  let denom := q - n.toFloat + 2
  let balaban := if denom == 0 then 0 else q/denom *
    orderedSum (edges.map fun (a,b,_) => delta (sums[a]! * sums[b]!))
  let signatures := d.map (fun row => row.toList.mergeSort.take 100)
  let types := g.atoms.toList.map Atom.z
  let atomEntropy := entropy (types.eraseDups.map fun z => (types.count z).toFloat)
  let mut connections : Std.HashMap (List Nat) Float := {}
  let order := fun (code : Nat) => if code == 12 then (1.5 : Float) else code.toFloat
  for i in [:n] do
    let center := signatures.toList.idxOf signatures[i]!
    let adj := g.neighbors[i]!
    for (j,b) in adj do
      let left := signatures.toList.idxOf signatures[j]!
      let v := order b
      if j > i && v > 1 then
        let key := [min center left, max center left]
        connections := connections.insert key (connections.getD key 0 + v * (v-1) / 2)
      for (k,c) in adj do
        if k > j then
          let right := signatures.toList.idxOf signatures[k]!
          let key := [min left right, center, max left right]
          connections := connections.insert key (connections.getD key 0 + v * order c)
  let weights := if connections.isEmpty then [1.0] else connections.toList.map Prod.snd
  let t := orderedSum weights
  let bertz := if n < 2 then 0 else n.toFloat * atomEntropy + t * (entropy weights + t.log / Float.log 2)
  return [("BalabanJ", floatJson balaban), ("BertzCT", floatJson bertz)]

def values (g : Graph) (edges : List (Nat × Nat × Nat)) (props : Array AtomProperties) :
    List (String × Json) := Id.run do
  let atoms := List.range g.atoms.size
  let degree := fun v => (g.atom v).degree.toFloat
  let path2 := paths g 2
  let path3 := paths g 3
  let path4 := paths g 4
  let a := (atoms.filter (fun i => (g.atom i).z > 1)).length
  let alpha := orderedSum (props.toList.map (fun x => decode x.alpha))
  let x := a.toFloat + alpha
  let k1 := divideSquare x (x-1) (edges.length.toFloat+alpha)
  let k2 := divideSquare (x-1) (x-2) (path2.length.toFloat+alpha)
  let k3 := divideSquare (if a % 2 == 1 then x-1 else x-2) (x-3) (path3.length.toFloat+alpha)
  let mut result := [("Kappa1", floatJson k1), ("Kappa2", floatJson k2), ("Kappa3", floatJson k3),
                     ("Phi", floatJson (if a == 0 then 0 else k1*k2/a.toFloat)),
                     ("Chi0", floatJson (orderedSum (atoms.map (fun i => delta (degree i))))),
                     ("Chi1", floatJson (orderedSum (edges.map (fun (i,j,_) => delta (degree i * degree j)))))]
  for variant in ["n", "v"] do
    let ds := Array.ofFn (n := g.atoms.size) fun i =>
      let prop := props[i.val]!
      let z := (g.atom i.val).z
      let raw := (prop.outer.toUInt32 - prop.attachedH.toUInt32).toNat.toFloat
      if variant == "n" then delta raw
      else if z <= 1 then 0 else if z <= 10 then delta raw
      else delta (raw / (z.toUInt32 - prop.outer.toUInt32 - 1).toNat.toFloat)
    result := ("Chi0" ++ variant, floatJson (orderedSum ds.toList)) ::
              ("Chi1" ++ variant, floatJson (orderedSum (edges.map fun (i,j,_) => ds[i]! * ds[j]!))) :: result
    for (k, ps) in [(2,path2),(3,path3),(4,path4)] do
      let vals := ps.map fun p =>
        let vs := if p.head? == p.getLast? then p.dropLast else p
        vs.foldl (fun prod i => prod * ds[i]!) 1.0
      result := (s!"Chi{k}{variant}", floatJson (orderedSum vals)) :: result
  return result ++ topology g edges

end Descriptors.Connectivity
