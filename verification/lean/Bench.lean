import Descriptors

def main (args : List String) : IO Unit := do
  let [path, count] := args | throw <| IO.userError "bench_counts ATOMS REPETITIONS"
  let some n := count.toNat? | throw <| IO.userError "invalid repetition count"
  let rows := ((← IO.FS.readFile path).trim.splitOn "\n").map (·.splitOn " ")
  let probes := #["C", "N", "O", "Cl", "F"]
  let start ← IO.monoNanosNow
  let mut checksum := 0
  for i in [:n] do
    let key := probes[i % probes.size]!
    for atoms in rows do
      checksum := checksum + Descriptors.countAtoms atoms key
  let stop ← IO.monoNanosNow
  IO.println s!"{stop - start},{checksum}"
