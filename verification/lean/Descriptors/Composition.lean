import Descriptors.Fragment
import Descriptors.Atomic

namespace Descriptors.Composition
open Lean Fragment Atomic

structure AtomProperties where
  mass : Nat
  exactMass : Nat
  attachedH : Nat
  outer : Nat
  radical : Nat
  tpsa : Nat
  alpha : Nat
  hybrid : String
  stereoPossible : Bool
  cipAssigned : Bool
  modernStereo : Bool
  deriving FromJson, Inhabited

/-- Atomwise molecular mass including attached hydrogens; separate hydrogen
vertices already carry an atomic mass in the initial sum. -/
def massSum (mass : AtomProperties → Nat) (h : Float) (xs : List AtomProperties) :=
  orderedSum (xs.map (fun a => decode (mass a))) + (xs.map AtomProperties.attachedH).sum.toFloat * h

def shapeScore (g : Graph) (xs : Array AtomProperties) (i : Nat) :=
  let a := g.atom i
  let props := xs[i]!
  let hybrid := match props.hybrid with | "SP" => 1 | "SP2" => 2 | "SP3" => 3 | _ => 4
  hybrid * (if props.modernStereo then 2 else 1) *
    (if a.rings > 0 && !a.aromatic then 2 else 1) * a.degree^2

def values (g : Graph) (props : Array AtomProperties) (expandedLogP expandedMR : List Nat) :
    List (String × Json) := Id.run do
  let ids := List.range g.atoms.size
  let atoms := props.toList
  let heavy := ids.filter fun i => (g.atom i).z > 1
  let carbon := ids.filter fun i => (g.atom i).z == 6
  let charge := (g.atoms.toList.map Atom.charge).sum
  let noH := ids.filter fun i => [7,8].contains (g.atom i).z
  let chargeF := if charge < 0 then -charge.natAbs.toFloat else charge.natAbs.toFloat
  let stereo := atoms.filter (fun a => a.stereoPossible || a.cipAssigned)
  return [
    ("MolWt", floatJson (massSum AtomProperties.mass 1.008 atoms)),
    ("ExactMolWt", floatJson (massSum AtomProperties.exactMass 1.007825032 atoms - 0.00054857991 * chargeF)),
    ("HeavyAtomMolWt", floatJson (orderedSum (heavy.map fun i => decode props[i]!.mass))),
    ("NumValenceElectrons", toJson ((Int.ofNat (atoms.map (fun a => a.outer + a.attachedH)).sum) - charge)),
    ("NHOHCount", toJson ((noH.map fun i => (g.atom i).hydrogens).sum)),
    ("NumRadicalElectrons", toJson ((atoms.map AtomProperties.radical).sum)),
    ("FractionCSP3", floatJson (if carbon.isEmpty then 0 else
      ((carbon.filter fun i => (g.atom i).totalDegree == 4).length.toFloat / carbon.length.toFloat))),
    ("MolLogP", floatJson (orderedSum (expandedLogP.map decode))),
    ("MolMR", floatJson (orderedSum (expandedMR.map decode))),
    ("TPSA", floatJson (orderedSum (atoms.map fun a => decode a.tpsa))),
    ("HallKierAlpha", floatJson (orderedSum (atoms.map fun a => decode a.alpha))),
    ("NumAtomStereoCenters", toJson stereo.length),
    ("NumUnspecifiedAtomStereoCenters", toJson (stereo.filter (fun a => !a.cipAssigned)).length),
    ("SPS", if heavy.isEmpty then Json.null else floatJson (((ids.map (shapeScore g props)).sum).toFloat / heavy.length.toFloat))]

/-- Endpoints of the common paths of two basis rings have degree one in
the subgraph of shared bonds. Two rings sharing one vertex count as spiro. -/
def ringJunctions (rings : List (List Nat)) (edges : Array (Nat × Nat × Nat)) : List (String × Json) := Id.run do
  let endpoints := fun (bs : List Nat) => bs.flatMap fun b =>
    let (i,j,_) := edges[b]!
    [i,j]
  let mut spiro := []
  let mut bridges := []
  for i in [:rings.length] do
    let r := rings[i]!
    let atoms := (endpoints r).eraseDups
    for s in rings.drop (i+1) do
      let common := atoms.filter (endpoints s).contains
      if let [v] := common then spiro := v :: spiro
      let shared := r.filter s.contains
      if shared.length >= 2 then
        let ends := endpoints shared
        bridges := (ends.filter fun v => ends.count v == 1) ++ bridges
  return [("NumSpiroAtoms", toJson spiro.eraseDups.length),
          ("NumBridgeheadAtoms", toJson bridges.eraseDups.length)]

end Descriptors.Composition
