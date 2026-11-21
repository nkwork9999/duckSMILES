import Descriptors.Atomic

namespace Descriptors.DrugLikeness
open Lean Atomic

/-- Default QED weights and fitted ADS coefficients from RDKit 2026.03.6.
Property order is MW, ALOGP, HBA, HBD, PSA, ROTB, AROM, ALERTS.
See RDKIT-DATA-LICENSE.txt. -/
structure Parameter where
  weight : Float
  a : Float
  b : Float
  c : Float
  d : Float
  e : Float
  f : Float
  scale : Float

def parameters : List Parameter := [
  ⟨0.66, 2.817065973, 392.5754953, 290.7489764, 2.419764353, 49.22325677, 65.37051707, 104.9805561⟩,
  ⟨0.46, 3.172690585, 137.8624751, 2.534937431, 4.581497897, 0.822739154, 0.576295591, 131.3186604⟩,
  ⟨0.05, 2.948620388, 160.4605972, 3.615294657, 4.435986202, 0.290141953, 1.300669958, 148.7763046⟩,
  ⟨0.61, 1.618662227, 1010.051101, 0.985094388, 1e-9, 0.713820843, 0.920922555, 258.1632616⟩,
  ⟨0.06, 1.876861559, 125.2232657, 62.90773554, 87.83366614, 12.01999824, 28.51324732, 104.5686167⟩,
  ⟨0.65, 0.01, 272.4121427, 2.55837997, 1.565547684, 1.271567166, 2.758063707, 105.4420403⟩,
  ⟨0.48, 3.21778897, 957.7374108, 2.274627939, 1e-9, 1.317690384, 0.375760881, 312.337261⟩,
  ⟨0.95, 0.01, 1199.094025, -0.09002883, 1e-9, 0.185904477, 0.875193782, 417.725314⟩]

def desirability (p : Parameter) (x : Float) :=
  let left := 1 + Float.exp (-(x - p.c + p.d/2)/p.e)
  let right := 1 + Float.exp (-(x - p.c - p.d/2)/p.f)
  (p.a + p.b / left * (1 - 1 / right)) / p.scale

def value (properties : List Nat) : Except String Json := do
  unless properties.length == parameters.length do throw "QED requires eight properties"
  let terms := (parameters.zip properties).map fun (p,bits) =>
    p.weight * Float.log (desirability p (decode bits))
  return floatJson (Float.exp (orderedSum terms / orderedSum (parameters.map Parameter.weight)))

end Descriptors.DrugLikeness
