import Lean

namespace Descriptors

/-- Encode a finite IEEE binary64 value as an exact JSON decimal. This avoids
rounding the reference to six decimal places before applying catalog tolerances. -/
def floatJson (x : Float) : Lean.Json := Id.run do
  if !x.isFinite then return .null
  let bits := x.toBits.toNat
  let exponent := bits / 2^52 % 2048
  let fraction := bits % 2^52
  let significand := if exponent == 0 then fraction else 2^52 + fraction
  if significand == 0 then return Lean.toJson (0 : Nat)
  let shift := max 1 exponent
  let places := 1075 - shift
  let unsigned : Nat := significand * 2^(shift - 1075) * 5^places
  let signed : Int := if bits >= 2^63 then -(unsigned : Int) else unsigned
  return .num ⟨signed, places⟩

end Descriptors
