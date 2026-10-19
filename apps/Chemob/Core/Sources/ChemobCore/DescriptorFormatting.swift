import Foundation

public enum DescriptorFormatting {
    public static func text(_ value: Double?, digits: Int, locale: Locale = .current) -> String {
        guard let value, value.isFinite else { return "—" }
        let precision = max(0, min(12, digits))
        var style = FloatingPointFormatStyle<Double>.number
            .grouping(.never).precision(.fractionLength(precision)).locale(locale)
        if abs(value) >= 1e8 { style = style.notation(.scientific) }
        let rendered = (value == 0 ? 0.0 : value).formatted(style)
        // Let the locale and Foundation rounding rules decide whether every
        // displayed digit is zero, including values such as -0.0014 at 2 digits.
        if value < 0 && !rendered.contains(where: { ($0.wholeNumberValue ?? 0) != 0 }) {
            return 0.0.formatted(style)
        }
        return rendered
    }
}
