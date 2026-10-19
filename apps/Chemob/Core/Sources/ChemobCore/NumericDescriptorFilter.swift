import Foundation

public enum NumericFilterIssue: String, LocalizedError, Sendable {
    case number = "Enter a valid number."
    case empty = "Enter a minimum or maximum."
    case order = "Minimum must not exceed maximum."
    case property = "Choose a property."
    public var errorDescription: String? { NSLocalizedString(rawValue, comment: "Numeric filter validation") }
}

/// A range over cached descriptor values. Missing/undefined values never
/// become zero; both endpoints are inclusive. No molecule calculation here.
public struct NumericDescriptorFilter: Equatable, Sendable {
    public let descriptorID: String
    public let minimum: Double?
    public let maximum: Double?

    public init(descriptorID: String, minimum: Double?, maximum: Double?) throws {
        guard !descriptorID.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty else { throw NumericFilterIssue.property }
        guard minimum != nil || maximum != nil else { throw NumericFilterIssue.empty }
        guard minimum?.isFinite != false, maximum?.isFinite != false else { throw NumericFilterIssue.number }
        if let minimum, let maximum, minimum > maximum { throw NumericFilterIssue.order }
        self.descriptorID = descriptorID; self.minimum = minimum; self.maximum = maximum
    }

    public func contains(_ value: Double?) -> Bool {
        guard let value, value.isFinite else { return false }
        if let minimum, value < minimum { return false }
        if let maximum, value > maximum { return false }
        return true
    }

    /// Scientific notation and the locale decimal separator are accepted.
    /// Group separators and partially parsed strings are rejected explicitly.
    public static func parseBound(_ text: String, locale: Locale = .current) throws -> Double? {
        var text = text.trimmingCharacters(in: .whitespacesAndNewlines)
        guard !text.isEmpty else { return nil }
        let separator = locale.decimalSeparator ?? "."
        if separator != "." {
            if text.contains(separator) && text.contains(".") { throw NumericFilterIssue.number }
            text = text.replacingOccurrences(of: separator, with: ".")
        }
        let grammar = #"\A[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?\z"#
        guard text.range(of: grammar, options: .regularExpression) != nil,
              let value = Double(text), value.isFinite else { throw NumericFilterIssue.number }
        return value
    }
}
