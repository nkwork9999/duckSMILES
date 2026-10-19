import Foundation

public struct CompoundInput: Identifiable, Codable, Equatable, Sendable {
    public var id: Int
    public var name: String
    public var smiles: String
    public init(id: Int, name: String, smiles: String) { self.id = id; self.name = name; self.smiles = smiles }
}

public enum CompoundImport {
    public static let maximumInputBytes = 32 * 1_048_576
    public static func appending(_ incoming: [CompoundInput], to existing: [CompoundInput]) throws -> [CompoundInput] {
        guard existing.count + incoming.count <= 10_000 else { throw ImportIssue.tooLarge }
        let lastID = existing.map(\.id).max() ?? -1
        guard lastID < Int.max - incoming.count else { throw ImportIssue.tooLarge }
        return existing + incoming.enumerated().map { offset, value in
            CompoundInput(id: lastID + offset + 1, name: value.name, smiles: value.smiles)
        }
    }
    public static func parse(_ text: String, format: String = "smi") throws -> [CompoundInput] {
        guard text.utf8.count <= maximumInputBytes else { throw ImportIssue.tooLarge }
        let text = text.replacingOccurrences(of: "\u{FEFF}", with: "")
        if format.lowercased() == "csv" || format.lowercased() == "tsv" {
            let rows = try fields(text, delimiter: format.lowercased() == "tsv" ? "\t" : ",")
            guard let header = rows.first else { return [] }
            let names = header.map { $0.trimmingCharacters(in: .whitespacesAndNewlines).lowercased() }
            guard let structure = names.firstIndex(where: { ["smiles", "canonical_smiles", "isomeric_smiles"].contains($0) }) else { throw ImportIssue.noStructure }
            let name = names.firstIndex(where: { ["name", "compound", "id"].contains($0) })
            return try bounded(rows.dropFirst().enumerated().filter { !$0.element.allSatisfy(\.isEmpty) }.map { index, row in
                CompoundInput(id: index, name: name.flatMap { $0 < row.count ? row[$0] : nil } ?? "Compound \(index + 1)", smiles: structure < row.count ? row[structure].trimmingCharacters(in: .whitespacesAndNewlines) : "")
            })
        }
        return try bounded(text.components(separatedBy: .newlines).filter { !$0.trimmingCharacters(in: .whitespaces).isEmpty }.enumerated().map { index, line in
            let parts = line.split(maxSplits: 1, whereSeparator: \.isWhitespace)
            let remainder = parts.count > 1 ? parts[1].drop(while: \.isWhitespace) : Substring()
            if remainder.first == "|" {
                // A CX block belongs to the structure, before the optional name.
                // Preserve malformed and unknown blocks for engine validation;
                // accepting their ordinary-SMILES prefix would lose chemistry.
                if let close = remainder.dropFirst().firstIndex(of: "|") {
                    let end = line.index(after: close)
                    let name = line[end...].trimmingCharacters(in: .whitespaces)
                    if !name.hasPrefix("|") {
                        return CompoundInput(id: index,
                            name: name.isEmpty ? "Compound \(index + 1)" : name,
                            smiles: String(line[parts[0].startIndex..<end]))
                    }
                }
                return CompoundInput(id: index, name: "Compound \(index + 1)",
                    smiles: String(line[parts[0].startIndex...]))
            }
            return CompoundInput(id: index, name: parts.count > 1 ? String(parts[1]) : "Compound \(index + 1)", smiles: String(parts[0]))
        })
    }

    private static func bounded(_ values: [CompoundInput]) throws -> [CompoundInput] {
        guard values.count <= 10_000 else { throw ImportIssue.tooLarge }
        return values
    }

    private static func fields(_ text: String, delimiter: Character) throws -> [[String]] {
        var rows: [[String]] = [], row: [String] = [], value: [UInt8] = [], quoted = false
        // CSV syntax is ASCII; scanning UTF-8 avoids a Character array whose
        // storage can be many times larger than a wide descriptor export.
        let bytes = Array(text.utf8)
        let separator: UInt8 = delimiter == "\t" ? 9 : 44
        var i = 0
        func finishField() {
            row.append(String(decoding: value, as: UTF8.self))
            value.removeAll(keepingCapacity: true)
        }
        while i < bytes.count {
            let c = bytes[i]
            if c == 34 {
                if quoted && i + 1 < bytes.count && bytes[i + 1] == 34 { value.append(34); i += 1 }
                else { quoted.toggle() }
            } else if c == separator && !quoted { finishField() }
            else if (c == 10 || c == 13) && !quoted {
                finishField(); rows.append(row); row = []
                if c == 13 && i + 1 < bytes.count && bytes[i + 1] == 10 { i += 1 }
            } else { value.append(c) }
            i += 1
        }
        guard !quoted else { throw ImportIssue.unclosedQuote }
        if !value.isEmpty || !row.isEmpty { finishField(); rows.append(row) }
        return rows
    }
}

public enum ImportIssue: LocalizedError {
    case noStructure, unclosedQuote, tooLarge, empty
    public var errorDescription: String? {
        switch self {
        case .empty: NSLocalizedString("No compounds were found in the input.", comment: "Empty import")
        case .noStructure: NSLocalizedString("No SMILES column found. Name the column smiles.", comment: "Import error")
        case .unclosedQuote: NSLocalizedString("Unclosed quotation mark in CSV.", comment: "Import error")
        case .tooLarge: NSLocalizedString("Import up to 32 MiB and 10,000 molecules.", comment: "Import error")
        }
    }
}
