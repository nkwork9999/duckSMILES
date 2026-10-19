import Foundation
import ChemobCore

/// Exercises the real disk store at capacity, with original inputs retained.
@main struct CacheBulk {
    static func main() async throws {
        guard CommandLine.arguments.count == 3 else {
            fatalError("CacheBulk input.csv summary.json")
        }
        let text = try String(contentsOfFile: CommandLine.arguments[1], encoding: .utf8)
        let inputs = try CompoundImport.parse(text, format: "csv")
        precondition(inputs.count == 10_000 && Set(inputs.map(\.smiles)).count == 10_000)
        let directory = URL.temporaryDirectory.appending(path: "chemob-bulk-cache-" + UUID().uuidString)
        defer { try? FileManager.default.removeItem(at: directory) }
        let first = ChemistryStore(directory: directory)
        let fresh = try await first.analyze(inputs)
        precondition(fresh.rows.count == inputs.count && fresh.reused == 0)
        precondition(fresh.rows.map(\.index) == Array(inputs.indices))
        await first.save()
        let url = directory.appending(path: "descriptors.json")
        let bytes = try Data(contentsOf: url).count
        let reopened = ChemistryStore(directory: directory)
        let restored = try await reopened.analyze(inputs)
        let statistics = await reopened.statistics
        precondition(restored.reused == 10_000 && statistics.computedStructures == 0)
        let encoder = JSONEncoder(); encoder.outputFormatting = [.sortedKeys]
        let freshData = try encoder.encode(fresh.rows)
        let restoredData = try encoder.encode(restored.rows)
        precondition(freshData == restoredData, "disk reload must preserve every field and source row")
        let summary: [String: Any] = [
            "platform": "macOS host, actual ChemistryStore, no UI",
            "revision": Chemistry.revision,
            "input_rows": inputs.count, "distinct_smiles": Set(inputs.map(\.smiles)).count,
            "invalid_rows": fresh.rows.filter { !$0.isValid }.count,
            "cache_bytes": bytes, "disk_reused_rows": restored.reused,
            "additional_native_calculations": statistics.computedStructures,
            "all_restored_values_and_row_order_equal": true
        ]
        let data = try JSONSerialization.data(withJSONObject: summary, options: [.prettyPrinted, .sortedKeys])
        try data.write(to: URL(fileURLWithPath: CommandLine.arguments[2]))
        print(String(decoding: data, as: UTF8.self))
    }
}
