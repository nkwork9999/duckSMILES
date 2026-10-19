import Foundation
import ChemobCore

@main struct BulkIO {
    @MainActor static func main() async throws {
        let arguments = CommandLine.arguments
        guard arguments.count >= 3 else { fatalError("BulkIO input.csv output-directory --ui-testing") }
        let source = URL(fileURLWithPath: arguments[1])
        let destination = URL(fileURLWithPath: arguments[2])
        let clock = ContinuousClock()
        let begin = clock.now
        let text = try String(contentsOf: source, encoding: .utf8)
        let expected = try CompoundImport.parse(text, format: "csv")
        precondition(expected.count == 10_000)
        let model = Workbench()
        model.importText(text, format: "csv")
        let deadline = clock.now.advanced(by: .seconds(300))
        while model.busy {
            guard clock.now < deadline else { fatalError("Timed out waiting for native workbench") }
            try await Task.sleep(for: .milliseconds(5))
        }
        guard model.error == nil else { fatalError(model.error!) }
        precondition(model.inputs == expected, "import must preserve every row, name and structure")
        precondition(model.results.count == expected.count)
        precondition(model.results.map(\.index) == expected.map(\.id), "results must retain source row order")
        for input in expected { precondition(model.result(for: input) != nil) }
        let computed = clock.now
        let coldStatistics = await ChemistryStore.shared.statistics
        let warmBegin = clock.now
        let warmModel = Workbench()
        warmModel.importText(text, format: "csv")
        while warmModel.busy { try await Task.sleep(for: .milliseconds(5)) }
        precondition(warmModel.error == nil && warmModel.reusedCount == 10_000)
        let warmEnd = clock.now
        let warmStatistics = await ChemistryStore.shared.statistics
        precondition(coldStatistics.computedStructures == warmStatistics.computedStructures)
        let stableEncoder = JSONEncoder(); stableEncoder.outputFormatting = [.sortedKeys]
        let coldRows = try stableEncoder.encode(model.results)
        let warmRows = try stableEncoder.encode(warmModel.results)
        precondition(coldRows == warmRows, "cached results must equal fresh results")
        let exported = try await model.exportCSV()
        defer { model.discardCSV(exported) }
        let exportData = try Data(contentsOf: exported)
        try exportData.write(to: destination.appendingPathComponent("bulk-output.csv"))
        try JSONEncoder().encode(model.results).write(to: destination.appendingPathComponent("bulk-native-results.json"))
        let finished = clock.now
        var roundtrip = "passed"
        do {
            let reloaded = try CompoundImport.parse(String(decoding: exportData, as: UTF8.self), format: "csv")
            precondition(reloaded == expected, "export/reimport must preserve IDs, names and SMILES")
        } catch { roundtrip = "failed: \(error.localizedDescription)" }
        func milliseconds(_ a: ContinuousClock.Instant, _ b: ContinuousClock.Instant) -> Double {
            let d = a.duration(to: b).components
            return Double(d.seconds)*1000 + Double(d.attoseconds)/1e15
        }
        let summary: [String: Any] = [
            "platform": "macOS host using actual App/Chemistry.swift; no UI rendering",
            "engine_revision": Chemistry.revision, "descriptor_count": Chemistry.descriptorCatalog.count,
            "rows": model.inputs.count, "result_rows": model.results.count,
            "invalid_rows": model.results.filter { !$0.isValid }.count,
            "input_bytes": text.utf8.count, "export_bytes": exportData.count,
            "source_row_order_preserved": true, "input_names_and_smiles_preserved": true,
            "csv_reimport": roundtrip, "read_import_calculate_ms": milliseconds(begin, computed),
            "export_and_result_snapshot_ms": milliseconds(warmEnd, finished),
            "cold_unique_native_calculations": coldStatistics.computedStructures,
            "warm_additional_native_calculations": warmStatistics.computedStructures - coldStatistics.computedStructures,
            "warm_read_import_calculate_ms": milliseconds(warmBegin, warmEnd),
            "warm_reused_rows": warmModel.reusedCount,
            "engine_and_swift_workbench_seconds": model.elapsed ?? -1
        ]
        let data = try JSONSerialization.data(withJSONObject: summary, options: [.prettyPrinted,.sortedKeys])
        try data.write(to: destination.appendingPathComponent("bulk-io-summary.json"))
        print(String(decoding:data,as:UTF8.self))
        guard roundtrip == "passed" else { throw ImportIssue.tooLarge }
    }
}
