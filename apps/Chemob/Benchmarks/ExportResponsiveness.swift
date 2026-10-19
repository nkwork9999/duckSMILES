import Foundation
import ChemobCore
import Darwin

@MainActor private final class Heartbeat {
    var started = false
    var stop = false
    var maximumGap = 0.0
    var ticks = 0
    func run() async {
        var previous = ContinuousClock.now
        started = true
        repeat {
            try? await Task.sleep(for: .milliseconds(20))
            let now = ContinuousClock.now
            maximumGap = max(maximumGap, seconds(previous.duration(to: now)))
            previous = now
            ticks += 1
        } while !stop
    }
}

private func seconds(_ duration: Duration) -> Double {
    Double(duration.components.seconds) + Double(duration.components.attoseconds) / 1e18
}

// This adapter remains async in the regression harness even when the old
// production export method is synchronous and blocks its caller's main actor.
@MainActor private func exportSnapshot(_ model: Workbench) async throws -> URL {
    try await model.exportCSV()
}

@main struct ExportResponsiveness {
    @MainActor static func main() async throws {
        precondition(CommandLine.arguments.contains("--ui-testing"))
        let destination = URL(fileURLWithPath: CommandLine.arguments[1], isDirectory: true)
        try FileManager.default.createDirectory(at: destination, withIntermediateDirectories: true)
        // Repetition makes setup quick through the real cache. This is a
        // 10k-row CSV/responsiveness test, not distinct-molecule throughput.
        let text = (["name,smiles"] + (0..<10_000).map { index in
            "\"row \(index), \"\"quoted\"\"\",CCO"
        }).joined(separator: "\n")
        let expected = try CompoundImport.parse(text, format: "csv")
        let model = Workbench()
        model.importText(text, format: "csv", replacing: true)
        let deadline = ContinuousClock.now.advanced(by: .seconds(120))
        while model.busy {
            precondition(ContinuousClock.now < deadline)
            try await Task.sleep(for: .milliseconds(10))
        }
        precondition(model.error == nil && model.inputs == expected)
        precondition(model.results.count == 10_000 && model.results.allSatisfy(\.isValid))
        let before = await ChemistryStore.shared.statistics
        precondition(before.computedStructures == 1 && model.reusedCount == 9_999)
        let monitor = Heartbeat()
        let heartbeat = Task { await monitor.run() }
        // Start the clock before export. Recording a final tick after stop is
        // essential: starting after the old blocking call would falsely pass.
        while !monitor.started { await Task.yield() }
        let begin = ContinuousClock.now
        let exported = try await exportSnapshot(model)
        defer { model.discardCSV(exported) }
        let elapsed = seconds(begin.duration(to: .now))
        monitor.stop = true
        await heartbeat.value
        let data = try Data(contentsOf: exported)
        try data.write(to: destination.appending(path: "output.csv"))
        let encoder = JSONEncoder()
        encoder.outputFormatting = [.sortedKeys]
        try encoder.encode(model.results).write(to: destination.appending(path: "results.json"))
        let restored = try CompoundImport.parse(String(decoding: data, as: UTF8.self), format: "csv")
        precondition(restored == expected, "quoted names/source order must survive export/reimport")
        precondition(!model.exporting)
        let secondTask = Task { try await exportSnapshot(model) }
        while !model.exporting { await Task.yield() }
        model.inputs.reverse()
        let second = try await secondTask.value
        model.inputs = expected
        precondition(second != exported, "new exports must not overwrite an existing share URL")
        let secondData = try Data(contentsOf: second)
        let retainedData = try Data(contentsOf: exported)
        precondition(secondData == data, "an export must retain its starting snapshot while UI state changes")
        precondition(retainedData == data)
        model.discardCSV(second)
        precondition(!FileManager.default.fileExists(atPath: second.path))
        let canceled = Task { try await exportSnapshot(model) }
        while !model.exporting { await Task.yield() }
        model.cancelExport()
        do {
            let unwanted = try await canceled.value
            model.discardCSV(unwanted)
            preconditionFailure("Canceled export returned a share URL")
        } catch is CancellationError { }
        precondition(!model.exporting)
        let foreignDirectory = destination.appending(path: "unowned-probe", directoryHint: .isDirectory)
        try FileManager.default.createDirectory(at: foreignDirectory, withIntermediateDirectories: true)
        let foreign = foreignDirectory.appending(path: "keep.txt")
        try Data("KEEP".utf8).write(to: foreign)
        model.discardCSV(foreign)
        let retainedForeign = try Data(contentsOf: foreign)
        precondition(retainedForeign == Data("KEEP".utf8), "cleanup may remove only registered export URLs")
        try FileManager.default.removeItem(at: foreignDirectory)
        let after = await ChemistryStore.shared.statistics
        precondition(after.computedStructures == before.computedStructures, "export must not calculate molecules")
        let summary: [String: Any] = [
            "scope": "actual App Workbench, 10000 repeated ethanol rows; main-actor heartbeat; not a chemistry-speed benchmark",
            "engine_revision": Chemistry.revision, "descriptor_count": Chemistry.descriptorCatalog.count,
            "rows": expected.count, "distinct_smiles": 1,
            "csv_bytes": data.count, "export_seconds_observation": elapsed,
            "maximum_main_actor_heartbeat_gap_seconds": monitor.maximumGap,
            "heartbeat_ticks": monitor.ticks, "allowed_gap_seconds": 0.5,
            "csv_reimport_preserved_inputs": true, "new_native_calculations_during_export": 0,
            "unique_share_urls": true, "prior_export_content_preserved": true,
            "cancellation": "passed", "discarded_export_removed": true,
            "snapshot_preserved_during_state_change": true, "unregistered_file_preserved": true,
            "gate": monitor.maximumGap < 0.5 ? "passed" : "failed"
        ]
        let report = try JSONSerialization.data(withJSONObject: summary, options: [.prettyPrinted, .sortedKeys])
        try report.write(to: destination.appending(path: "summary.json"))
        print(String(decoding: report, as: UTF8.self))
        if monitor.maximumGap >= 0.5 { exit(1) }
    }
}
