import Foundation
import ChemobCore

struct CachedAnalysis: Sendable {
    let rows: [AnalysisRow]
    let reused: Int
}

/// Native calls, cache decoding and disk I/O run on this actor, away from the UI.
actor ChemistryStore {
    struct Statistics: Codable, Sendable {
        var computedStructures = 0
        var graphBuilds = 0
        var drawingBuilds = 0
    }
    private struct GraphEntry: Codable, Sendable { let value: MoleculeGraph? }
    private struct DrawingEntry: Codable, Sendable { let value: ChemicalDepiction? }
    static let shared = ChemistryStore(directory: ProcessInfo.processInfo.arguments.contains("--ui-testing") ? nil : URL.applicationSupportDirectory.appending(path: "Chemob/chemistry-cache"))
    private let directory: URL?
    private var descriptors: MoleculeMemo<AnalysisRow>
    private var graphs: MoleculeMemo<GraphEntry>
    private var drawings: MoleculeMemo<DrawingEntry>
    private var loaded = false
    private var descriptorsDirty = false
    private var graphsDirty = false
    private var drawingsDirty = false
    private(set) var statistics = Statistics()

    init(directory: URL?) {
        self.directory = directory
        descriptors = MoleculeMemo(namespace: Chemistry.revision + ":descriptors", capacity: 10_000)
        graphs = MoleculeMemo(namespace: Chemistry.revision + ":drawing-graph-prepared-v1", capacity: 256)
        drawings = MoleculeMemo(namespace: Chemistry.revision + ":depiction-2d-v4", capacity: 256)
    }

    func analyze(_ inputs: [CompoundInput]) throws -> CachedAnalysis {
        try Task.checkCancellation()
        loadIfNeeded()
        var output = [AnalysisRow?](repeating: nil, count: inputs.count)
        var positions: [String: [Int]] = [:]
        var missing: [CompoundInput] = []
        var reused = 0
        for (index, input) in inputs.enumerated() {
            if var value = descriptors.value(for: input.smiles) {
                value.index = index; output[index] = value; reused += 1
            } else {
                if positions[input.smiles] == nil { missing.append(input) }
                else { reused += 1 }
                positions[input.smiles, default: []].append(index)
            }
        }
        for start in stride(from: 0, to: missing.count, by: 64) {
            try Task.checkCancellation()
            let batch = Array(missing[start..<min(start+64, missing.count)])
            let values = try Chemistry.analyze(batch)
            guard values.count == batch.count else { throw EngineIssue.rejected }
            for var value in values {
                guard batch.indices.contains(value.index) else { throw EngineIssue.rejected }
                let smiles = batch[value.index].smiles
                value.index = 0
                descriptors.insert(value, for: smiles)
                descriptorsDirty = true
                statistics.computedStructures += 1
                for index in positions[smiles] ?? [] {
                    value.index = index; output[index] = value
                }
            }
        }
        try Task.checkCancellation()
        guard output.allSatisfy({ $0 != nil }) else { throw EngineIssue.rejected }
        return CachedAnalysis(rows: output.compactMap { $0 }, reused: reused)
    }

    func graph(_ smiles: String) -> MoleculeGraph? {
        guard !Task.isCancelled else { return nil }
        loadIfNeeded()
        if let stored = graphs.value(for: smiles) { return stored.value }
        let value = Chemistry.graph(smiles)
        statistics.graphBuilds += 1
        graphs.insert(GraphEntry(value: value), for: smiles)
        graphsDirty = true
        save()
        return value
    }

    func drawing(_ smiles: String) -> ChemicalDepiction? {
        guard !Task.isCancelled else { return nil }
        loadIfNeeded()
        if let stored = drawings.value(for: smiles) { return stored.value }
        let value = graph(smiles).flatMap(ChemicalDepiction.make)
        guard !Task.isCancelled else { return nil }
        statistics.drawingBuilds += 1
        drawings.insert(DrawingEntry(value: value),for: smiles)
        drawingsDirty = true
        save()
        return value
    }

    func save() {
        guard let directory else { return }
        do {
            try FileManager.default.createDirectory(at: directory, withIntermediateDirectories: true)
            if descriptorsDirty {
                try descriptors.snapshot().write(to: directory.appending(path: "descriptors.json"), options: .atomic)
                descriptorsDirty = false
            }
            if graphsDirty {
                try graphs.snapshot().write(to: directory.appending(path: "graphs.json"), options: .atomic)
                graphsDirty = false
            }
            if drawingsDirty {
                try drawings.snapshot().write(to: directory.appending(path: "drawings.json"), options: .atomic)
                drawingsDirty = false
            }
        } catch { /* Cache failure must never prevent calculation or original-input saving. */ }
    }

    private func loadIfNeeded() {
        guard !loaded else { return }; loaded = true
        guard let directory else { return }
        func data(_ file: String) -> Data? {
            let url = directory.appending(path: file)
            guard let size = try? url.resourceValues(forKeys: [.fileSizeKey]).fileSize,
                  size <= 32 * 1_048_576 else { return nil }
            return try? Data(contentsOf: url)
        }
        if let bytes = data("descriptors.json") { try? descriptors.restore(bytes) }
        if let bytes = data("graphs.json") { try? graphs.restore(bytes) }
        if let bytes = data("drawings.json") { try? drawings.restore(bytes) }
    }
}
