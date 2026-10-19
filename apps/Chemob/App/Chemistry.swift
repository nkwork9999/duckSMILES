import Foundation
import Combine
import ChemobCore

struct AnalysisRow: Codable, Sendable {
    var index: Int
    let status: String
    let formula: String?
    let mw: Double?
    let exact_mass: Double?
    let heavy_atoms: Int?
    let bonds: Int?
    let formal_charge: Int?
    let fragments: Int?
    let hbd: Int?
    let hba: Int?
    let tpsa: Double?
    let logp: Double?
    let rotatable_bonds: Int?
    let fraction_csp3: Double?
    let descriptor_values: [Double?]?
    var isValid: Bool { status == "ok" }
}

struct AnalysisEnvelope: Decodable { let rows: [AnalysisRow] }
enum EngineIssue: LocalizedError {
    case rejected
    var errorDescription: String? { NSLocalizedString("Unable to complete the calculation. Check the input size and format.", comment: "Engine error") }
}

struct DescriptorDefinition: Codable, Identifiable, Sendable {
    let id: String
    let title: String
    let group: String
    let unit: String
    let digits: Int
    let origin: String
    let verification: String
    // A trained RDKit model must not silently ingest the stable alternatives
    // as the legacy numerical features with identical CSV column names.
    var csvColumnID: String {
        ["Ipc", "AvgIpc"].contains(id) ? id + "_stable" : id
    }
}
private struct DescriptorCatalog: Decodable { let schema: Int; let descriptors: [DescriptorDefinition] }

enum Chemistry {
    static let descriptorCatalog: [DescriptorDefinition] = {
        let count = ds_workbench_catalog(nil, 0)
        guard count > 0, count <= 1_048_576 else { return [] }
        var bytes = [UInt8](repeating: 0, count: Int(count))
        let written = bytes.withUnsafeMutableBufferPointer { ds_workbench_catalog($0.baseAddress, $0.count) }
        guard written == count,
              let catalog = try? JSONDecoder().decode(DescriptorCatalog.self, from: Data(bytes)), catalog.schema == 1,
              Set(catalog.descriptors.map(\.id)).count == catalog.descriptors.count else { return [] }
        return catalog.descriptors
    }()

    static let revision: String = {
        var bytes = [UInt8](repeating: 0, count: 256)
        let count = bytes.withUnsafeMutableBufferPointer { ds_workbench_revision($0.baseAddress, $0.count) }
        precondition(count > 0 && count <= bytes.count, "Unsupported native cache revision")
        return String(decoding: bytes.prefix(Int(count)), as: UTF8.self)
    }()
    static func analyze(_ inputs: [CompoundInput]) throws -> [AnalysisRow] {
        // Keep a blank structure distinct from an empty (zero-row) batch.
        let bytes = Array(EngineInput.newlineDelimited(inputs).utf8)
        guard !bytes.isEmpty else { return [] }
        var output = [UInt8](repeating: 0, count: max(4096, inputs.count * 4096))
        var count = bytes.withUnsafeBufferPointer { source in
            output.withUnsafeMutableBufferPointer { target in
                ds_workbench_analyze(source.baseAddress, source.count, target.baseAddress, target.count)
            }
        }
        guard count > 0, count <= 32_000_000 else { throw EngineIssue.rejected }
        if Int(count) > output.count {
            output = [UInt8](repeating: 0, count: Int(count))
            count = bytes.withUnsafeBufferPointer { source in
                output.withUnsafeMutableBufferPointer { target in
                    ds_workbench_analyze(source.baseAddress, source.count, target.baseAddress, target.count)
                }
            }
        }
        guard count > 0, count <= output.count else { throw EngineIssue.rejected }
        let rows = try JSONDecoder().decode(AnalysisEnvelope.self, from: Data(output.prefix(Int(count)))).rows
        guard !descriptorCatalog.isEmpty,
              rows.allSatisfy({ !$0.isValid || $0.descriptor_values?.count == descriptorCatalog.count }) else { throw EngineIssue.rejected }
        return rows
    }

    static func graph(_ smiles: String) -> MoleculeGraph? {
        let bytes = Array(smiles.utf8)
        var output = [UInt8](repeating: 0, count: 262144)
        let count = bytes.withUnsafeBufferPointer { source in
            output.withUnsafeMutableBufferPointer { target in
                ds_workbench_drawing_graph(source.baseAddress, source.count, target.baseAddress, target.count)
            }
        }
        guard count > 0, count <= output.count else { return nil }
        guard let graph=try? JSONDecoder().decode(MoleculeGraph.self,from:Data(output.prefix(Int(count)))),
              graph.schema=="ducksmiles-workbench-drawing-v1" else { return nil }
        return graph
    }
}

typealias MoleculeGraph = MolecularDrawingGraph

@MainActor final class Workbench: ObservableObject {
    @Published var inputs: [CompoundInput] = []
    @Published var results: [AnalysisRow] = []
    private var resultsByID: [Int: AnalysisRow] = [:]
    @Published var busy = false
    @Published private(set) var exporting = false
    @Published var error: String?
    @Published var elapsed: Double?
    @Published var query = ""
    @Published var numericFilter: NumericDescriptorFilter?
    @Published var reusedCount = 0
    private let chemistry = ChemistryStore.shared
    private var revision = UUID()
    private var task: Task<Void, Never>?
    private var exportTask: Task<URL, Error>?
    private var exportedURLs: Set<URL> = []
    private let testing = ProcessInfo.processInfo.arguments.contains("--ui-testing")
    private var storageURL: URL { URL.documentsDirectory.appending(path: "chemob-inputs.json") }

    var visible: [CompoundInput] {
        guard !query.isEmpty || numericFilter != nil else { return inputs }
        let descriptorIndex = numericFilter.flatMap { filter in Chemistry.descriptorCatalog.firstIndex { $0.id == filter.descriptorID } }
        return inputs.filter { input in
            guard query.isEmpty || input.name.localizedCaseInsensitiveContains(query) || input.smiles.localizedCaseInsensitiveContains(query) else { return false }
            guard let filter = numericFilter else { return true }
            guard let descriptorIndex, let row = resultsByID[input.id], row.isValid,
                  let values = row.descriptor_values, values.indices.contains(descriptorIndex) else { return false }
            return filter.contains(values[descriptorIndex])
        }
    }
    func result(for input: CompoundInput) -> AnalysisRow? {
        resultsByID[input.id]
    }
    func importText(_ text: String, format: String = "smi", replacing: Bool = false) {
        cancel(); let token = revision; busy = true
        task = Task {
            do {
                let values = try await Task.detached { try CompoundImport.parse(text, format: format) }.value
                try Task.checkCancellation()
                guard revision == token else { return }
                // Parse and validate before changing the current library.
                let next: [CompoundInput]
                if replacing {
                    guard !values.isEmpty else { throw ImportIssue.empty }
                    next = values
                } else { next = try CompoundImport.appending(values, to: inputs) }
                if replacing { query = ""; numericFilter = nil }
                run(next)
            } catch is CancellationError { }
            catch { guard revision == token else { return }; busy = false; self.error = error.localizedDescription }
        }
    }
    func examples() { importText("CCO Ethanol\nCC(=O)Oc1ccccc1C(=O)O Aspirin\nCn1c(=O)c2c(ncn2C)n(C)c1=O Caffeine\nc1ccccc1 Benzene\n[Na+].[Cl-] Sodium chloride\ninvalid Review this input") }
    func restore() {
        guard !testing else { return }
        if let data = try? Data(contentsOf: storageURL), let values = try? JSONDecoder().decode([CompoundInput].self, from: data) { run(values) }
    }
    func cancel() { task?.cancel(); cancelExport(); revision = UUID(); busy = false }
    func cancelExport() { exportTask?.cancel() }

    func discardCSV(_ url: URL) {
        guard exportedURLs.remove(url) != nil else { return }
        Self.removeCSV(url)
    }
    private func run(_ values: [CompoundInput]) {
        cancel()
        inputs = values; results = []; resultsByID = [:]; error = nil; elapsed = nil; reusedCount = 0
        guard !values.isEmpty else { return }
        busy = true
        let token = revision
        task = Task {
            let start = ContinuousClock.now
            do {
                var collected: [AnalysisRow] = []
                // Small chunks make cancellation responsive between native calls.
                for startIndex in stride(from: 0, to: values.count, by: 64) {
                    try Task.checkCancellation()
                    let batch = Array(values[startIndex..<min(values.count, startIndex + 64)])
                    let computed = try await chemistry.analyze(batch)
                    let local = computed.rows
                    try Task.checkCancellation()
                    guard revision == token else { return }
                    reusedCount += computed.reused
                    for var row in local {
                        row.index = batch[row.index].id
                        collected.append(row)
                        resultsByID[row.index] = row
                    }
                    results = collected
                }
                guard revision == token else { return }
                await chemistry.save()
                guard revision == token else { return }
                busy = false
                let duration = start.duration(to: .now)
                elapsed = Double(duration.components.seconds) + Double(duration.components.attoseconds) / 1e18
                if !testing { try JSONEncoder().encode(values).write(to: storageURL, options: .atomic) }
            } catch is CancellationError { }
            catch { guard revision == token else { return }; busy = false; self.error = error.localizedDescription }
        }
    }
    func exportCSV() async throws -> URL {
        if let existing = exportTask {
            let url = try await existing.value
            try Task.checkCancellation()
            guard !existing.isCancelled else { throw CancellationError() }
            exportedURLs.insert(url)
            return url
        }
        let inputSnapshot = inputs
        let resultSnapshot = resultsByID
        let catalogSnapshot = Chemistry.descriptorCatalog
        let worker = Task.detached(priority: .userInitiated) {
            try Self.writeCSV(inputs: inputSnapshot, results: resultSnapshot, catalog: catalogSnapshot)
        }
        exportTask = worker
        exporting = true
        defer { exportTask = nil; exporting = false }
        let url = try await withTaskCancellationHandler {
            try await worker.value
        } onCancel: {
            worker.cancel()
        }
        guard !worker.isCancelled && !Task.isCancelled else {
            Self.removeCSV(url)
            throw CancellationError()
        }
        exportedURLs.insert(url)
        return url
    }

    nonisolated private static func removeCSV(_ url: URL) {
        // Only URLs returned by our own writer reach this function.
        try? FileManager.default.removeItem(at: url.deletingLastPathComponent())
    }

    nonisolated private static func writeCSV(inputs: [CompoundInput], results: [Int: AnalysisRow], catalog: [DescriptorDefinition]) throws -> URL {
        func quote(_ s: String) -> String { "\"" + s.replacingOccurrences(of: "\"", with: "\"\"") + "\"" }
        let added = Array(catalog.dropFirst(12))
        var lines = ["name,smiles,status,formula,mw,exact_mass,logp,tpsa,hbd,hba,heavy_atoms,bonds,formal_charge,fragments,rotatable_bonds,fraction_csp3" + added.map { "," + $0.csvColumnID }.joined()]
        for input in inputs {
            try Task.checkCancellation()
            let row = results[input.id]
            var numbers = [row?.mw, row?.exact_mass, row?.logp, row?.tpsa, row?.hbd.map(Double.init), row?.hba.map(Double.init), row?.heavy_atoms.map(Double.init), row?.bonds.map(Double.init), row?.formal_charge.map(Double.init), row?.fragments.map(Double.init), row?.rotatable_bonds.map(Double.init), row?.fraction_csp3]
            numbers += added.indices.map { index in
                guard let values = row?.descriptor_values, values.indices.contains(index + 12) else { return nil }
                return values[index + 12]
            }
            lines.append(([quote(input.name), quote(input.smiles), quote(row?.status ?? "not_computed"), quote(row?.formula ?? "")] + numbers.map { $0.map { String($0) } ?? "" }).joined(separator: ","))
        }
        try Task.checkCancellation()
        let directory = URL.temporaryDirectory.appending(path: "Chemob-export-\(UUID().uuidString)", directoryHint: .isDirectory)
        try FileManager.default.createDirectory(at: directory, withIntermediateDirectories: true)
        let url = directory.appending(path: "Chemob-results.csv")
        do {
            try lines.joined(separator: "\n").write(to: url, atomically: true, encoding: .utf8)
            try Task.checkCancellation()
            return url
        } catch {
            removeCSV(url)
            throw error
        }
    }
}
