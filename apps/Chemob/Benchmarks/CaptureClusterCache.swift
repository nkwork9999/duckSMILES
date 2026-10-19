// Capture real frozen-v41 unsupported metal-cluster results for migration tests.
import Foundation
import ChemobCore

@main struct CaptureClusterCache {
    struct Source: Decodable { let name: String; let smiles: String }
    static func main() async throws {
        precondition(CommandLine.arguments.count == 3)
        precondition(Chemistry.revision == "ducksmiles-workbench-v41-harary291-rdkit2026036")
        let source = URL(fileURLWithPath: CommandLine.arguments[1])
        let directory = URL(fileURLWithPath: CommandLine.arguments[2], isDirectory: true)
        precondition(!FileManager.default.fileExists(atPath: directory.path))
        let sources = try JSONDecoder().decode([Source].self, from: Data(contentsOf: source))
        precondition(sources.count == 5)
        let inputs = sources.enumerated().map { CompoundInput(id: $0.offset, name: $0.element.name, smiles: $0.element.smiles) }
        let store = ChemistryStore(directory: directory)
        let result = try await store.analyze(inputs)
        precondition(result.reused == 0 && result.rows.count == 5 && result.rows.allSatisfy { !$0.isValid })
        let counts = await store.statistics
        precondition(counts.computedStructures == 5)
        await store.save()
        let reopened = ChemistryStore(directory: directory)
        let restored = try await reopened.analyze(inputs)
        precondition(restored.reused == 5 && restored.rows.allSatisfy { !$0.isValid })
        print("Captured and reopened five real v41 unsupported metal-cluster results")
    }
}
