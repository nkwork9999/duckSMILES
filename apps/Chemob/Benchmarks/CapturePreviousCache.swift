// Run against the frozen v38 library, then retain its actual Store snapshot.
import Foundation
import ChemobCore

@main struct CapturePreviousCache {
    static func main() async throws {
        precondition(CommandLine.arguments.count == 2)
        precondition(Chemistry.revision == "ducksmiles-workbench-v38-cxwedge291-rdkit2026036")
        let directory = URL(fileURLWithPath: CommandLine.arguments[1], isDirectory: true)
        precondition(!FileManager.default.fileExists(atPath: directory.path))
        let inputs = ["F[Pt@SP1](Cl)(Br)I", String(repeating: "C", count: 257)].enumerated().map {
            CompoundInput(id: $0.offset, name: "Previous unsupported", smiles: $0.element)
        }
        let store = ChemistryStore(directory: directory)
        let result = try await store.analyze(inputs)
        precondition(result.reused == 0 && result.rows.count == 2 && result.rows.allSatisfy { !$0.isValid })
        let counts = await store.statistics
        precondition(counts.computedStructures == 2)
        await store.save()
        let reopened = ChemistryStore(directory: directory)
        let restored = try await reopened.analyze(inputs)
        precondition(restored.reused == 2 && restored.rows.allSatisfy { !$0.isValid })
        print("Captured and reopened two real v38 unsupported results")
    }
}
