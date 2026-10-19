/// Encodes one source row per native ABI line without changing the source data.
public enum EngineInput {
    public static func newlineDelimited(_ inputs: [CompoundInput]) -> String {
        inputs.map { input in
            // Validated CX blocks may contain spaces; newlines must never add
            // ABI rows. A blank sentinel keeps the source row as an issue.
            input.smiles.isEmpty || input.smiles.utf8.count > 4096 || input.smiles.contains(where: \.isNewline) ? " " : input.smiles
        }.joined(separator: "\n")
    }
}
