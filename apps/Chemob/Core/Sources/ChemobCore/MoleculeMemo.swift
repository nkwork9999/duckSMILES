import Foundation

/// Bounded exact-input cache. Chemical equivalence is deliberately not inferred.
public struct MoleculeMemo<Value: Codable & Sendable>: Sendable {
    private struct Entry: Codable, Sendable {
        let smiles: String
        let value: Value
        var accessed: UInt64
    }
    private struct Snapshot: Codable {
        let namespace: String
        let entries: [Entry]
    }
    private let namespace: String
    private let capacity: Int
    private var entries: [String: Entry] = [:]
    private var clock: UInt64 = 0

    public init(namespace: String, capacity: Int) {
        self.namespace = namespace
        self.capacity = min(max(capacity, 1), 10_000)
    }
    public mutating func value(for smiles: String) -> Value? {
        guard var entry = entries[smiles] else { return nil }
        clock &+= 1; entry.accessed = clock; entries[smiles] = entry
        return entry.value
    }
    public mutating func insert(_ value: Value, for smiles: String) {
        clock &+= 1
        entries[smiles] = Entry(smiles: smiles, value: value, accessed: clock)
        if entries.count > capacity, let oldest = entries.min(by: { $0.value.accessed < $1.value.accessed })?.key {
            entries.removeValue(forKey: oldest)
        }
    }
    public func snapshot() throws -> Data {
        try JSONEncoder().encode(Snapshot(namespace: namespace, entries: entries.values.sorted { $0.accessed < $1.accessed }))
    }
    public mutating func restore(_ data: Data) throws {
        guard data.count <= 32 * 1_048_576 else { throw CocoaError(.fileReadTooLarge) }
        let snapshot = try JSONDecoder().decode(Snapshot.self, from: data)
        guard snapshot.namespace == namespace else { return }
        guard snapshot.entries.count <= 10_000,
              Set(snapshot.entries.map(\.smiles)).count == snapshot.entries.count else {
            throw CocoaError(.fileReadCorruptFile)
        }
        let retained = snapshot.entries.sorted { $0.accessed < $1.accessed }.suffix(capacity)
        entries = [:]; clock = 0
        for entry in retained { insert(entry.value, for: entry.smiles) }
    }
}
