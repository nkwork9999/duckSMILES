import XCTest
@testable import ChemobCore

final class MoleculeMemoTests: XCTestCase {
    func testReusesExactInputsAndRetainsChargeIsotopeAndStereoDifferences() {
        var memo = MoleculeMemo<Int>(namespace: "engine-v1", capacity: 10)
        memo.insert(42, for: "CCO")
        XCTAssertEqual(memo.value(for: "CCO"), 42)
        for other in ["OCC", "[13CH3]CO", "CC[O-]", "N[C@H](C)O", "N[C@@H](C)O"] {
            XCTAssertNil(memo.value(for: other))
        }
    }
    func testDiskRoundtripRejectsStaleEngineResultsAndEvictsLeastRecentEntry() throws {
        var memo = MoleculeMemo<Int>(namespace: "engine-v1", capacity: 2)
        memo.insert(1, for: "C"); memo.insert(2, for: "CC")
        XCTAssertEqual(memo.value(for: "C"), 1)
        memo.insert(3, for: "CCC")
        XCTAssertNil(memo.value(for: "CC"))
        let data = try memo.snapshot()
        var restored = MoleculeMemo<Int>(namespace: "engine-v1", capacity: 2)
        try restored.restore(data)
        XCTAssertEqual(restored.value(for: "C"), 1)
        XCTAssertEqual(restored.value(for: "CCC"), 3)
        var updated = MoleculeMemo<Int>(namespace: "engine-v2", capacity: 2)
        try updated.restore(data)
        XCTAssertNil(updated.value(for: "C"))
        XCTAssertNil(updated.value(for: "CCC"))
        XCTAssertThrowsError(try updated.restore(Data("broken".utf8)))
    }
}
