import XCTest
@testable import ChemobCore

final class CompoundImportTests: XCTestCase {
    func testCXAnnotationAndOptionalNameSurviveSMIImport() throws {
        let structure = "C[C@H](O)F \t|wU:1.0,wD:1.1|"
        let values = try CompoundImport.parse(structure + " Research sample\nCCO Ethanol\nCC(O)F |wU:1.0|")
        XCTAssertEqual(values.map(\.smiles), [structure, "CCO", "CC(O)F |wU:1.0|"])
        XCTAssertEqual(values.map(\.name), ["Research sample", "Ethanol", "Compound 3"])
        XCTAssertEqual(EngineInput.newlineDelimited(values), structure + "\nCCO\nCC(O)F |wU:1.0|")
    }
    func testCXInvalidOrUnknownMeaningIsNotDiscardedAtImport() throws {
        for source in ["CC |wU:0.0", "CC |^1:0|", "CC |w:0.0|", "CC |wU:0.0| |^1:0|"] {
            let values = try CompoundImport.parse(source)
            XCTAssertEqual(values.count, 1)
            XCTAssertEqual(values[0].smiles, source, "the engine must see the original unsupported input")
            XCTAssertEqual(values[0].name, "Compound 1")
        }
    }
    func testCXCommaListSurvivesCSVAndTSV() throws {
        let source = "CCC |wU:0.0,wD:1.1|"
        let csv = try CompoundImport.parse("name,smiles\nSample,\"" + source + "\"\n", format: "csv")
        let tsv = try CompoundImport.parse("name\tsmiles\nSample\t" + source + "\n", format: "tsv")
        XCTAssertEqual(csv, tsv)
        XCTAssertEqual(csv[0].smiles, source)
    }
    func testMultilineStructureCannotInjectExtraEngineRows() throws {
        let values = try CompoundImport.parse("name,smiles\nBroken,\"CCO\nC\"\nGood,CCO\nEmpty,", format: "csv")
        XCTAssertEqual(values.count, 3)
        XCTAssertEqual(EngineInput.newlineDelimited(values), " \nCCO\n ")
        XCTAssertEqual(values[0].smiles, "CCO\nC", "original input stays available for correction")
        XCTAssertEqual(EngineInput.newlineDelimited([]), "")
        for newline in ["\r", "\r\n", "\u{2028}"] {
            XCTAssertEqual(EngineInput.newlineDelimited([CompoundInput(id: 9, name: "Broken", smiles: "CCO" + newline + "C")]), " ")
        }
    }
    func testTenThousandResultRowsCanBeReimported() throws {
        let header = "name,smiles,status,formula,mw,exact_mass,logp,tpsa,hbd,hba,heavy_atoms,bonds,formal_charge,fragments,rotatable_bonds,fraction_csp3\n"
        let row = "\"試料, \"\"quoted\"\"\nsecond line\",CC(=O)Oc1ccccc1C(=O)O,ok,C9H8O4,180.15899999999996,180.042258736,1.3100999999999998,63.60000000000001,1.0,3.0,13.0,13.0,0.0,1.0,2.0,0.1111111111111111\n"
        let text = header + String(repeating: row, count: 10_000)
        XCTAssertGreaterThan(text.utf8.count, 1_048_576)
        let values = try CompoundImport.parse(text, format: "csv")
        XCTAssertEqual(values.count, 10_000)
        XCTAssertEqual(values.first?.name, "試料, \"quoted\"\nsecond line")
        XCTAssertEqual(values.last?.smiles, "CC(=O)Oc1ccccc1C(=O)O")
        XCTAssertEqual(values.last?.id, 9999)
    }
    func testWideDescriptorCSVPreservesTenThousandSourceRows() throws {
        // A full descriptor export can exceed the old 8 MiB limit even though
        // the molecule count is within the promised 10,000-row contract.
        let columns = (0..<267).map { "property_\($0)" }.joined(separator: ",")
        let numeric = Array(repeating: "0.123456789", count: 267).joined(separator: ",")
        let name = "試料 α नमूना\n\"quoted\""
        let text = "name,smiles," + columns + "\n" + String(repeating: "\"試料 α नमूना\n\"\"quoted\"\"\",CCO," + numeric + "\n", count: 10_000)
        XCTAssertGreaterThan(text.utf8.count, 8_388_608)
        let inputs = try CompoundImport.parse(text, format: "csv")
        XCTAssertEqual(inputs.count, 10_000)
        XCTAssertEqual(inputs.first?.name, name)
        XCTAssertEqual(inputs.last?.name, name)
        XCTAssertEqual(inputs.last?.smiles, "CCO")
        XCTAssertEqual(inputs.last?.id, 9999)
    }
    func testOversizeStructureDoesNotRejectOtherRowsInItsBatch() {
        let inputs = [CompoundInput(id: 0, name: "Too long", smiles: String(repeating: "C", count: 4097)), CompoundInput(id: 1, name: "Ethanol", smiles: "CCO")]
        XCTAssertEqual(EngineInput.newlineDelimited(inputs), " \nCCO")
        XCTAssertEqual(inputs[0].smiles.utf8.count, 4097)
    }
    func testOversizeFilesAndExcessRowsRemainBounded() {
        XCTAssertThrowsError(try CompoundImport.parse(String(repeating: "x", count: 33_554_433)))
        XCTAssertThrowsError(try CompoundImport.parse(String(repeating: "CCO\n", count: 10_001)))
    }
    func testAddingCompoundsPreservesExistingRowsAndAssignsUniqueIDs() throws {
        let existing = [CompoundInput(id: 0, name: "Ethanol", smiles: "CCO"), CompoundInput(id: 3, name: "Salt", smiles: "[Na+].[Cl-]")]
        let incoming = [CompoundInput(id: 0, name: "Methane", smiles: "C")]
        let result = try CompoundImport.appending(incoming, to: existing)
        XCTAssertEqual(result.map(\.smiles), ["CCO", "[Na+].[Cl-]", "C"])
        XCTAssertEqual(result.map(\.id), [0, 3, 4])
    }
    func testPastedSMILESPreserveInvalidRowsAndNames() throws {
        let values = try CompoundImport.parse("CCO Ethanol\ninvalid Broken\n\nc1ccccc1 Benzene")
        XCTAssertEqual(values.map(\.smiles), ["CCO", "invalid", "c1ccccc1"])
        XCTAssertEqual(values.map(\.name), ["Ethanol", "Broken", "Benzene"])
        XCTAssertEqual(Set(values.map(\.id)).count, 3)
    }
    func testCSVRespectsQuotedCommaAndMultilineFields() throws {
        let values = try CompoundImport.parse("name,SMILES\r\n\"Salt, sample\",[Na+].[Cl-]\r\n\"two\nlines\",CCO", format: "csv")
        XCTAssertEqual(values.map(\.smiles), ["[Na+].[Cl-]", "CCO"])
        XCTAssertEqual(values.first?.name, "Salt, sample")
        XCTAssertEqual(values.last?.name, "two\nlines")
    }
    func testMissingStructureColumnAndUnclosedQuoteAreReported() {
        XCTAssertThrowsError(try CompoundImport.parse("name,mass\nWater,18", format: "csv"))
        XCTAssertThrowsError(try CompoundImport.parse("name,smiles\n\"unclosed,CCO", format: "csv"))
    }
}
