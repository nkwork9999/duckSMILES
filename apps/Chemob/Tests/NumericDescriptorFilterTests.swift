import XCTest
import ChemobCore

final class NumericDescriptorFilterTests: XCTestCase {
    func testInclusiveBoundsKeepZeroAndExcludeUndefinedOrNonfiniteValues() throws {
        let range = try NumericDescriptorFilter(descriptorID: "SPS", minimum: 0, maximum: 10)
        XCTAssertTrue(range.contains(0)); XCTAssertTrue(range.contains(-0.0))
        XCTAssertTrue(range.contains(10))
        XCTAssertFalse(range.contains(-Double.leastNonzeroMagnitude))
        XCTAssertFalse(range.contains(Double(10).nextUp))
        let undefined: [Double?] = [nil, .nan, .infinity, -.infinity]
        for value in undefined { XCTAssertFalse(range.contains(value)) }
    }

    func testOneSidedRangesAndInvalidRangeDefinitions() throws {
        XCTAssertTrue(try NumericDescriptorFilter(descriptorID: "MolWt", minimum: 100, maximum: nil).contains(100))
        XCTAssertTrue(try NumericDescriptorFilter(descriptorID: "MolLogP", minimum: nil, maximum: -1).contains(-2))
        XCTAssertThrowsError(try NumericDescriptorFilter(descriptorID: "SPS", minimum: nil, maximum: nil))
        XCTAssertThrowsError(try NumericDescriptorFilter(descriptorID: "SPS", minimum: 2, maximum: 1))
        XCTAssertThrowsError(try NumericDescriptorFilter(descriptorID: "SPS", minimum: .nan, maximum: 1))
        XCTAssertThrowsError(try NumericDescriptorFilter(descriptorID: "SPS", minimum: 1, maximum: .infinity))
        XCTAssertThrowsError(try NumericDescriptorFilter(descriptorID: " ", minimum: 0, maximum: nil))
    }

    func testLocaleDecimalsAndScientificNotationWithoutGrouping() throws {
        for language in ["en_US", "hi_IN", "de_CH", "fr_CH", "it_CH"] {
            XCTAssertEqual(try NumericDescriptorFilter.parseBound(" -1.25e2 ", locale: Locale(identifier: language)), -125)
        }
        for language in ["fr_FR", "de_DE", "id_ID"] {
            XCTAssertEqual(try NumericDescriptorFilter.parseBound("-1,25e2", locale: Locale(identifier: language)), -125)
        }
        XCTAssertNil(try NumericDescriptorFilter.parseBound(" \n ", locale: Locale(identifier: "en_US")))
    }

    func testInvalidNumbersCannotBeSilentlyPartiallyParsedOrGrouped() {
        for text in ["1,234.5", "1,25", "1 234", "NaN", "+inf", "--1", "1x", "0x10", "1e10000"] {
            XCTAssertThrowsError(try NumericDescriptorFilter.parseBound(text, locale: Locale(identifier: "en_US")), text)
        }
        XCTAssertThrowsError(try NumericDescriptorFilter.parseBound("1.234,5", locale: Locale(identifier: "de_DE")))
    }
}
