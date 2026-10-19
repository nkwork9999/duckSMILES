import XCTest
import ChemobCore

final class DescriptorFormattingTests: XCTestCase {
    func testRoundedSmallNegativeValuesDisplayAsZero() {
        let english = Locale(identifier: "en_US")
        let german = Locale(identifier: "de_DE")
        XCTAssertEqual(DescriptorFormatting.text(-0.0014, digits: 2, locale: english), "0.00")
        XCTAssertEqual(DescriptorFormatting.text(-0.0014, digits: 2, locale: german), "0,00")
        XCTAssertEqual(DescriptorFormatting.text(-0.00049, digits: 3, locale: english), "0.000")
        XCTAssertEqual(DescriptorFormatting.text(-0.4, digits: 0, locale: english), "0")
        XCTAssertEqual(DescriptorFormatting.text(-0.0051, digits: 2, locale: english), "-0.01")
        XCTAssertEqual(DescriptorFormatting.text(-0.000001, digits: 6, locale: english), "-0.000001")
        XCTAssertEqual(Double(DescriptorFormatting.text(-1e9, digits: 4, locale: english)), -1e9)
    }
    func testLargeInformationValueKeepsMagnitudeAndFitsPropertyCell() {
        let text = DescriptorFormatting.text(1.23456789e53, digits: 4, locale: Locale(identifier: "en_US"))
        XCTAssertLessThanOrEqual(text.count, 16)
        let represented = Double(text)
        XCTAssertNotNil(represented)
        if let represented {
            XCTAssertEqual(represented / 1e53, 1.2346, accuracy: 1e-12)
        }
    }

    func testDisplayedDecimalUsesTheInterfaceLocale() {
        XCTAssertEqual(DescriptorFormatting.text(1.25, digits: 2, locale: Locale(identifier: "de_CH")), "1.25")
        XCTAssertEqual(DescriptorFormatting.text(1.25, digits: 2, locale: Locale(identifier: "de_DE")), "1,25")
        // Foundation's Swiss French locale uses a decimal point.
        XCTAssertEqual(DescriptorFormatting.text(1.25, digits: 2, locale: Locale(identifier: "fr_CH")), "1.25")
        XCTAssertEqual(DescriptorFormatting.text(1.25, digits: 2, locale: Locale(identifier: "fr_FR")), "1,25")
    }

    func testZeroAndUndefinedRemainDistinct() {
        let locale = Locale(identifier: "en_US")
        XCTAssertEqual(DescriptorFormatting.text(nil, digits: 2, locale: locale), "—")
        XCTAssertEqual(DescriptorFormatting.text(-0.0, digits: 2, locale: locale), "0.00")
        XCTAssertEqual(DescriptorFormatting.text(42, digits: 0, locale: locale), "42")
    }
}
