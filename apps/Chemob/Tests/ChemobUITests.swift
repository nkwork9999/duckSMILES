import XCTest
import UIKit

final class ChemobUITests: XCTestCase {
    @MainActor func testDoubleBondGeometryAndBridgeCarbonHaveVisibleResults() throws {
        continueAfterFailure=false
        XCUIDevice.shared.orientation = .portrait
        let app=XCUIApplication()
        app.launchArguments=["--ui-testing","-AppleLanguages","(en)","-AppleLocale","en_US"]
        app.launch();app.buttons["pasteSMILES"].tap()
        let editor=app.textViews["smilesEditor"]
        XCTAssertTrue(editor.waitForExistence(timeout:10))
        editor.tap();editor.typeText("F/C=C/F Trans alkene\nF/C=C\\F Cis alkene\nC/C1=C(/C)CCCCCCC1 Ring geometry limit\nC1CC2CCC1C2 Norbornane")
        app.buttons["analyzeInput"].tap()
        XCTAssertTrue(app.buttons["compound-0"].waitForExistence(timeout:30))
        let pad=app.windows.firstMatch.frame.width>=700
        for i in 0..<4 {
            let row=app.buttons["compound-\(i)"]
            if pad && !row.exists {
                let next=app.buttons["Next page"]
                XCTAssertTrue(next.isHittable,"The selected input must be reachable through pagination")
                next.tap()
            }
            XCTAssertTrue(row.waitForExistence(timeout:10))
            row.tap()
            XCTAssertTrue(app.otherElements["chemicalDepiction"].waitForExistence(timeout:10))
            if i==2 { XCTAssertTrue(app.staticTexts["2D structure · some stereochemistry not shown"].exists) }
            else { XCTAssertTrue(app.staticTexts["2D structure"].exists) }
            if i==3 { XCTAssertTrue(app.staticTexts["C7H12"].firstMatch.exists) }
            let shot=XCTAttachment(screenshot:app.screenshot());shot.name="E Z and bridge example \(i)";shot.lifetime = .keepAlways;add(shot)
            if !pad { app.buttons["Done"].tap() }
        }
    }

    @MainActor func testChemical2DPreviewAndStereoNoticeInSelectedResults() throws {
        continueAfterFailure = false
        XCUIDevice.shared.orientation = .portrait
        let app=XCUIApplication()
        app.launchArguments=["--ui-testing","-AppleLanguages","(en)","-AppleLocale","en_US"]
        app.launch()
        app.buttons["pasteSMILES"].tap()
        let editor=app.textViews["smilesEditor"]
        XCTAssertTrue(editor.waitForExistence(timeout: 10))
        editor.tap(); editor.typeText("c1ccccc1 Benzene\nCC(=O)Oc1ccccc1C(=O)O Aspirin\nCn1c(=O)c2c(ncn2C)n(C)c1=O Caffeine\nN[C@@H](C)C(=O)O Alanine\nF/C=C/F Difluoroethene")
        app.buttons["analyzeInput"].tap()
        XCTAssertTrue(app.buttons["compound-0"].waitForExistence(timeout: 30))
        let pad=app.windows.firstMatch.frame.width >= 700
        for i in 0..<5 {
            let row=app.buttons["compound-\(i)"]
            if pad && !row.exists {
                let next=app.buttons["Next page"]
                XCTAssertTrue(next.isHittable,"The selected input must be reachable through pagination")
                next.tap()
            }
            XCTAssertTrue(row.waitForExistence(timeout:10))
            row.tap()
            XCTAssertTrue(app.otherElements["chemicalDepiction"].waitForExistence(timeout: 10))
            XCTAssertTrue(app.staticTexts["2D structure"].exists)
            if pad { XCTAssertTrue(app.textFields["compoundSearch"].isHittable) }
            else {
                // Values and their units must remain separate when the 2D preview loads.
                let weight=app.staticTexts["descriptor-MolWt"].frame
                let exact=app.staticTexts["descriptor-ExactMolWt"].frame
                let unit=app.staticTexts["g/mol"].firstMatch.frame
                XCTAssertFalse(weight.insetBy(dx:-2,dy:0).intersects(exact),"Property numbers must not overlap")
                XCTAssertFalse(unit.insetBy(dx:-2,dy:0).intersects(exact),"A property unit must not overlap its neighboring value")
            }
            let shot=XCTAttachment(screenshot: app.screenshot());shot.name="Chemical 2D example \(i)";shot.lifetime = .keepAlways;add(shot)
            if !pad { app.buttons["Done"].tap() }
        }
    }

    @MainActor func testPhoneDescriptorLabelsAndValuesRemainInTheirColumns() throws {
        continueAfterFailure = false
        let app = XCUIApplication()
        app.launchArguments = ["--ui-testing", "-AppleLanguages", "(en)", "-AppleLocale", "en_US"]
        app.launch()
        guard app.windows.firstMatch.frame.width < 700 else { throw XCTSkip("Phone property columns") }
        app.buttons["loadExamples"].tap()
        XCTAssertTrue(app.buttons["compound-0"].waitForExistence(timeout: 30))
        app.buttons["compound-0"].tap()
        XCTAssertTrue(app.staticTexts["descriptor-ExactMolWt"].waitForExistence(timeout: 10))
        let columnX = app.staticTexts["descriptor-ExactMolWt"].frame.minX
        print("Initial phone right column frame: \(app.staticTexts["descriptor-ExactMolWt"].frame)")
        let initial = XCTAttachment(screenshot: app.screenshot())
        initial.name = "Phone property columns after presentation"; initial.lifetime = .keepAlways; add(initial)
        let bonds = app.staticTexts["descriptor-ds_BondCount"]
        // The 2D preview places these rows below the compact phone viewport.
        // Keep the same column/separation assertions after reaching the actual rows.
        for _ in 0..<4 {
            if bonds.isHittable && app.staticTexts["Fragments"].isHittable { break }
            app.scrollViews.firstMatch.swipeUp()
        }
        XCTAssertTrue(bonds.waitForExistence(timeout: 10))
        XCTAssertTrue(bonds.isHittable)
        XCTAssertTrue(app.staticTexts["Fragments"].isHittable)
        for (left, right) in [("Heavy atoms", "Bonds"), ("Formal charge", "Fragments")] {
            let first = app.staticTexts[left], second = app.staticTexts[right]
            XCTAssertTrue(first.exists && second.exists)
            print("Phone label frames \(left): \(first.frame); \(right): \(second.frame)")
            XCTAssertLessThanOrEqual(first.frame.maxX + 8, second.frame.minX, "Property labels must stay separate")
        }
        XCTAssertEqual(bonds.frame.minX, columnX, accuracy: 2, "Integer properties must align with the right column")
        let lower = XCTAttachment(screenshot: app.screenshot())
        lower.name = "Phone lower property labels and values remain separate"; lower.lifetime = .keepAlways; add(lower)
    }

    @MainActor func testMetalClustersReachInlineResultsWithOriginalGeometryPreserved() throws {
        continueAfterFailure = false
        XCUIDevice.shared.orientation = .portrait
        let app = XCUIApplication()
        app.launchArguments = ["--ui-testing", "-AppleLanguages", "(en)", "-AppleLocale", "en_US"]
        app.launch()
        guard app.windows.firstMatch.frame.width >= 700 else { throw XCTSkip("iPad inline results") }
        app.buttons["pasteSMILES"].tap()
        let editor = app.textViews["smilesEditor"]
        XCTAssertTrue(editor.waitForExistence(timeout: 10))
        let iron = "[S]1[Fe]2[S-][Fe]3[S-][Fe]1[S]23"
        let geometry = "[OH][Fe@SP1]([OH])([OH])[O]1[O][Fe@TB14]1([OH])([OH])[O][Fe@TB8]([OH])([OH])([OH])[OH]"
        let gold = "[Cl][Au]1([Cl])[Cl][Au]([Cl])([Cl])[Cl]1"
        editor.tap(); editor.typeText(iron + " Iron cluster\n" + geometry + " Original geometry\n" + gold + " Gold cluster\ninvalid Review source")
        app.buttons["analyzeInput"].tap()
        XCTAssertTrue(app.buttons["compound-0"].waitForExistence(timeout: 30))
        app.buttons["compound-0"].tap()
        XCTAssertTrue(app.staticTexts["Fe3S4-2"].exists)
        XCTAssertTrue(app.staticTexts[iron].firstMatch.exists)
        app.buttons["descriptorGroup"].tap(); app.buttons["Stereo"].tap()
        XCTAssertEqual(app.staticTexts["descriptor-SPS"].label, "40.2857")
        XCTAssertTrue(app.staticTexts["descriptor-SPS"].isHittable)
        XCTAssertTrue(app.textFields["compoundSearch"].isHittable)
        let first = XCTAttachment(screenshot: app.screenshot())
        first.name = "Metal cluster native results in fixed inline panel"
        first.lifetime = .keepAlways; add(first)
        app.buttons["compound-1"].tap()
        XCTAssertTrue(app.staticTexts["H9Fe3O12"].exists)
        XCTAssertTrue(app.staticTexts[geometry].firstMatch.exists)
        XCTAssertEqual(app.staticTexts["descriptor-SPS"].label, "22.8000")
        let second = XCTAttachment(screenshot: app.screenshot())
        second.name = "Original non-tetrahedral source retained with corrected calculation"
        second.lifetime = .keepAlways; add(second)
        app.segmentedControls["compoundCategory"].buttons["Valid"].tap()
        XCTAssertTrue(app.staticTexts["3 matches"].exists)
        XCTAssertFalse(app.buttons["compound-3"].exists)
    }

    @MainActor func testCoordinateFreeCXWedgesPreserveSourceAndReachInlineStereoValues() throws {
        continueAfterFailure = false
        XCUIDevice.shared.orientation = .portrait
        let app = XCUIApplication()
        app.launchArguments = ["--ui-testing", "-AppleLanguages", "(en)", "-AppleLocale", "en_US"]
        app.launch()
        guard app.windows.firstMatch.frame.width >= 700 else { throw XCTSkip("iPad inline results") }
        app.buttons["pasteSMILES"].tap()
        let editor = app.textViews["smilesEditor"]
        XCTAssertTrue(editor.waitForExistence(timeout: 10))
        let source = "CC(O)F |wU:1.0|"
        editor.tap(); editor.typeText(source + " Wedge source\nCC |^1:0| Unsupported radical\nCCO Ethanol")
        app.buttons["analyzeInput"].tap()
        XCTAssertTrue(app.buttons["compound-0"].waitForExistence(timeout: 30))
        app.buttons["compound-0"].tap()
        XCTAssertTrue(app.staticTexts["Wedge source"].firstMatch.exists)
        XCTAssertTrue(app.staticTexts["C2H5FO"].exists)
        XCTAssertTrue(app.staticTexts[source].firstMatch.exists)
        XCTAssertFalse(app.staticTexts["-0.00"].exists)
        XCTAssertTrue(app.staticTexts["0.00"].exists)
        app.buttons["descriptorGroup"].tap(); app.buttons["Stereo"].tap()
        XCTAssertEqual(app.staticTexts["descriptor-SPS"].label, "15.7500")
        XCTAssertEqual(app.staticTexts["descriptor-NumAtomStereoCenters"].label, "1")
        XCTAssertEqual(app.staticTexts["descriptor-NumUnspecifiedAtomStereoCenters"].label, "1")
        XCTAssertTrue(app.staticTexts["descriptor-SPS"].isHittable)
        XCTAssertTrue(app.textFields["compoundSearch"].isHittable)
        let result = XCTAttachment(screenshot: app.screenshot())
        result.name = "Original CX wedge source and validated descriptor values inline"
        result.lifetime = .keepAlways; add(result)
        app.segmentedControls["compoundCategory"].buttons["Needs review"].tap()
        XCTAssertTrue(app.buttons["compound-1"].exists)
        XCTAssertFalse(app.buttons["compound-0"].exists)
        app.buttons["compound-1"].tap()
        XCTAssertTrue(app.staticTexts["CC |^1:0|"].firstMatch.exists)
        XCTAssertTrue(app.staticTexts["Invalid or unsupported structure. The original SMILES is preserved."].exists)
    }

    @MainActor func testCSVShareSheetShowsTheFileAndCancelsWithoutChangingTheLibrary() throws {
        continueAfterFailure = false
        XCUIDevice.shared.orientation = .portrait
        let app = XCUIApplication()
        app.launchArguments = ["--ui-testing", "-AppleLanguages", "(en)", "-AppleLocale", "en_US"]
        app.launch(); app.buttons["loadExamples"].tap()
        XCTAssertTrue(app.buttons["compound-0"].waitForExistence(timeout: 30))
        app.buttons["exportCSV"].tap()
        let file = app.otherElements["LP.CaptionBar.TopCaption"].firstMatch
        guard file.waitForExistence(timeout: 30) else {
            print(app.debugDescription)
            let failed = XCTAttachment(screenshot: XCUIScreen.main.screenshot())
            failed.name = "Share presentation diagnostic"; failed.lifetime = .keepAlways; add(failed)
            XCTFail("A loaded activity controller is insufficient: the file must be visible in the actual share sheet")
            return
        }
        XCTAssertTrue(file.label.contains("Chemob-results"))
        XCTAssertTrue(file.isHittable)
        let shown = XCTAttachment(screenshot: XCUIScreen.main.screenshot())
        shown.name = "CSV file visible in actual system share sheet"; shown.lifetime = .keepAlways; add(shown)
        if app.buttons["header.closeButton"].firstMatch.exists { app.buttons["header.closeButton"].firstMatch.tap() }
        else { print(app.debugDescription); XCTFail("The share sheet must offer cancellation"); return }
        XCTAssertTrue(app.textFields["compoundSearch"].isHittable)
        XCTAssertTrue(app.staticTexts["6 matches"].exists)
        if !app.buttons["compound-5"].exists { app.buttons["Next page"].tap() }
        XCTAssertTrue(app.buttons["compound-5"].waitForExistence(timeout: 5))
    }

    @MainActor func testFrenchNumericDisplayUsesTheSameDecimalSeparatorInTableAndDetails() throws {
        continueAfterFailure = false
        XCUIDevice.shared.orientation = .portrait
        let app = XCUIApplication()
        app.launchArguments = ["--ui-testing", "-AppleLanguages", "(fr)", "-AppleLocale", "fr_FR"]
        app.launch()
        guard app.windows.firstMatch.frame.width >= 700 else { throw XCTSkip("iPad inline table") }
        app.buttons["pasteSMILES"].tap()
        let editor = app.textViews["smilesEditor"]
        XCTAssertTrue(editor.waitForExistence(timeout: 10))
        editor.tap(); editor.typeText("CCO Ethanol")
        app.buttons["analyzeInput"].tap()
        let row = app.buttons["compound-0"]
        XCTAssertTrue(row.waitForExistence(timeout: 30))
        XCTAssertTrue(row.staticTexts["46,07"].exists)
        row.tap()
        XCTAssertEqual(app.staticTexts["descriptor-MolWt"].label, "46,07")
        app.buttons["descriptorGroup"].tap(); app.buttons["Graphe"].tap()
        XCTAssertEqual(app.staticTexts["descriptor-Ipc"].label, "2,7549")
        XCTAssertEqual(app.staticTexts["descriptor-AvgIpc"].label, "0,9183")
        XCTAssertTrue(app.staticTexts["descriptor-AvgIpc"].isHittable)
        let shot = XCTAttachment(screenshot: app.screenshot())
        shot.name = "French table and graph descriptors use locale decimals"; shot.lifetime = .keepAlways; add(shot)
    }

    @MainActor func testOvervalentInputsAreReviewRowsWhileHypervalentPhosphorusRemainsValid() throws {
        continueAfterFailure = false
        XCUIDevice.shared.orientation = .portrait
        let app = XCUIApplication()
        app.launchArguments = ["--ui-testing", "-AppleLanguages", "(en)", "-AppleLocale", "en_US"]
        app.launch()
        guard app.windows.firstMatch.frame.width >= 700 else { throw XCTSkip("iPad inline categories") }
        app.buttons["pasteSMILES"].tap()
        let editor = app.textViews["smilesEditor"]
        XCTAssertTrue(editor.waitForExistence(timeout: 10))
        editor.tap(); editor.typeText("[P-](F)(F)(F)(F)(F)F Hexafluorophosphate\nC(C)(C)(C)(C)C Overvalent carbon\nCCO Ethanol")
        app.buttons["analyzeInput"].tap()
        XCTAssertTrue(app.buttons["compound-0"].waitForExistence(timeout: 30))
        app.buttons["compound-1"].tap()
        XCTAssertTrue(app.staticTexts["Invalid or unsupported structure. The original SMILES is preserved."].exists)
        XCTAssertTrue(app.staticTexts["C(C)(C)(C)(C)C"].firstMatch.exists)
        XCTAssertFalse(app.staticTexts["descriptor-MolWt"].exists)
        app.segmentedControls["compoundCategory"].buttons["Needs review"].tap()
        XCTAssertTrue(app.buttons["compound-1"].exists)
        XCTAssertFalse(app.buttons["compound-0"].exists)
        XCTAssertTrue(app.staticTexts["1 match"].exists)
        let rejected = XCTAttachment(screenshot: app.screenshot())
        rejected.name = "Overvalent source retained as a review row"; rejected.lifetime = .keepAlways; add(rejected)
        app.segmentedControls["compoundCategory"].buttons["Valid"].tap()
        XCTAssertFalse(app.buttons["compound-1"].exists)
        app.buttons["compound-0"].tap()
        XCTAssertTrue(app.staticTexts["F6P-"].exists)
        XCTAssertTrue(app.staticTexts["descriptor-MolWt"].isHittable)
        let preserved = XCTAttachment(screenshot: app.screenshot())
        preserved.name = "Valid charged hypervalent phosphorus preserved"; preserved.lifetime = .keepAlways; add(preserved)
    }

    @MainActor func testStableGraphInformationIsInlineAndExplainsModelCompatibility() throws {
        continueAfterFailure = false
        XCUIDevice.shared.orientation = .portrait
        let app = XCUIApplication()
        app.launchArguments = ["--ui-testing", "-AppleLanguages", "(en)", "-AppleLocale", "en_US"]
        app.launch()
        guard app.windows.firstMatch.frame.width >= 700 else { throw XCTSkip("iPad inline properties") }
        app.buttons["pasteSMILES"].tap()
        let editor = app.textViews["smilesEditor"]
        XCTAssertTrue(editor.waitForExistence(timeout: 10))
        editor.tap(); editor.typeText(String(repeating: "C", count: 128) + " Long carbon chain\n[HH] Hydrogen")
        app.buttons["analyzeInput"].tap()
        XCTAssertTrue(app.buttons["compound-0"].waitForExistence(timeout: 45))
        app.buttons["compound-0"].tap()
        app.buttons["descriptorGroup"].tap(); app.buttons["Graph"].tap()
        for _ in 0..<4 {
            if app.staticTexts["descriptor-Ipc"].exists { break }
            app.buttons["Next properties"].tap()
        }
        XCTAssertEqual(app.staticTexts["descriptor-Ipc"].label, "1.5501E27")
        XCTAssertEqual(app.staticTexts["descriptor-AvgIpc"].label, "3.8058")
        XCTAssertTrue(app.staticTexts["descriptor-Ipc"].isHittable)
        XCTAssertTrue(app.staticTexts["descriptor-AvgIpc"].isHittable)
        XCTAssertTrue(app.textFields["compoundSearch"].isHittable)
        XCTAssertFalse(app.buttons["Done"].exists)
        let numbers = XCTAttachment(screenshot: app.screenshot())
        numbers.name = "Stable graph information scientific notation inline"; numbers.lifetime = .keepAlways; add(numbers)
        app.buttons["graphInformationInfo"].tap()
        let explanation = app.staticTexts["graphInformationExplanation"]
        XCTAssertTrue(explanation.waitForExistence(timeout: 5))
        XCTAssertTrue(explanation.label.contains("RDKit Ipc or AvgIpc"))
        let note = XCTAttachment(screenshot: app.screenshot())
        note.name = "Explicit model feature compatibility explanation"; note.lifetime = .keepAlways; add(note)
        // A known visible row outside the popover dismisses it on iPad.
        app.buttons["compound-1"].tap()
        app.buttons["compound-1"].tap()
        XCTAssertEqual(app.staticTexts["descriptor-Ipc"].label, "0.0000")
        XCTAssertEqual(app.staticTexts["descriptor-AvgIpc"].label, "0.0000")
        XCTAssertTrue(app.staticTexts["H2"].exists)
    }

    @MainActor func testReplaceLibraryIsExplicitAndMalformedCSVKeepsPreviousInputs() throws {
        continueAfterFailure = false
        XCUIDevice.shared.orientation = .portrait
        let app = XCUIApplication()
        app.launchArguments = ["--ui-testing", "-AppleLanguages", "(en)", "-AppleLocale", "en_US"]
        app.launch(); app.buttons["loadExamples"].tap()
        XCTAssertTrue(app.buttons["compound-0"].waitForExistence(timeout: 30))
        app.buttons["pasteSMILES"].tap()
        let mode = app.segmentedControls["inputMode"]
        XCTAssertTrue(mode.waitForExistence(timeout: 10), "a full10k library needs an explicit replacement path")
        mode.buttons["Replace library"].tap()
        app.segmentedControls["inputFormat"].buttons["CSV"].tap()
        let editor = app.textViews["smilesEditor"]
        editor.tap(); editor.typeText("smiles,name\nCCO,Replacement ethanol\nCCCC,Replacement butane")
        app.buttons["analyzeInput"].tap()
        XCTAssertTrue(app.buttons["compound-1"].waitForExistence(timeout: 30))
        XCTAssertTrue(app.staticTexts["Replacement butane"].firstMatch.exists)
        XCTAssertFalse(app.buttons["compound-2"].exists)
        XCTAssertTrue(app.staticTexts["2 matches"].exists)
        app.buttons["pasteSMILES"].tap()
        app.segmentedControls["inputMode"].buttons["Replace library"].tap()
        app.segmentedControls["inputFormat"].buttons["CSV"].tap()
        let broken = app.textViews["smilesEditor"]
        broken.tap(); broken.typeText("smiles,name\n\"CCO")
        app.buttons["analyzeInput"].tap()
        XCTAssertTrue(app.alerts["Check your input"].waitForExistence(timeout: 10))
        app.alerts.buttons["OK"].tap()
        XCTAssertTrue(app.staticTexts["Replacement butane"].firstMatch.exists)
        XCTAssertTrue(app.staticTexts["2 matches"].exists)
        XCTAssertFalse(app.buttons["compound-2"].exists)
        app.buttons["pasteSMILES"].tap()
        app.segmentedControls["inputMode"].buttons["Replace library"].tap()
        app.segmentedControls["inputFormat"].buttons["CSV"].tap()
        let empty = app.textViews["smilesEditor"]
        empty.tap(); empty.typeText("smiles,name")
        app.buttons["analyzeInput"].tap()
        XCTAssertTrue(app.alerts["Check your input"].waitForExistence(timeout: 10))
        XCTAssertTrue(app.alerts.staticTexts["No compounds were found in the input."].exists)
        app.alerts.buttons["OK"].tap()
        XCTAssertTrue(app.staticTexts["Replacement butane"].firstMatch.exists)
        XCTAssertTrue(app.staticTexts["2 matches"].exists)
        let shot = XCTAttachment(screenshot: app.screenshot())
        shot.name = "Explicit replacement preserves previous list on parse failure"; shot.lifetime = .keepAlways; add(shot)
    }

    @MainActor func testNumericFiltersUseCachedDescriptorsAndExcludeUndefinedValues() throws {
        continueAfterFailure = false
        XCUIDevice.shared.orientation = .portrait
        let app = XCUIApplication()
        app.launchArguments = ["--ui-testing", "-AppleLanguages", "(en)", "-AppleLocale", "en_US"]
        app.launch()
        guard app.windows.firstMatch.frame.width >= 700 else { throw XCTSkip("iPad inline numeric filtering") }
        app.buttons["pasteSMILES"].tap()
        let editor = app.textViews["smilesEditor"]
        XCTAssertTrue(editor.waitForExistence(timeout: 10))
        editor.tap(); editor.typeText("C Methane\nCC(O)Cl Chlorohydrin\n[H][H] Hydrogen\ninvalid Broken")
        app.buttons["analyzeInput"].tap()
        XCTAssertTrue(app.buttons["compound-0"].waitForExistence(timeout: 30))
        let filter = app.buttons["numericFilter"]
        XCTAssertTrue(filter.waitForExistence(timeout: 10), "descriptor ranges must be available below search")
        filter.tap()
        app.buttons["filterDescriptor"].tap(); app.buttons["Stereo"].tap(); app.buttons["Spacial score (normalized)"].tap()
        let lower = app.textFields["filterMinimum"]
        lower.tap(); lower.typeText("0")
        app.buttons["applyNumericFilter"].tap()
        XCTAssertTrue(app.staticTexts["2 matches"].exists)
        XCTAssertTrue(app.buttons["compound-0"].exists && app.buttons["compound-1"].exists)
        XCTAssertFalse(app.buttons["compound-2"].exists || app.buttons["compound-3"].exists)
        filter.tap()
        let upper = app.textFields["filterMaximum"]
        upper.tap(); upper.typeText("0")
        app.buttons["applyNumericFilter"].tap()
        XCTAssertTrue(app.staticTexts["1 match"].exists)
        XCTAssertTrue(app.buttons["compound-0"].exists)
        XCTAssertFalse(app.buttons["compound-1"].exists)
        filter.tap()
        let bounds = XCTAttachment(screenshot: app.screenshot())
        bounds.name = "Cached descriptor range editor inline with zero bounds"; bounds.lifetime = .keepAlways; add(bounds)
        app.buttons["filterDescriptor"].tap(); app.buttons["Molecular distances"].tap(); app.buttons["Wiener index"].tap()
        app.buttons["applyNumericFilter"].tap()
        XCTAssertTrue(app.staticTexts["1 match"].exists)
        XCTAssertTrue(app.buttons["compound-0"].exists)
        filter.tap(); app.buttons["clearNumericFilter"].tap()
        XCTAssertTrue(app.staticTexts["4 matches"].exists)
        XCTAssertTrue(app.textFields["compoundSearch"].isHittable)
        let shot = XCTAttachment(screenshot: app.screenshot())
        shot.name = "Numeric filters preserve zero and exclude undefined descriptors"; shot.lifetime = .keepAlways; add(shot)
    }

    @MainActor func testLegacyCountsAndModernSPSReachInlinePortraitAndLandscape() throws {
        continueAfterFailure = false
        XCUIDevice.shared.orientation = .portrait
        defer { XCUIDevice.shared.orientation = .portrait }
        let app = XCUIApplication()
        app.launchArguments = ["--ui-testing", "-AppleLanguages", "(en)", "-AppleLocale", "en_US"]
        app.launch()
        guard app.windows.firstMatch.frame.width >= 700 else { throw XCTSkip("iPad inline stereo properties") }
        app.buttons["pasteSMILES"].tap()
        let editor = app.textViews["smilesEditor"]
        XCTAssertTrue(editor.waitForExistence(timeout: 10))
        editor.tap(); editor.typeText("CC(O)Cl Unspecified\nC[C@H](O)Cl Assigned\nCC1CCC(C)CC1 Dependent ring\n[H][H] Hydrogen\nB1OBOBO1 Boroxine")
        app.buttons["analyzeInput"].tap()
        XCTAssertTrue(app.buttons["compound-0"].waitForExistence(timeout: 30))
        app.buttons["compound-0"].tap()
        app.buttons["descriptorGroup"].tap(); app.buttons["Stereo"].tap()
        XCTAssertEqual(app.staticTexts["descriptor-NumAtomStereoCenters"].label, "1")
        XCTAssertEqual(app.staticTexts["descriptor-NumUnspecifiedAtomStereoCenters"].label, "1")
        XCTAssertEqual(app.staticTexts["descriptor-SPS"].label, "15.7500")
        app.buttons["compound-1"].tap()
        XCTAssertEqual(app.staticTexts["descriptor-NumUnspecifiedAtomStereoCenters"].label, "0")
        XCTAssertEqual(app.staticTexts["descriptor-SPS"].label, "15.7500")
        app.buttons["compound-2"].tap()
        XCTAssertEqual(app.staticTexts["descriptor-NumAtomStereoCenters"].label, "0")
        XCTAssertEqual(app.staticTexts["descriptor-SPS"].label, "39.7500")
        XCTAssertTrue(app.staticTexts["descriptor-SPS"].isHittable)
        XCTAssertTrue(app.textFields["compoundSearch"].isHittable)
        let portrait = XCTAttachment(screenshot: app.screenshot())
        portrait.name = "Legacy counts and modern ring SPS inline"; portrait.lifetime = .keepAlways; add(portrait)
        let search = app.textFields["compoundSearch"]
        search.tap(); search.typeText("Boroxine\n")
        app.buttons["compound-4"].tap()
        XCTAssertEqual(app.staticTexts["descriptor-SPS"].label, "20.0000")
        app.buttons["Clear search"].tap()
        app.buttons["compound-3"].tap()
        XCTAssertEqual(app.staticTexts["descriptor-SPS"].label, "—")
        XCUIDevice.shared.orientation = .landscapeLeft
        let rotated = XCTNSPredicateExpectation(predicate: NSPredicate { _, _ in
            let frame = app.windows.firstMatch.frame; return frame.width > frame.height
        }, object: nil)
        let orientationResult = XCTWaiter.wait(for: [rotated], timeout: 10)
        if orientationResult != .completed {
            let diagnostic = XCTAttachment(screenshot: app.screenshot())
            diagnostic.name = "Rotation diagnostic with actual rendered screen"; diagnostic.lifetime = .keepAlways; add(diagnostic)
            print("Rotation diagnostic device=\(XCUIDevice.shared.orientation.rawValue) window=\(app.windows.firstMatch.frame) screenshot=\(app.screenshot().image.size)")

        }
        XCTAssertEqual(orientationResult, .completed)
        let renderedLandscape = XCTNSPredicateExpectation(predicate: NSPredicate { _, _ in
            let size = XCUIScreen.main.screenshot().image.size
            return size.width > size.height
        }, object: nil)
        XCTAssertEqual(XCTWaiter.wait(for: [renderedLandscape], timeout: 10), .completed)
        XCTAssertTrue(app.staticTexts["descriptor-SPS"].isHittable)
        XCTAssertTrue(app.textFields["compoundSearch"].isHittable)
        XCTAssertFalse(app.buttons["Done"].exists)
        let landscape = XCTAttachment(screenshot: XCUIScreen.main.screenshot())
        landscape.name = "Undefined SPS remains inline in landscape"; landscape.lifetime = .keepAlways; add(landscape)
    }

    @MainActor func testTenThousandOriginalCSVInputsThroughPasteSearchAndInlineProperties() throws {
        struct Expected: Decodable { let lastName: String; let lastFormula: String; let lastWiener: Int; let lastID: Int }
        continueAfterFailure = false
        XCUIDevice.shared.orientation = .portrait
        let bundle = Bundle(for: ChemobUITests.self)
        let csv = try String(contentsOf: XCTUnwrap(bundle.url(forResource: "chebi-10000", withExtension: "csv")), encoding: .utf8)
        let expected = try JSONDecoder().decode(Expected.self, from: Data(contentsOf: XCTUnwrap(bundle.url(forResource: "chebi-10000-expectation", withExtension: "json"))))
        let app = XCUIApplication()
        app.launchArguments = ["--ui-testing", "-AppleLanguages", "(en)", "-AppleLocale", "en_US"]
        app.launch()
        guard app.windows.firstMatch.frame.width >= 700 else { throw XCTSkip("iPad fixed work surface") }
        app.buttons["loadExamples"].tap()
        XCTAssertTrue(app.buttons["compound-0"].waitForExistence(timeout: 30))
        app.buttons["pasteSMILES"].tap()
        let mode = app.segmentedControls["inputMode"]
        XCTAssertTrue(mode.waitForExistence(timeout: 10))
        mode.buttons["Replace library"].tap()
        let format = app.segmentedControls["inputFormat"]
        XCTAssertTrue(format.waitForExistence(timeout: 10), "pasted CSV must use its actual named-column parser")
        format.buttons["CSV"].tap()
        UIPasteboard.general.string = csv
        defer { UIPasteboard.general.string = "" }
        let permission = addUIInterruptionMonitor(withDescription: "Public test-fixture clipboard") { alert in
            let allow = alert.buttons["Allow Paste"]
            if allow.exists { allow.tap(); return true }
            return false
        }
        defer { removeUIInterruptionMonitor(permission) }
        let editor = app.textViews["smilesEditor"]
        editor.tap(); editor.press(forDuration: 1.2)
        let menuPaste = app.menuItems["Paste"]
        if menuPaste.waitForExistence(timeout: 3) { menuPaste.tap() }
        else { app.buttons["Paste"].tap() }
        let analyze = app.buttons["analyzeInput"]
        XCTAssertTrue(analyze.waitForExistence(timeout: 10))
        analyze.tap()
        XCTAssertTrue(app.buttons["Cancel"].waitForExistence(timeout: 10))
        let finished = XCTNSPredicateExpectation(predicate: NSPredicate(format: "exists == false"), object: app.buttons["Cancel"])
        XCTAssertEqual(XCTWaiter.wait(for: [finished], timeout: 240), .completed)
        XCTAssertTrue(app.buttons["exportCSV"].isHittable)
        XCTAssertTrue(app.staticTexts["10,000 matches"].exists)
        app.segmentedControls["compoundCategory"].buttons["Needs review"].tap()
        XCTAssertTrue(app.staticTexts["0 matches"].exists)
        app.segmentedControls["compoundCategory"].buttons["Valid"].tap()
        XCTAssertTrue(app.staticTexts["10,000 matches"].exists)
        let search = app.textFields["compoundSearch"]
        XCTAssertTrue(search.isHittable)
        search.tap(); search.typeText(expected.lastName + "\n")
        let last = app.buttons["compound-\(expected.lastID)"]
        XCTAssertTrue(last.waitForExistence(timeout: 10))
        last.tap()
        XCTAssertTrue(app.staticTexts[expected.lastFormula].exists)
        app.buttons["descriptorGroup"].tap(); app.buttons["Molecular distances"].tap()
        XCTAssertEqual(app.staticTexts["descriptor-ds_Wiener"].label, String(expected.lastWiener))
        XCTAssertTrue(app.staticTexts["descriptor-ds_Wiener"].isHittable)
        XCTAssertTrue(search.isHittable)
        let result = XCTAttachment(screenshot: app.screenshot())
        result.name = "Ten thousand original CSV inputs with last row selected inline"; result.lifetime = .keepAlways; add(result)
        app.buttons["exportCSV"].tap()
        XCTAssertTrue(app.otherElements["csvShareSheet"].waitForExistence(timeout: 30))
        let sharedFile = app.otherElements["LP.CaptionBar.TopCaption"].firstMatch
        XCTAssertTrue(sharedFile.waitForExistence(timeout: 30))
        XCTAssertTrue(sharedFile.label.contains("Chemob-results"))
        XCTAssertTrue(sharedFile.isHittable)
        let share = XCTAttachment(screenshot: XCUIScreen.main.screenshot())
        share.name = "Actual ten thousand row CSV ready for sharing"; share.lifetime = .keepAlways; add(share)
    }

    @MainActor func testAdditionalFiftyDescriptorsReachInlineCategoriesAndUndefinedValues() throws {
        continueAfterFailure = false
        XCUIDevice.shared.orientation = .portrait
        let app = XCUIApplication()
        app.launchArguments = ["--ui-testing", "-AppleLanguages", "(en)", "-AppleLocale", "en_US"]
        app.launch()
        guard app.windows.firstMatch.frame.width >= 700 else { throw XCTSkip("iPad inline properties") }
        app.buttons["pasteSMILES"].tap()
        let editor = app.textViews["smilesEditor"]
        XCTAssertTrue(editor.waitForExistence(timeout: 10))
        editor.tap(); editor.typeText("CCCC Butane\n[H][H] Hydrogen\n[NH3]->[Cu+2]<-[NH3] Copper complex")
        app.buttons["analyzeInput"].tap()
        XCTAssertTrue(app.buttons["compound-0"].waitForExistence(timeout: 30))
        app.buttons["compound-0"].tap()
        app.buttons["descriptorGroup"].tap(); app.buttons["Molecular distances"].tap()
        XCTAssertEqual(app.staticTexts["descriptor-ds_Wiener"].label, "10")
        XCTAssertEqual(app.staticTexts["descriptor-ds_Harary"].label, "4.3333")
        XCTAssertEqual(app.staticTexts["descriptor-ds_EccentricConnectivity"].label, "14")
        XCTAssertEqual(app.staticTexts["descriptor-ds_DegreeDistance"].label, "28")
        XCTAssertEqual(app.staticTexts["descriptor-ds_Gutman"].label, "19")
        XCTAssertEqual(app.staticTexts["descriptor-ds_MDEC11"].label, "0.3333")
        XCTAssertTrue(app.staticTexts["descriptor-ds_MDEC11"].isHittable)
        let distances = XCTAttachment(screenshot: app.screenshot())
        distances.name = "Added molecular distance descriptors inline"; distances.lifetime = .keepAlways; add(distances)
        app.buttons["descriptorGroup"].tap(); app.buttons["Walks"].tap()
        XCTAssertEqual(app.staticTexts["descriptor-ds_MWC3"].label, "2.8332")
        if !app.staticTexts["descriptor-ds_SRW3"].exists { app.buttons["Next properties"].tap() }
        XCTAssertEqual(app.staticTexts["descriptor-ds_SRW3"].label, "0.0000")
        app.buttons["compound-1"].tap()
        XCTAssertEqual(app.staticTexts["descriptor-ds_SRW3"].label, "—")
        app.buttons["compound-0"].tap()
        if !app.staticTexts["descriptor-ds_SRW10"].exists { app.buttons["Next properties"].tap() }
        XCTAssertEqual(app.staticTexts["descriptor-ds_SRW10"].label, "5.5094")
        XCTAssertTrue(app.staticTexts["descriptor-ds_SRW10"].isHittable)
        app.buttons["compound-2"].tap()
        app.buttons["descriptorGroup"].tap(); app.buttons["Molecular distances"].tap()
        XCTAssertEqual(app.staticTexts["descriptor-ds_Wiener"].label, "4")
        XCTAssertTrue(app.staticTexts["Preview unavailable"].waitForExistence(timeout: 10), "a completed unsupported preview must not show an endless loading indicator")
        XCTAssertTrue(app.textFields["compoundSearch"].isHittable)
        XCTAssertFalse(app.buttons["Done"].exists)
        let directed = XCTAttachment(screenshot: app.screenshot())
        directed.name = "Directed coordination input and added descriptor"; directed.lifetime = .keepAlways; add(directed)
    }

    @MainActor func testMacrocycleBridgeRulePreservesSmallAromaticRingsInline() throws {
        continueAfterFailure = false
        XCUIDevice.shared.orientation = .portrait
        let app = XCUIApplication()
        app.launchArguments = ["--ui-testing", "-AppleLanguages", "(en)", "-AppleLocale", "en_US"]
        app.launch()
        guard app.windows.firstMatch.frame.width >= 700 else { throw XCTSkip("iPad inline properties") }
        app.buttons["pasteSMILES"].tap()
        let editor = app.textViews["smilesEditor"]
        XCTAssertTrue(editor.waitForExistence(timeout: 10))
        editor.tap(); editor.typeText("O1C=CC=CC=CC=C1 Oxonine\nO1C=CC=C1 Furan")
        app.buttons["analyzeInput"].tap()
        XCTAssertTrue(app.buttons["compound-0"].waitForExistence(timeout: 30))
        app.buttons["compound-0"].tap()
        app.buttons["descriptorGroup"].tap(); app.buttons["Rings"].tap()
        XCTAssertEqual(app.staticTexts["descriptor-NumAromaticRings"].label, "0")
        XCTAssertEqual(app.staticTexts["descriptor-NumAliphaticRings"].label, "1")
        let macro = XCTAttachment(screenshot: app.screenshot())
        macro.name = "Nine-member ether ring remains nonaromatic"; macro.lifetime = .keepAlways; add(macro)
        app.buttons["compound-1"].tap()
        XCTAssertEqual(app.staticTexts["descriptor-NumAromaticRings"].label, "1")
        XCTAssertEqual(app.staticTexts["descriptor-NumAliphaticRings"].label, "0")
        XCTAssertTrue(app.staticTexts["O1C=CC=C1"].exists)
        XCTAssertTrue(app.textFields["compoundSearch"].isHittable)
        XCTAssertFalse(app.buttons["Done"].exists)
        let small = XCTAttachment(screenshot: app.screenshot())
        small.name = "Small furan ring remains aromatic"; small.lifetime = .keepAlways; add(small)
    }
    @MainActor func testQEDUsesNativeHydrogenAndElementRulesInline() throws {
        continueAfterFailure = false
        XCUIDevice.shared.orientation = .portrait
        let app = XCUIApplication()
        app.launchArguments = ["--ui-testing", "-AppleLanguages", "(en)", "-AppleLocale", "en_US"]
        app.launch()
        guard app.windows.firstMatch.frame.width >= 700 else { throw XCTSkip("iPad inline properties") }
        app.buttons["pasteSMILES"].tap()
        let editor = app.textViews["smilesEditor"]
        XCTAssertTrue(editor.waitForExistence(timeout: 10))
        editor.tap(); editor.typeText("[H]OOC Hydroperoxide\nc1cc[asH]c1 Arsenic ring\n[Hg] Mercury\nc1cc[te]c1 Tellurophene")
        app.buttons["analyzeInput"].tap()
        XCTAssertTrue(app.buttons["compound-0"].waitForExistence(timeout: 30))
        app.buttons["compound-0"].tap()
        app.buttons["descriptorGroup"].tap(); app.buttons["Screening"].tap()
        XCTAssertEqual(app.staticTexts["descriptor-qed"].label, "0.3106")
        XCTAssertTrue(app.staticTexts["descriptor-qed"].isHittable)
        XCTAssertTrue(app.staticTexts["[H]OOC"].exists)
        let hydrogen = XCTAttachment(screenshot: app.screenshot())
        hydrogen.name = "QED default hydrogen removal with original input"; hydrogen.lifetime = .keepAlways; add(hydrogen)
        app.buttons["compound-2"].tap()
        XCTAssertEqual(app.staticTexts["descriptor-qed"].label, "0.4831")
        app.buttons["compound-3"].tap()
        XCTAssertEqual(app.staticTexts["descriptor-qed"].label, "0.5159")
        app.buttons["compound-1"].tap()
        XCTAssertEqual(app.staticTexts["descriptor-qed"].label, "0.4156")
        app.buttons["descriptorGroup"].tap(); app.buttons["Rings"].tap()
        XCTAssertEqual(app.staticTexts["descriptor-NumAromaticRings"].label, "0")
        XCTAssertTrue(app.staticTexts["c1cc[asH]c1"].exists)
        XCTAssertTrue(app.textFields["compoundSearch"].isHittable)
        XCTAssertFalse(app.buttons["Done"].exists)
        let arsenic = XCTAttachment(screenshot: app.screenshot())
        arsenic.name = "Arsenic ring uses prepared nonaromatic values"; arsenic.lifetime = .keepAlways; add(arsenic)
    }
    @MainActor func testMorganDensitiesAndRadicalPreparationUseNativeInlineValues() throws {
        continueAfterFailure = false
        XCUIDevice.shared.orientation = .portrait
        let app = XCUIApplication()
        app.launchArguments = ["--ui-testing", "-AppleLanguages", "(en)", "-AppleLocale", "en_US"]
        app.launch()
        guard app.windows.firstMatch.frame.width >= 700 else { throw XCTSkip("iPad inline properties") }
        app.buttons["pasteSMILES"].tap()
        let editor = app.textViews["smilesEditor"]
        XCTAssertTrue(editor.waitForExistence(timeout: 10))
        editor.tap(); editor.typeText("CCO Ethanol\n[n+]1ccccc1 Pyridinyl cation")
        app.buttons["analyzeInput"].tap()
        XCTAssertTrue(app.buttons["compound-0"].waitForExistence(timeout: 30))
        app.buttons["compound-0"].tap()
        app.buttons["descriptorGroup"].tap(); app.buttons["Fingerprints"].tap()
        for radius in 1...3 {
            let value = app.staticTexts["descriptor-FpDensityMorgan\(radius)"]
            XCTAssertEqual(value.label, "2.0000"); XCTAssertTrue(value.isHittable)
        }
        app.buttons["compound-1"].tap()
        for (radius, expected) in [(1,"1.0000"),(2,"2.0000"),(3,"2.1667")] {
            XCTAssertEqual(app.staticTexts["descriptor-FpDensityMorgan\(radius)"].label, expected)
        }
        let keys = XCTAttachment(screenshot: app.screenshot())
        keys.name = "Three native unfolded Morgan densities"; keys.lifetime = .keepAlways; add(keys)
        app.buttons["descriptorGroup"].tap(); app.buttons["Essentials"].tap()
        XCTAssertEqual(app.staticTexts["descriptor-NumHAcceptors"].label, "1")
        XCTAssertEqual(app.staticTexts["descriptor-MolLogP"].label, "0.48")
        XCTAssertTrue(app.staticTexts["[n+]1ccccc1"].exists)
        XCTAssertTrue(app.textFields["compoundSearch"].isHittable)
        XCTAssertFalse(app.buttons["Done"].exists)
        let preparation = XCTAttachment(screenshot: app.screenshot())
        preparation.name = "Corrected radical properties with original input"; preparation.lifetime = .keepAlways; add(preparation)
    }
    @MainActor func testCoordinationCountsAndUndefinedBCUTReachInlineResults() throws {
        continueAfterFailure = false
        XCUIDevice.shared.orientation = .portrait
        let app = XCUIApplication()
        app.launchArguments = ["--ui-testing", "-AppleLanguages", "(en)", "-AppleLocale", "en_US"]
        app.launch()
        guard app.windows.firstMatch.frame.width >= 700 else { throw XCTSkip("iPad inline properties") }
        let cerium = "O=N1=[O][Ce-2]23456([O]1)([O]N(=O)=[O]2)([O]N(=O)=[O]3)([O]N(=O)=[O]4)([O]N(=O)=[O]5)[O]N(=O)=[O]6"
        app.buttons["pasteSMILES"].tap()
        let editor = app.textViews["smilesEditor"]
        XCTAssertTrue(editor.waitForExistence(timeout: 10))
        editor.tap(); editor.typeText(cerium + " Cerium nitrate\n[OH2][Mg+2] Coordinated water")
        app.buttons["analyzeInput"].tap()
        XCTAssertTrue(app.buttons["compound-0"].waitForExistence(timeout: 30))
        app.buttons["compound-0"].tap()
        XCTAssertEqual(app.staticTexts["descriptor-NumHAcceptors"].label, "18")
        XCTAssertEqual(app.staticTexts["descriptor-NumRotatableBonds"].label, "12")
        app.buttons["descriptorGroup"].tap(); app.buttons["Rings"].tap()
        XCTAssertEqual(app.staticTexts["descriptor-RingCount"].label, "0")
        XCTAssertTrue(app.staticTexts["descriptor-RingCount"].isHittable)
        let ring = XCTAttachment(screenshot: app.screenshot())
        ring.name = "Coordination edges excluded from default ring count"; ring.lifetime = .keepAlways; add(ring)
        app.buttons["compound-1"].tap()
        app.buttons["descriptorGroup"].tap(); app.buttons["BCUT"].tap()
        for id in ["MWHI", "MWLOW", "CHGHI", "CHGLO", "LOGPHI", "LOGPLOW", "MRHI", "MRLOW"] {
            let value = app.staticTexts["descriptor-BCUT2D_" + id]
            XCTAssertEqual(value.label, "—"); XCTAssertTrue(value.isHittable)
        }
        XCTAssertTrue(app.staticTexts["[OH2][Mg+2]"].exists)
        XCTAssertTrue(app.textFields["compoundSearch"].isHittable)
        XCTAssertFalse(app.buttons["Done"].exists)
        let undefined = XCTAttachment(screenshot: app.screenshot())
        undefined.name = "Undefined BCUT retained for dative bonds"; undefined.lifetime = .keepAlways; add(undefined)
    }
    @MainActor func testCleanedChemistryAndExplicitHydrogenRotorReachInlineResults() throws {
        continueAfterFailure = false
        XCUIDevice.shared.orientation = .portrait
        let app = XCUIApplication()
        app.launchArguments = ["--ui-testing", "-AppleLanguages", "(en)", "-AppleLocale", "en_US"]
        app.launch()
        guard app.windows.firstMatch.frame.width >= 700 else { throw XCTSkip("iPad inline properties") }
        app.buttons["pasteSMILES"].tap()
        let editor = app.textViews["smilesEditor"]
        XCTAssertTrue(editor.waitForExistence(timeout: 10))
        editor.tap(); editor.typeText("CN(=O)=O Nitromethane\n[2H]ON(=O)=O Deuterated nitric acid")
        app.buttons["analyzeInput"].tap()
        XCTAssertTrue(app.buttons["compound-0"].waitForExistence(timeout: 30))
        app.buttons["compound-0"].tap()
        XCTAssertEqual(app.staticTexts["descriptor-MolLogP"].label, "-0.11")
        XCTAssertEqual(app.staticTexts["descriptor-TPSA"].label, "43.14")
        XCTAssertTrue(app.staticTexts["CN(=O)=O"].exists)
        let nitro = XCTAttachment(screenshot: app.screenshot())
        nitro.name = "Corrected nitromethane properties with original source"; nitro.lifetime = .keepAlways; add(nitro)
        app.buttons["compound-1"].tap()
        XCTAssertEqual(app.staticTexts["descriptor-NumRotatableBonds"].label, "1")
        XCTAssertEqual(app.staticTexts["descriptor-TPSA"].label, "63.37")
        XCTAssertTrue(app.staticTexts["[2H]ON(=O)=O"].exists)
        XCTAssertTrue(app.textFields["compoundSearch"].isHittable)
        XCTAssertFalse(app.buttons["Done"].exists)
        let hydrogen = XCTAttachment(screenshot: app.screenshot())
        hydrogen.name = "Explicit isotope H strict rotor with original source"; hydrogen.lifetime = .keepAlways; add(hydrogen)
    }
    @MainActor func testStructuralCountsUseNativeValuesAcrossInlinePages() throws {
        continueAfterFailure = false
        XCUIDevice.shared.orientation = .portrait
        let app = XCUIApplication()
        app.launchArguments = ["--ui-testing", "-AppleLanguages", "(en)", "-AppleLocale", "en_US"]
        app.launch()
        guard app.windows.firstMatch.frame.width >= 700 else { throw XCTSkip("iPad inline properties") }
        func reveal(_ id: String) {
            for _ in 0..<40 {
                if app.staticTexts["descriptor-" + id].isHittable { return }
                let next = app.buttons["Next properties"]
                guard next.exists && next.isEnabled else { return }
                next.tap()
            }
        }
        app.buttons["pasteSMILES"].tap()
        let editor = app.textViews["smilesEditor"]
        XCTAssertTrue(editor.waitForExistence(timeout: 10))
        editor.tap(); editor.typeText("NC(=O)N Urea\nC1CC2CCC1C2 Norbornane\nC1CCC2(CC1)CCC2 Spiro compound")
        app.buttons["analyzeInput"].tap()
        XCTAssertTrue(app.buttons["compound-0"].waitForExistence(timeout: 30))
        app.buttons["compound-0"].tap()
        app.buttons["descriptorGroup"].tap()
        app.buttons["Functional groups"].tap()
        reveal("NumAmideBonds")
        XCTAssertEqual(app.staticTexts["descriptor-NumAmideBonds"].label, "2")
        XCTAssertTrue(app.staticTexts["descriptor-NumAmideBonds"].isHittable)
        let amide = XCTAttachment(screenshot: app.screenshot())
        amide.name = "Native two amide bonds"; amide.lifetime = .keepAlways; add(amide)
        app.buttons["compound-1"].tap()
        app.buttons["descriptorGroup"].tap()
        app.buttons["Rings"].tap()
        reveal("NumSpiroAtoms")
        XCTAssertEqual(app.staticTexts["descriptor-NumSpiroAtoms"].label, "0")
        reveal("NumBridgeheadAtoms")
        XCTAssertEqual(app.staticTexts["descriptor-NumBridgeheadAtoms"].label, "2")
        XCTAssertTrue(app.staticTexts["descriptor-NumBridgeheadAtoms"].isHittable)
        let bridge = XCTAttachment(screenshot: app.screenshot())
        bridge.name = "Native two bridgehead atoms"; bridge.lifetime = .keepAlways; add(bridge)
        app.buttons["compound-2"].tap()
        app.buttons["descriptorGroup"].tap(); app.buttons["Essentials"].tap()
        app.buttons["descriptorGroup"].tap(); app.buttons["Rings"].tap()
        reveal("NumSpiroAtoms")
        XCTAssertEqual(app.staticTexts["descriptor-NumSpiroAtoms"].label, "1")
        XCTAssertTrue(app.staticTexts["descriptor-NumSpiroAtoms"].isHittable)
        XCTAssertTrue(app.textFields["compoundSearch"].isHittable)
        XCTAssertFalse(app.buttons["Done"].exists)
        let spiro = XCTAttachment(screenshot: app.screenshot())
        spiro.name = "Native one spiro atom"; spiro.lifetime = .keepAlways; add(spiro)
    }
    @MainActor func testBCUTDescriptorsUseNativeValuesInline() throws {
        continueAfterFailure = false
        XCUIDevice.shared.orientation = .portrait
        let app = XCUIApplication()
        app.launchArguments = ["--ui-testing", "-AppleLanguages", "(en)", "-AppleLocale", "en_US"]
        app.launch()
        guard app.windows.firstMatch.frame.width >= 700 else { throw XCTSkip("iPad inline properties") }
        app.buttons["loadExamples"].tap()
        XCTAssertTrue(app.buttons["compound-0"].waitForExistence(timeout: 30))
        app.buttons["compound-0"].tap()
        app.buttons["descriptorGroup"].tap()
        app.buttons["BCUT"].tap()
        let references = [
            "BCUT2D_MWHI": "16.2490", "BCUT2D_MWLOW": "10.9082",
            "BCUT2D_CHGHI": "1.3403", "BCUT2D_CHGLO": "-1.5224",
            "BCUT2D_LOGPHI": "1.2962", "BCUT2D_LOGPLOW": "-1.5684",
            "BCUT2D_MRHI": "3.8358", "BCUT2D_MRLOW": "0.3184"
        ]
        for (id, expected) in references {
            XCTAssertEqual(app.staticTexts["descriptor-" + id].label, expected)
            XCTAssertTrue(app.staticTexts["descriptor-" + id].isHittable)
        }
        XCTAssertTrue(app.textFields["compoundSearch"].isHittable)
        XCTAssertFalse(app.buttons["Done"].exists)
        let shot = XCTAttachment(screenshot: app.screenshot())
        shot.name = "Native BCUT eight descriptors"; shot.lifetime = .keepAlways; add(shot)
    }
    @MainActor func testChargeDescriptorsUseNativeValuesAcrossInlinePages() throws {
        continueAfterFailure = false
        XCUIDevice.shared.orientation = .portrait
        let app = XCUIApplication()
        app.launchArguments = ["--ui-testing", "-AppleLanguages", "(en)", "-AppleLocale", "en_US"]
        app.launch()
        guard app.windows.firstMatch.frame.width >= 700 else { throw XCTSkip("iPad inline properties") }
        app.buttons["loadExamples"].tap()
        XCTAssertTrue(app.buttons["compound-0"].waitForExistence(timeout: 30))
        app.buttons["compound-0"].tap()
        app.buttons["descriptorGroup"].tap()
        app.buttons["Partial charge"].tap()
        XCTAssertEqual(app.staticTexts["descriptor-MaxPartialCharge"].label, "0.0402")
        XCTAssertEqual(app.staticTexts["descriptor-MinPartialCharge"].label, "-0.3967")
        XCTAssertEqual(app.staticTexts["descriptor-MaxAbsPartialCharge"].label, "0.3967")
        XCTAssertEqual(app.staticTexts["descriptor-MinAbsPartialCharge"].label, "0.0402")
        let charge = XCTAttachment(screenshot: app.screenshot())
        charge.name = "Native partial charge extrema"; charge.lifetime = .keepAlways; add(charge)
        app.buttons["descriptorGroup"].tap()
        app.buttons["Charge surface"].tap()
        XCTAssertEqual(app.staticTexts["descriptor-PEOE_VSA1"].label, "5.1065")
        XCTAssertEqual(app.staticTexts["descriptor-PEOE_VSA7"].label, "6.9237")
        XCTAssertEqual(app.staticTexts["descriptor-PEOE_VSA8"].label, "6.6069")
        let firstSurface = XCTAttachment(screenshot: app.screenshot())
        firstSurface.name = "Native charge surface values"; firstSurface.lifetime = .keepAlways; add(firstSurface)
        app.buttons["Next properties"].tap()
        XCTAssertEqual(app.staticTexts["descriptor-PEOE_VSA13"].label, "0.0000")
        XCTAssertEqual(app.staticTexts["descriptor-PEOE_VSA14"].label, "0.0000")
        XCTAssertTrue(app.textFields["compoundSearch"].isHittable)
        XCTAssertFalse(app.buttons["Done"].exists)
        let surface = XCTAttachment(screenshot: app.screenshot())
        surface.name = "Native charge surface last page"; surface.lifetime = .keepAlways; add(surface)
    }
    @MainActor func testRadicalElectronCountUsesNativeContext() throws {
        continueAfterFailure = false
        XCUIDevice.shared.orientation = .portrait
        let app = XCUIApplication()
        app.launchArguments = ["--ui-testing", "-AppleLanguages", "(en)", "-AppleLocale", "en_US"]
        app.launch()
        guard app.windows.firstMatch.frame.width >= 700 else { throw XCTSkip("iPad inline properties") }
        app.buttons["pasteSMILES"].tap()
        let editor = app.textViews["smilesEditor"]
        XCTAssertTrue(editor.waitForExistence(timeout: 10))
        editor.tap(); editor.typeText("[CH3] Methyl radical")
        app.buttons["analyzeInput"].tap()
        XCTAssertTrue(app.buttons["compound-0"].waitForExistence(timeout: 30))
        app.buttons["compound-0"].tap()
        XCTAssertEqual(app.staticTexts["descriptor-FractionCSP3"].label, "0.00")
        app.buttons["descriptorGroup"].tap()
        app.buttons["Composition"].tap()
        app.buttons["Next properties"].tap()
        XCTAssertEqual(app.staticTexts["descriptor-NumRadicalElectrons"].label, "1")
        XCTAssertTrue(app.staticTexts["descriptor-NumRadicalElectrons"].isHittable)
        XCTAssertTrue(app.textFields["compoundSearch"].isHittable)
        let shot = XCTAttachment(screenshot: app.screenshot())
        shot.name = "Native radical electron count"; shot.lifetime = .keepAlways; add(shot)
    }
    @MainActor func testTopologyDescriptorsUseNativeValuesAcrossInlinePages() throws {
        continueAfterFailure = false
        XCUIDevice.shared.orientation = .portrait
        let app = XCUIApplication()
        app.launchArguments = ["--ui-testing", "-AppleLanguages", "(en)", "-AppleLocale", "en_US"]
        app.launch()
        guard app.windows.firstMatch.frame.width >= 700 else { throw XCTSkip("iPad inline properties") }
        app.buttons["loadExamples"].tap()
        XCTAssertTrue(app.buttons["compound-0"].waitForExistence(timeout: 30))
        app.buttons["compound-0"].tap()
        app.buttons["descriptorGroup"].tap()
        app.buttons["Topology"].tap()
        XCTAssertEqual(app.staticTexts["descriptor-BalabanJ"].label, "1.6330")
        XCTAssertEqual(app.staticTexts["descriptor-BertzCT"].label, "2.7549")
        XCTAssertEqual(app.staticTexts["descriptor-Chi1v"].label, "1.0233")
        XCTAssertEqual(app.staticTexts["descriptor-Chi3n"].label, "0.0000")
        XCTAssertEqual(app.staticTexts["descriptor-Chi3v"].label, "0.0000")
        let first = XCTAttachment(screenshot: app.screenshot())
        first.name = "Native topology connectivity"; first.lifetime = .keepAlways; add(first)
        app.buttons["Next properties"].tap()
        XCTAssertEqual(app.staticTexts["descriptor-HallKierAlpha"].label, "-0.0400")
        XCTAssertEqual(app.staticTexts["descriptor-Kappa1"].label, "2.9600")
        XCTAssertEqual(app.staticTexts["descriptor-Kappa2"].label, "1.9600")
        XCTAssertEqual(app.staticTexts["descriptor-Kappa3"].label, "1.9600")
        XCTAssertEqual(app.staticTexts["descriptor-Phi"].label, "1.9339")
        XCTAssertTrue(app.staticTexts["descriptor-Phi"].isHittable)
        XCTAssertTrue(app.textFields["compoundSearch"].isHittable)
        XCTAssertFalse(app.buttons["Done"].exists)
        let last = XCTAttachment(screenshot: app.screenshot())
        last.name = "Native topology shape and flexibility"; last.lifetime = .keepAlways; add(last)
    }
    @MainActor func testPartialExplicitHydrogenKeepsOriginalSurfaceBins() throws {
        continueAfterFailure = false
        XCUIDevice.shared.orientation = .portrait
        let app = XCUIApplication()
        app.launchArguments = ["--ui-testing", "-AppleLanguages", "(en)", "-AppleLocale", "en_US"]
        app.launch()
        guard app.windows.firstMatch.frame.width >= 700 else { throw XCTSkip("iPad inline properties") }
        app.buttons["pasteSMILES"].tap()
        let editor = app.textViews["smilesEditor"]
        XCTAssertTrue(editor.waitForExistence(timeout: 10))
        editor.tap(); editor.typeText("[2H]O Isotopic water")
        app.buttons["analyzeInput"].tap()
        XCTAssertTrue(app.buttons["compound-0"].waitForExistence(timeout: 30))
        app.buttons["compound-0"].tap()
        XCTAssertEqual(app.staticTexts["descriptor-MolLogP"].label, "-0.82")
        app.buttons["descriptorGroup"].tap()
        app.buttons["Surface area"].tap()
        XCTAssertEqual(app.staticTexts["descriptor-SMR_VSA1"].label, "6.9097")
        XCTAssertEqual(app.staticTexts["descriptor-SMR_VSA2"].label, "0.0000")
        XCTAssertTrue(app.textFields["compoundSearch"].isHittable)
        XCTAssertFalse(app.buttons["Done"].exists)
        let shot = XCTAttachment(screenshot: app.screenshot())
        shot.name = "Original explicit-H surface bins"; shot.lifetime = .keepAlways; add(shot)
        app.buttons["descriptorGroup"].tap()
        app.buttons["BCUT"].tap()
        XCTAssertEqual(app.staticTexts["descriptor-BCUT2D_MWHI"].label, "15.9990")
        XCTAssertEqual(app.staticTexts["descriptor-BCUT2D_MWLOW"].label, "15.9990")
        XCTAssertEqual(app.staticTexts["descriptor-BCUT2D_CHGHI"].label, "-0.4115")
        XCTAssertTrue(app.textFields["compoundSearch"].isHittable)
        let bcut = XCTAttachment(screenshot: app.screenshot())
        bcut.name = "BCUT removes isotope hydrogen without changing original surface"; bcut.lifetime = .keepAlways; add(bcut)
    }
    @MainActor func testSurfaceDescriptorsUseNativeValuesAndInlinePages() throws {
        continueAfterFailure = false
        XCUIDevice.shared.orientation = .portrait
        let app = XCUIApplication()
        app.launchArguments = ["--ui-testing", "-AppleLanguages", "(en)", "-AppleLocale", "en_US"]
        app.launch()
        guard app.windows.firstMatch.frame.width >= 700 else { throw XCTSkip("iPad inline properties") }
        app.buttons["loadExamples"].tap()
        XCTAssertTrue(app.buttons["compound-0"].waitForExistence(timeout: 30))
        app.buttons["compound-0"].tap()
        app.buttons["descriptorGroup"].tap()
        app.buttons["Surface area"].tap()
        XCTAssertEqual(app.staticTexts["descriptor-LabuteASA"].label, "19.8984")
        XCTAssertEqual(app.staticTexts["descriptor-SMR_VSA5"].label, "6.9237")
        XCTAssertTrue(app.staticTexts["descriptor-LabuteASA"].isHittable)
        let area = XCTAttachment(screenshot: app.screenshot())
        area.name = "Native surface area bins"; area.lifetime = .keepAlways; add(area)
        app.buttons["Next properties"].tap()
        XCTAssertEqual(app.staticTexts["descriptor-SlogP_VSA2"].label, "11.7134")
        XCTAssertEqual(app.staticTexts["descriptor-SlogP_VSA5"].label, "6.9237")
        XCTAssertFalse(app.staticTexts["descriptor-LabuteASA"].exists)
        app.buttons["descriptorGroup"].tap()
        app.buttons["Electronic surface"].tap()
        XCTAssertEqual(app.staticTexts["descriptor-EState_VSA2"].label, "6.6069")
        XCTAssertEqual(app.staticTexts["descriptor-EState_VSA9"].label, "5.1065")
        app.buttons["Next properties"].tap()
        XCTAssertEqual(app.staticTexts["descriptor-VSA_EState3"].label, "7.5694")
        XCTAssertEqual(app.staticTexts["descriptor-VSA_EState8"].label, "1.9306")
        XCTAssertTrue(app.staticTexts["descriptor-VSA_EState8"].isHittable)
        XCTAssertTrue(app.textFields["compoundSearch"].isHittable)
        XCTAssertFalse(app.buttons["Done"].exists)
        let electronic = XCTAttachment(screenshot: app.screenshot())
        electronic.name = "Native surface EState bins"; electronic.lifetime = .keepAlways; add(electronic)
    }
    @MainActor func testFusedAromaticityCorrectionReachesInlinePanel() throws {
        continueAfterFailure = false
        XCUIDevice.shared.orientation = .portrait
        let app = XCUIApplication()
        app.launchArguments = ["--ui-testing", "-AppleLanguages", "(en)", "-AppleLocale", "en_US"]
        app.launch()
        guard app.windows.firstMatch.frame.width >= 700 else { throw XCTSkip("iPad inline properties") }
        app.buttons["pasteSMILES"].tap()
        let editor = app.textViews["smilesEditor"]
        XCTAssertTrue(editor.waitForExistence(timeout: 10))
        editor.tap()
        editor.typeText("ClC1=CC2=NC3=C(CCCC3)C(=C2C=C1)Cl Acridine regression")
        app.buttons["analyzeInput"].tap()
        XCTAssertTrue(app.buttons["compound-0"].waitForExistence(timeout: 30))
        app.buttons["compound-0"].tap()
        app.buttons["descriptorGroup"].tap()
        app.buttons["Graph"].tap()
        // Independent RDKit NCI-1218 references. The old engine showed 6/6.
        XCTAssertEqual(app.staticTexts["descriptor-ds_AromaticBonds"].label, "11")
        XCTAssertEqual(app.staticTexts["descriptor-ds_AromaticAtoms"].label, "10")
        XCTAssertTrue(app.textFields["compoundSearch"].isHittable)
        XCTAssertFalse(app.buttons["Done"].exists)
        let shot = XCTAttachment(screenshot: app.screenshot())
        shot.name = "Corrected fused aromaticity"; shot.lifetime = .keepAlways; add(shot)
        app.buttons["descriptorGroup"].tap()
        app.buttons["Essentials"].tap()
        XCTAssertEqual(app.staticTexts["descriptor-MolLogP"].label, "4.42")
    }
    @MainActor func testFunctionalGroupsAndElectronicIndices() throws {
        continueAfterFailure = false
        XCUIDevice.shared.orientation = .portrait
        let app = XCUIApplication()
        app.launchArguments = ["--ui-testing", "-AppleLanguages", "(en)", "-AppleLocale", "en_US"]
        app.launch()
        guard app.windows.firstMatch.frame.width >= 700 else { throw XCTSkip("iPad inline properties") }
        app.buttons["loadExamples"].tap()
        XCTAssertTrue(app.buttons["compound-0"].waitForExistence(timeout: 30))
        app.buttons["compound-0"].tap()
        app.buttons["descriptorGroup"].tap()
        app.buttons["Functional groups"].tap()
        let hydroxyl = app.staticTexts["descriptor-fr_Al_OH"]
        XCTAssertTrue(hydroxyl.waitForExistence(timeout: 5))
        XCTAssertEqual(hydroxyl.label, "1")
        XCTAssertTrue(hydroxyl.isHittable)
        let functional = XCTAttachment(screenshot: app.screenshot())
        functional.name = "Native functional groups"; functional.lifetime = .keepAlways; add(functional)
        app.buttons["Next properties"].tap()
        XCTAssertTrue(app.staticTexts["descriptor-fr_C_S"].isHittable)
        XCTAssertFalse(app.staticTexts["descriptor-fr_Al_OH"].exists)
        app.buttons["descriptorGroup"].tap()
        app.buttons["Electronic"].tap()
        let maximum = app.staticTexts["descriptor-MaxEStateIndex"]
        XCTAssertTrue(maximum.waitForExistence(timeout: 5))
        XCTAssertEqual(maximum.label, "7.5694")
        XCTAssertEqual(app.staticTexts["descriptor-MinEStateIndex"].label, "0.2500")
        XCTAssertTrue(app.staticTexts["descriptor-MinAbsEStateIndex"].isHittable)
        XCTAssertTrue(app.textFields["compoundSearch"].isHittable)
        XCTAssertFalse(app.buttons["Done"].exists)
        let electronic = XCTAttachment(screenshot: app.screenshot())
        electronic.name = "Native electronic indices"; electronic.lifetime = .keepAlways; add(electronic)
    }
    @MainActor func testIPadSelectionStaysBelowSearch() throws {
        continueAfterFailure = false
        defer { XCUIDevice.shared.orientation = .portrait }
        XCUIDevice.shared.orientation = .portrait
        let app = XCUIApplication()
        app.launchArguments = ["--ui-testing", "-AppleLanguages", "(en)", "-AppleLocale", "en_US"]
        app.launch()
        let portrait = XCTNSPredicateExpectation(predicate: NSPredicate { _, _ in
            let frame = app.windows.firstMatch.frame
            return frame.height > frame.width
        }, object: nil)
        XCTAssertEqual(XCTWaiter.wait(for: [portrait], timeout: 10), .completed)
        guard app.windows.firstMatch.frame.width >= 700 else { throw XCTSkip("iPad layout") }
        app.buttons["loadExamples"].tap()
        XCTAssertTrue(app.buttons["compound-0"].waitForExistence(timeout: 30))
        app.buttons["compound-0"].tap()
        let detail = app.otherElements["inlineCompoundDetail"]
        XCTAssertTrue(detail.waitForExistence(timeout: 3))
        XCTAssertFalse(app.buttons["Done"].exists)
        XCTAssertTrue(app.textFields["compoundSearch"].isHittable)
        XCTAssertTrue(app.staticTexts["Molecular weight"].waitForExistence(timeout: 10))
        XCTAssertTrue(app.staticTexts["Molecular weight"].isHittable)
        XCTAssertTrue(app.staticTexts["Fraction sp³ carbons"].isHittable)
        app.buttons["descriptorGroup"].tap()
        app.buttons["Graph"].tap()
        let refractivity = app.staticTexts["descriptor-MolMR"]
        XCTAssertTrue(refractivity.waitForExistence(timeout: 5))
        XCTAssertEqual(refractivity.label, "12.76", "Ethanol MR must come from the extended native vector")
        XCTAssertTrue(refractivity.isHittable)
        let graphProperties = XCTAttachment(screenshot: app.screenshot())
        graphProperties.name = "Extended graph properties"
        graphProperties.lifetime = .keepAlways
        add(graphProperties)
        app.buttons["descriptorGroup"].tap()
        app.buttons["Composition"].tap()
        XCTAssertTrue(app.staticTexts["descriptor-NumValenceElectrons"].isHittable)
        app.buttons["Next properties"].tap()
        XCTAssertTrue(app.staticTexts["descriptor-ds_CountCl"].isHittable)
        app.buttons["descriptorGroup"].tap()
        app.buttons["Essentials"].tap()
        XCTAssertGreaterThan(detail.frame.minY, app.windows.firstMatch.frame.midY - 70)
        XCTAssertTrue(app.staticTexts["On device"].isHittable)
        app.buttons["Needs review"].tap()
        XCTAssertTrue(app.staticTexts["Invalid or unsupported structure. The original SMILES is preserved."].waitForExistence(timeout: 5))
        app.buttons["All"].tap()
        let search = app.textFields["compoundSearch"]
        search.tap(); search.typeText("Caffeine\n")
        XCTAssertTrue(app.staticTexts["194.19"].firstMatch.waitForExistence(timeout: 5))
        XCTAssertFalse(app.buttons["compound-0"].exists)
        XCTAssertEqual(app.staticTexts["selectedCompoundName"].label, "Caffeine")
        let shot = XCTAttachment(screenshot: app.screenshot()); shot.name = "Inline iPad workbench"; shot.lifetime = .keepAlways; add(shot)
        XCUIDevice.shared.orientation = .landscapeLeft
        let rotated = XCTNSPredicateExpectation(predicate: NSPredicate { _, _ in
            let frame = app.windows.firstMatch.frame
            return frame.width > frame.height
        }, object: nil)
        XCTAssertEqual(XCTWaiter.wait(for: [rotated], timeout: 10), .completed)
        // The window bounds update before the Simulator rotation animation completes.
        Thread.sleep(forTimeInterval: 1)
        app.buttons["Next properties"].tap()
        XCTAssertTrue(app.staticTexts["Fraction sp³ carbons"].isHittable)
        XCTAssertTrue(app.textFields["compoundSearch"].isHittable)
        XCTAssertLessThanOrEqual(detail.frame.maxY, app.windows.firstMatch.frame.maxY)
        let landscape = XCTAttachment(screenshot: XCUIScreen.main.screenshot()); landscape.name = "Inline iPad landscape"; landscape.lifetime = .keepAlways; add(landscape)
    }
    @MainActor func testSampleSelectionShowsNativeDescriptors() throws {
        let app = XCUIApplication()
        app.launchArguments = ["--ui-testing", "-AppleLanguages", "(en)", "-AppleLocale", "en_US"]
        app.launch()
        guard app.windows.firstMatch.frame.width < 700 else { throw XCTSkip("iPhone flow; iPad has its inline selection test") }
        let home = XCTAttachment(screenshot: app.screenshot()); home.name = "Chemob Home"; home.lifetime = .keepAlways; add(home)
        app.buttons["loadExamples"].tap()
        XCTAssertTrue(app.buttons["compound-0"].waitForExistence(timeout: 30))
        XCTAssertTrue(app.staticTexts["46.07"].firstMatch.waitForExistence(timeout: 30))
        XCTAssertTrue(app.buttons["compound-5"].isHittable, "Six sample inputs should fit on the phone")
        let search = app.textFields["compoundSearch"]
        search.tap(); search.typeText("Caffeine")
        XCTAssertTrue(app.buttons["compound-2"].exists)
        XCTAssertFalse(app.buttons["compound-0"].exists)
        app.buttons["Clear search"].tap()
        search.typeText("\n")
        app.buttons["compound-0"].tap()
        XCTAssertTrue(app.staticTexts["C2H6O"].waitForExistence(timeout: 10))
        XCTAssertTrue(app.staticTexts["46.07"].exists)
        let detail = XCTAttachment(screenshot: app.screenshot()); detail.name = "Chemob Native Analysis"; detail.lifetime = .keepAlways; add(detail)
        app.buttons["Done"].tap()
        let invalid = app.buttons["compound-5"]
        if !invalid.isHittable { app.swipeUp() }
        XCTAssertTrue(invalid.waitForExistence(timeout: 10))
        invalid.tap()
        XCTAssertTrue(app.staticTexts["Invalid or unsupported structure. The original SMILES is preserved."].waitForExistence(timeout: 10))
        app.buttons["Done"].tap()
        app.swipeDown()
        app.buttons["pasteSMILES"].tap()
        let editor = app.textViews["smilesEditor"]
        XCTAssertTrue(editor.waitForExistence(timeout: 10))
        editor.tap(); editor.typeText("CCO Reused ethanol")
        app.buttons["analyzeInput"].tap()
        XCTAssertTrue(app.staticTexts["Reused 7 results"].waitForExistence(timeout: 30))
        app.swipeUp()
        XCTAssertTrue(app.buttons["compound-6"].waitForExistence(timeout: 10))
    }
}
