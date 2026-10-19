import Foundation
import ChemobCore

@main struct CacheIntegration {
    static func main() async throws {
        precondition(Chemistry.revision == "ducksmiles-workbench-v42-cluster291-rdkit2026036")
        precondition(Chemistry.descriptorCatalog.count == 291)
        precondition(Chemistry.descriptorCatalog.filter { $0.origin == "rdkit-standard" }.count == 217)
        guard CommandLine.arguments.count == 6 else {
            fatalError("CacheIntegration v36-cache.json v37-cache.json v38-cache.json v41-cluster-cache.json cluster-sources.json")
        }
        let directory = URL.temporaryDirectory.appendingPathComponent("chemob-cache-test-" + UUID().uuidString)
        defer { try? FileManager.default.removeItem(at: directory) }
        // A real previously accepted invalid molecule must not survive revision
        // invalidation. This fixture was captured from frozen v36, not mocked.
        try FileManager.default.createDirectory(at: directory, withIntermediateDirectories: true)
        let oldCache = try Data(contentsOf: URL(fileURLWithPath: CommandLine.arguments[1]))
        try oldCache.write(to: directory.appending(path: "descriptors.json"))
        let store = ChemistryStore(directory: directory)
        let oldFalsePositive = try await store.analyze([CompoundInput(id: 0, name: "Previously accepted", smiles: "C(C)(C)(C)(C)C")])
        precondition(oldFalsePositive.reused == 0 && !oldFalsePositive.rows[0].isValid)
        let revisionCounts = await store.statistics
        precondition(revisionCounts.computedStructures == 1)
        // A real v37 negative must also be invalidated now that the annotated
        // structure is supported. Its original source stays the cache key.
        let cxDirectory = directory.appending(path: "cx-revision")
        try FileManager.default.createDirectory(at: cxDirectory, withIntermediateDirectories: true)
        try Data(contentsOf: URL(fileURLWithPath: CommandLine.arguments[2]))
            .write(to: cxDirectory.appending(path: "descriptors.json"))
        let cxStore = ChemistryStore(directory: cxDirectory)
        let cxSource = "CC(O)F |wU:1.0|"
        let cxInputs = try CompoundImport.parse(cxSource + " First\n" + cxSource + " Again")
        precondition(cxInputs.map(\.smiles) == [cxSource, cxSource])
        let cx = try await cxStore.analyze(cxInputs)
        precondition(cx.reused == 1 && cx.rows.allSatisfy(\.isValid))
        precondition(cx.rows[0].formula == "C2H5FO" && abs(cx.rows[0].mw! - 64.059) <= 1e-8)
        let cxCounts = await cxStore.statistics
        precondition(cxCounts.computedStructures == 1)
        await cxStore.save()
        let cxReopened = ChemistryStore(directory: cxDirectory)
        let cxRestored = try await cxReopened.analyze(cxInputs)
        precondition(cxRestored.reused == 2 && cxRestored.rows.allSatisfy(\.isValid))
        let cxRestoredCounts = await cxReopened.statistics
        precondition(cxRestoredCounts.computedStructures == 0)
        let inputs = [CompoundInput(id: 20, name: "First", smiles: "CCO"), CompoundInput(id: 21, name: "Again", smiles: "CCO"), CompoundInput(id: 22, name: "Broken", smiles: "invalid")]
        let first = try await store.analyze(inputs)
        precondition(first.rows.count == 3 && first.reused == 1)
        precondition(first.rows.map(\.index) == [0,1,2])
        precondition(first.rows[0].formula == "C2H6O" && !first.rows[2].isValid)
        let second = try await store.analyze(inputs)
        precondition(second.reused == 3)
        let encoder = JSONEncoder(); encoder.outputFormatting = [.sortedKeys]
        let firstData = try encoder.encode(first.rows); let secondData = try encoder.encode(second.rows)
        precondition(firstData == secondData)
        _ = await store.graph("CCO"); _ = await store.graph("CCO")
        let counts = await store.statistics
        precondition(counts.computedStructures == 3 && counts.graphBuilds == 1)
        let correctedInputs = [CompoundInput(id: 30, name: "Nitromethane", smiles: "CN(=O)=O"), CompoundInput(id: 31, name: "Again", smiles: "CN(=O)=O")]
        let corrected = try await store.analyze(correctedInputs)
        precondition(corrected.reused == 1)
        precondition(abs(corrected.rows[0].logp! - -0.10710000000000008) <= 1e-9)
        precondition(abs(corrected.rows[0].tpsa! - 43.14) <= 1e-9)
        precondition(corrected.rows[0].descriptor_values?.count == Chemistry.descriptorCatalog.count)
        let cerium = "O=N1=[O][Ce-2]23456([O]1)([O]N(=O)=[O]2)([O]N(=O)=[O]3)([O]N(=O)=[O]4)([O]N(=O)=[O]5)[O]N(=O)=[O]6"
        let coordinatedInputs = [cerium, cerium, "[OH2][Mg+2]", "[OH2][Mg+2]"].enumerated().map {
            CompoundInput(id: 40 + $0.offset, name: "Coordination", smiles: $0.element)
        }
        let coordinated = try await store.analyze(coordinatedInputs)
        precondition(coordinated.reused == 2)
        let indices = Dictionary(uniqueKeysWithValues: Chemistry.descriptorCatalog.enumerated().map { ($0.element.id, $0.offset) })
        func value(_ row: AnalysisRow, _ id: String) -> Double? { row.descriptor_values![indices[id]!] }
        // Recompute actual old negative entries when geometry/domain support
        // changes, then reuse every descriptor across duplicates and disk load.
        let largeDirectory = directory.appending(path: "geometry-large-revision")
        try FileManager.default.createDirectory(at: largeDirectory, withIntermediateDirectories: true)
        try Data(contentsOf: URL(fileURLWithPath: CommandLine.arguments[3]))
            .write(to: largeDirectory.appending(path: "descriptors.json"))
        let largeSources = ["F[Pt@SP1](Cl)(Br)I", String(repeating: "C", count: 257)]
        let largeInputs = (largeSources + largeSources).enumerated().map {
            CompoundInput(id: 200 + $0.offset, name: "Newly supported", smiles: $0.element)
        }
        let largeStore = ChemistryStore(directory: largeDirectory)
        let large = try await largeStore.analyze(largeInputs)
        precondition(large.reused == 2 && large.rows.allSatisfy(\.isValid))
        precondition(large.rows[0].formula == "BrClFIPt" && large.rows[1].formula == "C257H516")
        precondition(abs(value(large.rows[0], "SPS")! - 15.2) < 1e-9)
        precondition(value(large.rows[0], "NumUnspecifiedAtomStereoCenters") == 0)
        precondition(value(large.rows[0], "NumAtomStereoCenters") == 1)
        precondition(value(large.rows[1], "ds_Wiener") == 2_829_056)
        let largeCounts = await largeStore.statistics
        precondition(largeCounts.computedStructures == 2)
        await largeStore.save()
        let largeReopened = ChemistryStore(directory: largeDirectory)
        let largeRestored = try await largeReopened.analyze(largeInputs)
        precondition(largeRestored.reused == 4)
        let largeData = try encoder.encode(large.rows)
        let largeRestoredData = try encoder.encode(largeRestored.rows)
        precondition(largeRestoredData == largeData)
        let restoredLargeCounts = await largeReopened.statistics
        precondition(restoredLargeCounts.computedStructures == 0)
        // Five actual frozen-v41 negative entries must be recomputed, then
        // reused across duplicate input and disk reopen with all291 intact.
        struct ClusterSource: Decodable { let name: String; let smiles: String }
        let clusterSources = try JSONDecoder().decode([ClusterSource].self, from:
            Data(contentsOf: URL(fileURLWithPath: CommandLine.arguments[5])))
        precondition(clusterSources.count == 5)
        let clusterDirectory = directory.appending(path: "cluster-revision")
        try FileManager.default.createDirectory(at: clusterDirectory, withIntermediateDirectories: true)
        try Data(contentsOf: URL(fileURLWithPath: CommandLine.arguments[4]))
            .write(to: clusterDirectory.appending(path: "descriptors.json"))
        let clusterInputs = (clusterSources + clusterSources).enumerated().map {
            CompoundInput(id: 300 + $0.offset, name: $0.element.name, smiles: $0.element.smiles)
        }
        let clusterStore = ChemistryStore(directory: clusterDirectory)
        let clusters = try await clusterStore.analyze(clusterInputs)
        precondition(clusters.reused == 5 && clusters.rows.allSatisfy(\.isValid))
        precondition(Array(clusters.rows.prefix(5)).map(\.formula) == ["H9Fe3O12", "Fe3S4-2", "H9Fe3O12", "Au2Cl6", "Br12Ta6"])
        precondition(abs(value(clusters.rows[1], "SPS")! - 40.285714285714285) < 1e-9)
        precondition(value(clusters.rows[4], "RingCount") == 18)
        let clusterCounts = await clusterStore.statistics
        precondition(clusterCounts.computedStructures == 5)
        await clusterStore.save()
        let clusterReopened = ChemistryStore(directory: clusterDirectory)
        let clusterRestored = try await clusterReopened.analyze(clusterInputs)
        precondition(clusterRestored.reused == 10)
        let clusterData = try encoder.encode(clusters.rows)
        let clusterRestoredData = try encoder.encode(clusterRestored.rows)
        precondition(clusterData == clusterRestoredData)
        let clusterRestoredCounts = await clusterReopened.statistics
        precondition(clusterRestoredCounts.computedStructures == 0)
        precondition(coordinated.rows[0].hba == 18 && coordinated.rows[0].rotatable_bonds == 12)
        precondition(value(coordinated.rows[0], "RingCount") == 0 && value(coordinated.rows[0], "ds_SingleBonds") == 18)
        precondition(abs(value(coordinated.rows[0], "BertzCT")! - 83617.07236829486) <= 1e-9)
        for descriptor in Chemistry.descriptorCatalog where descriptor.id.hasPrefix("BCUT2D_") {
            precondition(value(coordinated.rows[2], descriptor.id) == nil, "unsupported BCUT bond type must stay undefined")
        }
        let aromaticInputs = ["[n+]1ccccc1", "[n+]1ccccc1", "[NH3+]"].enumerated().map {
            CompoundInput(id: 50 + $0.offset, name: "Radical", smiles: $0.element)
        }
        let aromatic = try await store.analyze(aromaticInputs)
        precondition(aromatic.reused == 1)
        precondition(aromatic.rows[0].hba == 1 && abs(aromatic.rows[0].logp! - 0.47650000000000003) <= 1e-9)
        for (id, expected) in [("FpDensityMorgan1", 1.0), ("FpDensityMorgan2", 2.0), ("FpDensityMorgan3", 2.1666666666666665)] {
            precondition(abs(value(aromatic.rows[0], id)! - expected) <= 1e-9)
        }
        precondition(aromatic.rows[2].hba == 1 && aromatic.rows[2].hbd == 1)
        let macroSource = "O=C1OC[C@H]2OC(=O)C3=C(OC(=O)C4=CC(O)=C(O)C(O)=C4C4=C(O)C(O)=C(O)C=C4C(=O)O3)[C@@H]2OC(=O)C2=CC(O)=C(O)C(O)=C2C2=C(O)C(O)=C(O)C=C12"
        let macroInputs = [macroSource, macroSource].enumerated().map {
            CompoundInput(id: 60 + $0.offset, name: "Macrocycle", smiles: $0.element)
        }
        let macro = try await store.analyze(macroInputs)
        precondition(macro.reused == 1)
        precondition(abs(value(macro.rows[0], "FpDensityMorgan1")! - 0.4107142857142857) <= 1e-9)
        precondition(abs(macro.rows[0].logp! - 1.3480000000000005) <= 1e-9)
        let qedInputs = ["c1cc[asH]c1", "[H]OOC", "[Hg]", "c1cc[se]c1", "c1cc[asH]c1", "c1cc[te]c1"].enumerated().map {
            CompoundInput(id: 70 + $0.offset, name: "QED", smiles: $0.element)
        }
        let qed = try await store.analyze(qedInputs)
        precondition(qed.reused == 1)
        precondition(abs(qed.rows[0].logp! - 0.4638999999999999) <= 1e-9, "As aromaticity correction must reach the cache")
        for (row, expected) in [(0, 0.41560724564829626), (1, 0.31061534304923694), (2, 0.48313597685210746), (5, 0.5159286568196739)] {
            precondition(abs(value(qed.rows[row], "qed")! - expected) <= 1e-9)
        }
        precondition(value(qed.rows[3], "NumRadicalElectrons") == 0)
        let valenceInputs = ["C:C", "C:C", "C1:C:C:C:C:C#1", "C1=CC#CC=C1", "C:S:C"].enumerated().map {
            CompoundInput(id: 80 + $0.offset, name: "Valence", smiles: $0.element)
        }
        let valence = try await store.analyze(valenceInputs)
        precondition(valence.reused == 1)
        precondition(valence.rows[0].formula == "C2H4", "partial aromatic bond valence must reach the cache")
        precondition(value(valence.rows[0], "ds_TotalHydrogens") == 4)
        for row in [valence.rows[2], valence.rows[3]] {
            precondition(row.formula == "C6H4")
            precondition(value(row, "NumAromaticRings") == 1 && value(row, "ds_TripleBonds") == 1)
            precondition(abs(value(row, "BertzCT")! - 71.96100505779535) <= 1e-9)
        }
        precondition(Chemistry.descriptorCatalog.filter { $0.origin == "native-extra" }.count == 50,
                     "the fifty added descriptors must reach the actual Swift catalog")
        let topologyInputs = ["CCCC", "CCCC", "[H][H]", "[NH3]->[Cu+2]<-[NH3]"].enumerated().map {
            CompoundInput(id: 90 + $0.offset, name: "Topology", smiles: $0.element)
        }
        let topology = try await store.analyze(topologyInputs)
        precondition(topology.reused == 1 && topology.rows.allSatisfy(\.isValid))
        for (id, expected) in [("ds_Wiener", 10.0), ("ds_Harary", 13.0/3.0),
                               ("ds_EccentricConnectivity", 14.0), ("ds_DegreeDistance", 28.0),
                               ("ds_Gutman", 19.0), ("ds_MDEC11", 1.0/3.0),
                               ("ds_MWC3", 2.833213344056216), ("ds_SRW3", 0.0)] {
            precondition(abs(value(topology.rows[0], id)! - expected) <= 1e-9)
        }
        for spec in Chemistry.descriptorCatalog where spec.origin == "native-extra" {
            precondition(value(topology.rows[2], spec.id) == nil, "a hydrogen-only graph has undefined extra50")
        }
        precondition(value(topology.rows[3], "ds_Wiener") == 4)
        guard indices["SPS"] != nil && indices["NumAtomStereoCenters"] != nil && indices["NumUnspecifiedAtomStereoCenters"] != nil else {
            FileHandle.standardError.write(Data("FAIL: all three newly added standards must reach the actual Swift catalog\n".utf8))
            preconditionFailure("missing standard descriptors")
        }
        let stereoInputs = ["CC(O)Cl", "C[C@H](O)Cl", "CC1CCC(C)CC1", "C[C@H](C)O", "[H][H]", "C[C@H](O)Cl", "[H]C([H])([H])[H]"].enumerated().map {
            CompoundInput(id: 100 + $0.offset, name: "Stereo", smiles: $0.element)
        }
        let stereo = try await store.analyze(stereoInputs)
        precondition(stereo.reused == 2 && stereo.rows.allSatisfy(\.isValid))
        for (index,total,unspecified,sps) in [(0,1.0,1.0,15.75),(1,1.0,0.0,15.75),(2,0.0,0.0,39.75),(3,0.0,0.0,9.0),(6,0.0,0.0,64.0)] {
            precondition(value(stereo.rows[index], "NumAtomStereoCenters") == total)
            precondition(value(stereo.rows[index], "NumUnspecifiedAtomStereoCenters") == unspecified)
            precondition(value(stereo.rows[index], "SPS") == sps)
        }
        precondition(value(stereo.rows[4], "SPS") == nil, "SPS no-heavy domain must remain null")
        precondition(value(stereo.rows[3], "ds_ExplicitHydrogens") == 0 && value(stereo.rows[3], "ds_ImplicitHydrogens") == 8)
        let donorInputs = ["B1OBOBO1", "B1OBOBO1", "B1SBSBS1",
            "C[C@@H]1C/C=C/CC2=CC=C(CC3=C(C(=O)O1)C(O)=CC(O)=C3Cl)O2"].enumerated().map {
            CompoundInput(id: 120 + $0.offset, name: "Aromatic preparation", smiles: $0.element)
        }
        let donor = try await store.analyze(donorInputs)
        precondition(donor.reused == 1 && donor.rows.allSatisfy(\.isValid))
        for (index,expected) in [(0,20.0),(1,20.0),(2,20.0),(3,21.916666666666668)] {
            guard value(donor.rows[index], "SPS").map({abs($0-expected)<=1e-9}) == true else {
                FileHandle.standardError.write(Data("FAIL: corrected aromatic preparation must reach Swift cache\n".utf8))
                preconditionFailure("stale aromatic preparation")
            }
        }
        precondition(value(donor.rows[0], "NumAromaticRings") == 0 && value(donor.rows[2], "NumAromaticRings") == 0)
        // Independent Decimal85 molecular controls: keep exact graph definitions
        // through the real Swift store, then through its on-disk representation.
        let informationInputs = ["CCO", "[HH]", "[nH]1cccccc1", String(repeating: "C", count: 128), "[nH]1cccccc1"].enumerated().map {
            CompoundInput(id: 140 + $0.offset, name: "Stable graph information", smiles: $0.element)
        }
        let information = try await store.analyze(informationInputs)
        precondition(information.reused == 2 && information.rows.allSatisfy(\.isValid))
        for (index,ipc,average) in [(0,2.7548875021634685,0.9182958340544895),
                                    (1,0.0,0.0), (2,58.974147804380216,1.9023918646574263),
                                    (3,1.5501297715622655e27,3.8058131928162322)] {
            precondition(abs(value(information.rows[index], "Ipc")! - ipc) <= max(1e-10, abs(ipc)*1e-12))
            precondition(abs(value(information.rows[index], "AvgIpc")! - average) <= max(1e-12, abs(average)*1e-12))
        }
        precondition(information.rows[1].formula == "H2" && value(information.rows[1], "ds_TotalHydrogens") == 2)
        precondition(information.rows[2].formula == "C6H7N" && value(information.rows[2], "ds_ImplicitHydrogens") == 7)
        await store.save()
        let reopened = ChemistryStore(directory: directory)
        let disk = try await reopened.analyze(inputs)
        _ = await reopened.graph("CCO")
        let restoredCounts = await reopened.statistics
        precondition(disk.reused == 3 && restoredCounts.computedStructures == 0 && restoredCounts.graphBuilds == 0)
        let correctedDisk = try await reopened.analyze(correctedInputs)
        let correctedDiskCounts = await reopened.statistics
        precondition(correctedDisk.reused == 2 && correctedDiskCounts.computedStructures == 0)
        let correctedData = try encoder.encode(corrected.rows)
        let correctedDiskData = try encoder.encode(correctedDisk.rows)
        precondition(correctedData == correctedDiskData, "corrected descriptors must survive disk reload exactly")
        let coordinatedDisk = try await reopened.analyze(coordinatedInputs)
        let coordinatedDiskCounts = await reopened.statistics
        precondition(coordinatedDisk.reused == 4 && coordinatedDiskCounts.computedStructures == 0)
        let coordinatedData = try encoder.encode(coordinated.rows)
        let coordinatedDiskData = try encoder.encode(coordinatedDisk.rows)
        precondition(coordinatedData == coordinatedDiskData, "directional chemistry and undefined values must survive disk reload")
        let aromaticDisk = try await reopened.analyze(aromaticInputs)
        let aromaticDiskCounts = await reopened.statistics
        precondition(aromaticDisk.reused == 3 && aromaticDiskCounts.computedStructures == 0)
        let aromaticData = try encoder.encode(aromatic.rows)
        let aromaticDiskData = try encoder.encode(aromaticDisk.rows)
        precondition(aromaticData == aromaticDiskData, "Morgan densities and radical preparation must survive disk reload")
        let macroDisk = try await reopened.analyze(macroInputs)
        let macroDiskCounts = await reopened.statistics
        precondition(macroDisk.reused == 2 && macroDiskCounts.computedStructures == 0)
        let macroData = try encoder.encode(macro.rows)
        let macroDiskData = try encoder.encode(macroDisk.rows)
        precondition(macroData == macroDiskData, "corrected macrocycle properties must survive disk reload")
        let qedDisk = try await reopened.analyze(qedInputs)
        let qedDiskCounts = await reopened.statistics
        precondition(qedDisk.reused == 6 && qedDiskCounts.computedStructures == 0)
        let qedData = try encoder.encode(qed.rows)
        let qedDiskData = try encoder.encode(qedDisk.rows)
        precondition(qedData == qedDiskData, "QED and corrected As/Se values must survive disk reload")
        let valenceDisk = try await reopened.analyze(valenceInputs)
        let valenceDiskCounts = await reopened.statistics
        precondition(valenceDisk.reused == 5 && valenceDiskCounts.computedStructures == 0)
        let valenceData = try encoder.encode(valence.rows)
        let valenceDiskData = try encoder.encode(valenceDisk.rows)
        precondition(valenceData == valenceDiskData, "partial aromatic valence and triple aromatic flags must survive disk reload")
        let topologyDisk = try await reopened.analyze(topologyInputs)
        let topologyDiskCounts = await reopened.statistics
        precondition(topologyDisk.reused == 4 && topologyDiskCounts.computedStructures == 0)
        let topologyData = try encoder.encode(topology.rows)
        let topologyDiskData = try encoder.encode(topologyDisk.rows)
        precondition(topologyData == topologyDiskData, "all50 added fields, nulls and directed inputs must survive disk reload")
        let stereoDisk = try await reopened.analyze(stereoInputs)
        let stereoDiskCounts = await reopened.statistics
        precondition(stereoDisk.reused == stereoInputs.count && stereoDiskCounts.computedStructures == 0)
        let stereoData = try encoder.encode(stereo.rows)
        let stereoDiskData = try encoder.encode(stereoDisk.rows)
        precondition(stereoData == stereoDiskData, "legacy counts, modern SPS, nulls and H provenance must survive disk reload")
        let donorDisk = try await reopened.analyze(donorInputs)
        let donorDiskCounts = await reopened.statistics
        precondition(donorDisk.reused == donorInputs.count && donorDiskCounts.computedStructures == 0)
        let donorData = try encoder.encode(donor.rows)
        let donorDiskData = try encoder.encode(donorDisk.rows)
        precondition(donorData == donorDiskData, "corrected aromatic donors must survive disk reload")
        let informationDisk = try await reopened.analyze(informationInputs)
        let informationDiskCounts = await reopened.statistics
        precondition(informationDisk.reused == informationInputs.count && informationDiskCounts.computedStructures == 0)
        let informationData = try encoder.encode(information.rows)
        let informationDiskData = try encoder.encode(informationDisk.rows)
        precondition(informationData == informationDiskData, "stable information and corrected H/Kekulize provenance must survive disk reload")
        // Different spellings, isotope, stereo and charge must not share a result.
        let distinct = ["OCC", "[13CH4]", "[NH4+]", "N[C@H](C)O", "N[C@@H](C)O"]
        _ = try await reopened.analyze(distinct.enumerated().map { CompoundInput(id:$0.offset,name:"distinct",smiles:$0.element) })
        let distinctCounts = await reopened.statistics
        precondition(distinctCounts.computedStructures == 5)
        for file in ["descriptors.json", "graphs.json"] {
            let url = directory.appendingPathComponent(file)
            var object = try JSONSerialization.jsonObject(with: Data(contentsOf: url)) as! [String: Any]
            object["namespace"] = "previous-engine"
            try JSONSerialization.data(withJSONObject: object).write(to: url)
        }
        let invalidated = ChemistryStore(directory: directory)
        _ = try await invalidated.analyze(inputs)
        _ = await invalidated.graph("CCO")
        let refreshed = await invalidated.statistics
        precondition(refreshed.computedStructures == 2 && refreshed.graphBuilds == 1)
        print("PASS duplicate inputs, corrected chemistry/coordination/Morgan/radical/QED/valence/undefined cache, invalid-input cache, row positions, graph reuse, disk reload, stale-engine invalidation, charge/isotope/stereo separation")
    }
}
