import Foundation
import ChemobCore

@main struct DepictionCacheProbe {
    static func main() async throws {
        let directory=URL(fileURLWithPath: CommandLine.arguments[1])
        let store=ChemistryStore(directory: directory)
        let encoder=JSONEncoder(); encoder.outputFormatting=[.sortedKeys]
        let first=await store.drawing("CCO")
        guard let first else { fatalError("Ethanol drawing unavailable") }
        precondition(first.labels.map(\.text)==["OH"],"Drawing must retain implicit hydrogen metadata")
        let firstData=try encoder.encode(first)
        let second=await store.drawing("CCO")
        let secondData=try encoder.encode(second)
        precondition(firstData==secondData,"Duplicate drawing output changed")
        let counts=await store.statistics
        precondition(counts.graphBuilds==1 && counts.drawingBuilds==1,"Duplicate must reuse both graph and drawing")
        let reopened=ChemistryStore(directory: directory)
        let disk=await reopened.drawing("CCO")
        let diskData=try encoder.encode(disk)
        let diskCounts=await reopened.statistics
        precondition(firstData==diskData && diskCounts.graphBuilds==0 && diskCounts.drawingBuilds==0,"Disk reopen must reuse geometry without native work")
        _ = await reopened.drawing("invalid")
        _ = await reopened.drawing("invalid")
        let invalidCounts=await reopened.statistics
        precondition(invalidCounts.graphBuilds==1 && invalidCounts.drawingBuilds==1,"Unavailable drawings must also be reused")
        let chiral=await reopened.drawing("N[C@@H](C)C(=O)O")
        precondition(chiral?.stereochemistryOmitted==false && chiral?.wedges.count==1,"Tetrahedral stereo must have a real wedge")
        let alkene=await reopened.drawing("F/C=C/F")
        precondition(alkene?.stereochemistryOmitted==false && alkene?.alkenes.count==1,"E/Z geometry must be represented and retained")
        let stereoReopen=ChemistryStore(directory: directory)
        let stereoDisk=await stereoReopen.drawing("N[C@@H](C)C(=O)O")
        let originalStereo=try encoder.encode(chiral), diskStereo=try encoder.encode(stereoDisk)
        let stereoCounts=await stereoReopen.statistics
        precondition(originalStereo==diskStereo && stereoCounts.graphBuilds==0 && stereoCounts.drawingBuilds==0,"Wedge orientation and stereo metadata must survive disk reuse")
        print("PASS actual native graph→2D, OH label, duplicate/disk reuse without calculations, unavailable-result reuse, stereo notice")
    }
}
