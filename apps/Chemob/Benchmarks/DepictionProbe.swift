import Foundation

/// Compile alongside Core/.../ChemicalDepiction.swift; no App/UI state.
@main struct DepictionProbe {
    struct Record: Decodable { let name: String; let smiles: String; let graph: MolecularDrawingGraph? }
    static func main() throws {
        let directory = URL(fileURLWithPath: CommandLine.arguments[1])
        try FileManager.default.createDirectory(at: directory,withIntermediateDirectories: true)
        var count=0,available=0,graphUnavailable=0,qualityUnavailable=0,stereo=0
        var failures=[[String:String]](),curated=[[String:Any]]()
        var partitions=[String:[String:Int]](),latencies=[Double](),timings=[[String:Any]]()
        func bump(_ group:String,_ metric:String) { var p=partitions[group,default:[:]];p[metric,default:0] += 1;partitions[group]=p }
        let proofURL=directory.appending(path:"stereo-drawings.jsonl")
        _ = FileManager.default.createFile(atPath: proofURL.path,contents: nil)
        let proof=try FileHandle(forWritingTo: proofURL)
        defer { try? proof.close() }
        let clock=ContinuousClock(),start=clock.now
        while let line=readLine() {
            let record=try JSONDecoder().decode(Record.self,from: Data(line.utf8)); count += 1
            let group=record.name.hasPrefix("curated-") ? "curated" : "original"
            bump(group,"inputs")
            guard let graph=record.graph else { graphUnavailable += 1;bump(group,"native_graph_unavailable");continue }
            let t=clock.now
            let value=ChemicalDepiction.make(graph),elapsed=t.duration(to:clock.now)
            let milliseconds=Double(elapsed.components.seconds)*1000+Double(elapsed.components.attoseconds)/1e15
            latencies.append(milliseconds);timings.append(["name":record.name,"milliseconds":milliseconds,"atoms":graph.atoms.count,"available":value != nil])
            guard let drawing=value else {
                qualityUnavailable += 1
                bump(group,"quality_unavailable")
                failures.append(["name":record.name,"smiles":record.smiles]);continue
            }
            bump(group,"available")
            available += 1; if drawing.stereochemistryOmitted { stereo += 1 }
            if drawing.stereochemistryOmitted { bump(group,"stereo_partly_omitted") }
            if !drawing.wedges.isEmpty || !drawing.alkenes.isEmpty {
                struct Proof: Encodable { let smiles: String;let drawing: ChemicalDepiction;let screenPoints:[DrawingPoint] }
                try proof.write(contentsOf: JSONEncoder().encode(Proof(smiles: record.smiles,drawing: drawing,screenPoints:drawing.points.map(drawing.screenPoint))))
                try proof.write(contentsOf: Data([10]))
            }
            if record.name.hasPrefix("curated-") {
                try drawing.svg().write(to: directory.appending(path: record.name+".svg"),atomically: true,encoding: .utf8)
                let data=try JSONEncoder().encode(drawing)
                try data.write(to: directory.appending(path: record.name+".json"))
                curated.append(["name":record.name,"milliseconds":Double(t.duration(to: clock.now).components.attoseconds)/1e15])
            }
            if count%1000==0 { try FileHandle.standardError.write(contentsOf:Data("Processed \(count) inputs\n".utf8)) }
        }
        let d=start.duration(to: clock.now)
        let sorted=latencies.sorted()
        let summary:[String:Any] = ["inputs":count,"available":available,"native_graph_unavailable":graphUnavailable,"partitions":partitions,
            "layout_only_milliseconds":["median":sorted[sorted.count/2],"p95":sorted[min(sorted.count-1,Int(Double(sorted.count)*0.95))],"max":sorted.last!],
            "slowest_layouts":timings.sorted { ($0["milliseconds"] as! Double)>($1["milliseconds"] as! Double) }.prefix(10).map{$0},
            "timing_scope":"Single host sweep, pure layout excludes decode/SVG/cache; no physical device claim or controlled other-host workload",
            "quality_unavailable":qualityUnavailable,"stereo_omitted_and_flagged":stereo,"seconds_including_decode_and_fixture_output":Double(d.components.seconds)+Double(d.components.attoseconds)/1e18,"curated":curated]
        try JSONSerialization.data(withJSONObject: summary,options:[.prettyPrinted,.sortedKeys]).write(to:directory.appending(path:"summary.json"))
        try JSONSerialization.data(withJSONObject:failures,options:[.prettyPrinted,.sortedKeys]).write(to:directory.appending(path:"unavailable.json"))
        print(summary)
    }
}
