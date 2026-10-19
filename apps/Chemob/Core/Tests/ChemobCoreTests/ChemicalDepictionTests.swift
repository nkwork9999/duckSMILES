import XCTest
#if canImport(FoundationXML)
import FoundationXML
#endif
@testable import ChemobCore

final class ChemicalDepictionTests: XCTestCase {
    func testCoordinationArrowPointsTowardTheAcceptor() throws {
        let graph=MolecularDrawingGraph(atoms:[.init(index:0,symbol:"N"),.init(index:1,symbol:"Cu")],bonds:[.init(a:0,b:1,kind:"dative")])
        let drawing=try XCTUnwrap(ChemicalDepiction.make(graph))
        XCTAssertEqual(drawing.strokes.count,3,"Coordination needs a shaft and an actual arrowhead")
        if drawing.strokes.count==3 {
            let head=Array(drawing.strokes.suffix(2)),acceptor=drawing.points[1]
            XCTAssertEqual(head[0].b,head[1].b)
            for edge in head { XCTAssertLessThan(edge.b.distance(acceptor),edge.a.distance(acceptor)) }
            XCTAssertNotEqual(head[0].a,head[1].a)
        }
    }

    func testPreparationCannotSilentlyLoseNonTetrahedralStereo() throws {
        let graph=MolecularDrawingGraph(atoms:[.init(index:0,symbol:"Pt"),.init(index:1,symbol:"Cl"),.init(index:2,symbol:"F"),.init(index:3,symbol:"Br"),.init(index:4,symbol:"I")],bonds:[1,2,3,4].map{.init(a:0,b:$0,kind:"single")},unrepresented_stereo:true)
        let drawing=try XCTUnwrap(ChemicalDepiction.make(graph))
        XCTAssertTrue(drawing.stereochemistryOmitted,"Prepared graph's omitted geometry must stay visible to the user")
    }

    func testSlashAndBackslashProduceOppositeDoubleBondConfigurations() throws {
        for (direction,wantSameSide) in [(1,false),(-1,true)] {
            let graph=MolecularDrawingGraph(atoms:(0..<4).map{.init(index:$0,symbol:($0==0 || $0==3) ? "F":"C")},
                bonds:[.init(a:0,b:1,kind:"single",direction:1),.init(a:1,b:2,kind:"double"),.init(a:2,b:3,kind:"single",direction:direction)])
            let drawing=try XCTUnwrap(ChemicalDepiction.make(graph)),p=drawing.points
            let axis=p[2]-p[1],left=p[0]-p[1],right=p[3]-p[2]
            let product=(axis.x*left.y-axis.y*left.x)*(axis.x*right.y-axis.y*right.x)
            XCTAssertEqual(product>0,wantSameSide,"Slash parity must change geometric cis/trans")
            XCTAssertEqual(drawing.alkenes.count,1)
            XCTAssertFalse(drawing.stereochemistryOmitted)
        }
    }

    func testStraightCarbonVerticesAreNeverInvisibleInSkeletalNotation() throws {
        let graph=MolecularDrawingGraph.carbon(edges:[(0,1),(1,2),(2,3),(3,4),(4,5),(5,0),(5,6),(6,2)])
        let drawing=try XCTUnwrap(ChemicalDepiction.make(graph))
        for i in graph.atoms.indices {
            let neighbors=graph.bonds.compactMap { $0.a==i ? $0.b : $0.b==i ? $0.a : nil }
            guard neighbors.count==2,!drawing.labels.contains(where:{$0.point==drawing.points[i]}) else { continue }
            let a=drawing.points[neighbors[0]]-drawing.points[i],b=drawing.points[neighbors[1]]-drawing.points[i]
            XCTAssertGreaterThan(abs(a.x*b.y-a.y*b.x),0.02,"An unlabelled straight carbon would disappear into one long bond")
        }
    }

    func testSVGUsesUpwardMolecularYForCorrectWedgeHandedness() throws {
        let origin=DrawingPoint(x:0,y:0)
        let drawing=ChemicalDepiction(points:[origin,.init(x:1,y:0),.init(x:0,y:1)],strokes:[.init(a:origin,b:.init(x:1,y:0),dashed:false),.init(a:origin,b:.init(x:0,y:1),dashed:false)],labels:[],wedges:[],alkenes:[],stereochemistryOmitted:false)
        let probe=SVGStrokeProbe(),parser=XMLParser(data:Data(drawing.svg().utf8))
        parser.delegate=probe;XCTAssertTrue(parser.parse())
        XCTAssertEqual(probe.endpoints.count,2)
        XCTAssertLessThan(probe.endpoints[1].y,probe.endpoints[0].y,"Mirroring Y without reversing wedges changes a molecule's handedness")
    }

    func testBridgedBicyclicGraphCanUsePlanarGeometryWithoutCrossings() throws {
        // Norbornane topology: a regular outer hexagon plus its center provides
        // an independent crossing-free unit-bond embedding of this graph.
        let graph=MolecularDrawingGraph.carbon(edges: [(0,1),(1,2),(2,3),(3,4),(4,5),(5,0),(5,6),(6,2)])
        let drawing=try XCTUnwrap(ChemicalDepiction.make(graph))
        XCTAssertEqual(drawing.points.count,7)
        for b in graph.bonds { XCTAssertTrue((0.85...1.15).contains(drawing.points[b.a].distance(drawing.points[b.b]))) }
    }

    func testBranchedTreeDoesNotCrossItsExistingBranches() throws {
        // Original3431965 topology, elements/stereo suppressed for this pure
        // geometric regression. Old greedy placement crossed existing bonds.
        let edges=[(0,1),(1,2),(2,3),(2,4),(4,5),(4,6),(6,7),(6,8),(8,9),(9,10),(10,11),(11,12),(11,13),(9,14),(14,15),(14,16),(16,17),(17,18),(17,19),(19,20),(19,21),(21,22),(21,23),(23,24),(24,25),(24,26),(26,27),(26,28),(28,29),(28,30),(16,31),(31,32),(32,33),(32,34),(31,35),(35,36),(35,37)]
        XCTAssertNotNil(ChemicalDepiction.make(.carbon(edges: edges)))
    }

    func testTetrahedralEnantiomersUseOppositeWedgesAndKeepOtherBonds() throws {
        var drawings=[ChemicalDepiction]()
        for token in ["@","@@"] {
            let graph=MolecularDrawingGraph(atoms: [.init(index: 0,symbol: "N"),
                .init(index: 1,symbol: "C",bracket_h_count: 1,chirality_token: token,stereo_neighbor_order: [0,-1,2,3]),
                .init(index: 2,symbol: "C"),.init(index: 3,symbol: "O")],
                bonds: [.init(a: 0,b: 1,kind: "single"),.init(a: 1,b: 2,kind: "single"),.init(a: 1,b: 3,kind: "single")])
            let drawing=try XCTUnwrap(ChemicalDepiction.make(graph))
            XCTAssertEqual(drawing.wedges.count,1,"A specified tetrahedral center needs a real wedge")
            XCTAssertEqual(drawing.strokes.count,2,"Wedge replaces exactly one single bond")
            XCTAssertFalse(drawing.stereochemistryOmitted)
            drawings.append(drawing)
        }
        if let a=drawings[0].wedges.first,let b=drawings[1].wedges.first {
            XCTAssertEqual(a.center,b.center);XCTAssertEqual(a.neighbor,b.neighbor)
            XCTAssertNotEqual(a.hashed,b.hashed,"Enantiomer change must reverse the depicted configuration")
        }
    }

    // Catches the old circular placement: a chain must have 120-degree bends
    // and all bonds one unit long, rather than closing into a polygon.
    func testChainHasUniformBondsAndZigzagAngles() throws {
        let graph = MolecularDrawingGraph.carbon(edges: [(0,1),(1,2),(2,3),(3,4)])
        let drawing = try XCTUnwrap(ChemicalDepiction.make(graph))
        for bond in graph.bonds { XCTAssertEqual(drawing.points[bond.a].distance(drawing.points[bond.b]), 1, accuracy: 0.01) }
        for i in 1...3 {
            let a = drawing.points[i-1] - drawing.points[i], b = drawing.points[i+1] - drawing.points[i]
            XCTAssertEqual((a.x*b.x+a.y*b.y)/(a.length*b.length), -0.5, accuracy: 0.01)
        }
        XCTAssertGreaterThan(drawing.points[0].distance(drawing.points[4]), 2.9)
    }

    // Rejects ring collapse, wrong fused-ring side, or overlapping atoms.
    func testRegularBenzeneAndFusedRings() throws {
        for edges in [[(0,1),(1,2),(2,3),(3,4),(4,5),(5,0)],
                      [(0,1),(1,2),(2,3),(3,4),(4,5),(5,0),(2,6),(6,7),(7,8),(8,9),(9,3)]] {
            let graph = MolecularDrawingGraph.carbon(edges: edges)
            let drawing = try XCTUnwrap(ChemicalDepiction.make(graph))
            for bond in graph.bonds { XCTAssertEqual(drawing.points[bond.a].distance(drawing.points[bond.b]), 1, accuracy: 0.02) }
            for i in drawing.points.indices { for j in 0..<i { XCTAssertGreaterThan(drawing.points[i].distance(drawing.points[j]), 0.65) } }
        }
    }

    func testTripleBondIsLinearAndDrawsThreeSeparateStrokes() throws {
        let graph = MolecularDrawingGraph(atoms: (0..<3).map { .init(index: $0, symbol: $0 == 2 ? "N" : "C") },
            bonds: [.init(a: 0,b: 1,kind: "single"), .init(a: 1,b: 2,kind: "triple")])
        let drawing = try XCTUnwrap(ChemicalDepiction.make(graph))
        XCTAssertEqual(drawing.points[0].distance(drawing.points[2]), 2, accuracy: 0.01)
        XCTAssertEqual(drawing.strokes.count, 4)
    }

    func testDisconnectedIonsKeepChargeAndIsotopeAndDoNotGainBond() throws {
        let graph = MolecularDrawingGraph(atoms: [.init(index: 0,symbol: "Na",formal_charge: 1),
            .init(index: 1,symbol: "Cl",formal_charge: -1), .init(index: 2,symbol: "C",isotope: 13,implicit_h_count: 4)], bonds: [])
        let drawing = try XCTUnwrap(ChemicalDepiction.make(graph))
        XCTAssertTrue(drawing.strokes.isEmpty)
        XCTAssertEqual(drawing.labels.map(\.text), ["Na+", "Cl−", "13CH4"])
        XCTAssertGreaterThan(drawing.points[0].distance(drawing.points[1]), 1.5)
    }

    func testBondOrdersAromaticityAndAttachedHydrogensRemainVisible() throws {
        let graph = MolecularDrawingGraph(atoms: [.init(index: 0,symbol: "N",aromatic: true,bracket_h_count: 1),
            .init(index: 1,symbol: "C"),.init(index: 2,symbol: "O")],
            bonds: [.init(a: 0,b: 1,kind: "aromatic"),.init(a: 1,b: 2,kind: "double")])
        let drawing = try XCTUnwrap(ChemicalDepiction.make(graph))
        XCTAssertEqual(drawing.strokes.count, 4)
        XCTAssertEqual(drawing.strokes.filter(\.dashed).count, 1)
        XCTAssertEqual(drawing.labels.map(\.text), ["NH", "O"])
    }

    func testStereoIsNeverSilentlyPresentedAsSpecified2DGeometry() throws {
        let graph = MolecularDrawingGraph(atoms: [.init(index: 0,symbol: "C",chirality_token: "@@"),.init(index: 1,symbol: "F")],bonds: [.init(a: 0,b: 1,kind: "single")])
        let drawing = try XCTUnwrap(ChemicalDepiction.make(graph))
        XCTAssertTrue(drawing.stereochemistryOmitted)
        XCTAssertFalse(drawing.svg().contains("NaN"))
    }

    func testInvalidGraphAndUnrepresentableGeometryAreUnavailable() {
        XCTAssertNil(ChemicalDepiction.make(.init(atoms: [.init(index: 0,symbol: "C")],bonds: [.init(a: 0,b: 5,kind: "single")])))
        XCTAssertNil(ChemicalDepiction.make(.init(atoms: [],bonds: [])))
        // K5 cannot have a crossing-free chemical drawing.
        let edges = (0..<5).flatMap { i in ((i+1)..<5).map { (i,$0) } }
        XCTAssertNil(ChemicalDepiction.make(.carbon(edges: edges)))
    }
}

private final class SVGStrokeProbe: NSObject,XMLParserDelegate {
    var scale=[DrawingPoint(x:1,y:1)],endpoints=[DrawingPoint]()
    func parser(_ parser:XMLParser,didStartElement name:String,namespaceURI:String?,qualifiedName:String?,attributes:[String:String]) {
        if name=="g" {
            var value=scale.last!
            if let t=attributes["transform"],t.hasPrefix("scale(") {
                let numbers=t.dropFirst(6).dropLast().split(whereSeparator:{$0=="," || $0==" "}).compactMap{Double($0)}
                if numbers.count==2 { value.x *= numbers[0];value.y *= numbers[1] }
            }
            scale.append(value)
        }
        if name=="path",let d=attributes["d"],d.hasPrefix("M") {
            let numbers=d.split(whereSeparator:{"MLZ ,".contains($0)}).compactMap{Double($0)}
            if numbers.count==4 { endpoints.append(.init(x:numbers[2]*scale.last!.x,y:numbers[3]*scale.last!.y)) }
        }
    }
    func parser(_ parser:XMLParser,didEndElement name:String,namespaceURI:String?,qualifiedName:String?) { if name=="g" { scale.removeLast() } }
}

private extension MolecularDrawingGraph {
    static func carbon(edges: [(Int, Int)]) -> Self {
        .init(atoms: (0...(edges.flatMap { [$0.0,$0.1] }.max() ?? 0)).map { .init(index: $0,symbol: "C") },
              bonds: edges.map { .init(a: $0.0,b: $0.1,kind: "single") })
    }
}
