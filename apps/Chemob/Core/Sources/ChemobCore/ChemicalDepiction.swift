import Foundation

/// Versioned drawing input from the descriptor preparation pipeline. Retains
/// input atom order, attached hydrogen counts and supported stereo annotations.
/// Unsupported configurations are explicitly flagged by the native producer.
public struct MolecularDrawingGraph: Codable, Sendable {
    public struct Atom: Codable, Sendable {
        public let index: Int
        public let symbol: String
        public let formal_charge: Int
        public let isotope: Int?
        public let aromatic: Bool?
        public let bracket_h_count: Int?
        public let implicit_h_count: Int?
        public let chirality_token: String?
        public let stereo_neighbor_order: [Int]?
        public init(index: Int, symbol: String, formal_charge: Int = 0, isotope: Int? = nil, aromatic: Bool? = nil, bracket_h_count: Int? = nil, implicit_h_count: Int? = nil, chirality_token: String? = nil, stereo_neighbor_order: [Int]? = nil) {
            self.index=index; self.symbol=symbol; self.formal_charge=formal_charge; self.isotope=isotope
            self.aromatic=aromatic; self.bracket_h_count=bracket_h_count; self.implicit_h_count=implicit_h_count; self.chirality_token=chirality_token
            self.stereo_neighbor_order=stereo_neighbor_order
        }
    }
    public struct Bond: Codable, Sendable {
        public let a: Int; public let b: Int; public let kind: String; public let direction: Int?
        enum CodingKeys: String, CodingKey { case a="begin", b="end", kind="bond_type", direction }
        public init(a: Int,b: Int,kind: String,direction: Int? = nil) { self.a=a; self.b=b; self.kind=kind; self.direction=direction }
    }
    public let atoms: [Atom]; public let bonds: [Bond]
    public let schema:String?
    public let unrepresented_stereo:Bool?
    public init(atoms: [Atom], bonds: [Bond],schema:String? = nil,unrepresented_stereo:Bool? = nil) {
        self.atoms=atoms; self.bonds=bonds;self.schema=schema;self.unrepresented_stereo=unrepresented_stereo
    }
}

public struct DrawingPoint: Codable, Sendable, Equatable {
    public var x: Double; public var y: Double
    public init(x: Double,y: Double) { self.x=x; self.y=y }
    public var length: Double { hypot(x,y) }
    public func distance(_ p: Self) -> Double { (self-p).length }
    public static func - (a: Self,b: Self) -> Self { .init(x: a.x-b.x,y: a.y-b.y) }
    public static func + (a: Self,b: Self) -> Self { .init(x: a.x+b.x,y: a.y+b.y) }
    public static func * (a: Self,b: Double) -> Self { .init(x: a.x*b,y: a.y*b) }
}

public struct ChemicalDepiction: Codable, Sendable {
    public struct Stroke: Codable, Sendable { public let a: DrawingPoint; public let b: DrawingPoint; public let dashed: Bool }
    public struct Label: Codable, Sendable { public let point: DrawingPoint; public let text: String; public let element: String }
    public struct Wedge: Codable, Sendable {
        public let center: Int; public let neighbor: Int
        public let a: DrawingPoint; public let b: DrawingPoint; public let hashed: Bool
    }
    public struct Alkene: Codable, Sendable { public let a:Int;public let b:Int }
    public let points: [DrawingPoint]
    public let strokes: [Stroke]
    public let labels: [Label]
    public let wedges: [Wedge]
    public let alkenes: [Alkene]
    public let stereochemistryOmitted: Bool
    /// Bounded, deterministic template layout. A failed quality gate returns
    /// unavailable rather than a distorted or crossing chemical structure.
    public static func make(_ graph: MolecularDrawingGraph) -> Self? {
        let n=graph.atoms.count
        guard n>0,n<=256 else { return nil }
        if let first=attempt(graph,root: 0) { return first }
        if let bridge=attempt(graph,root: 0,theta: true) { return bridge }
        var degrees=[Int](repeating: 0,count: n)
        for bond in graph.bonds where (0..<n).contains(bond.a) && (0..<n).contains(bond.b) { degrees[bond.a] += 1;degrees[bond.b] += 1 }
        let hubs=graph.atoms.indices.sorted { degrees[$0] == degrees[$1] ? $0 < $1 : degrees[$0] > degrees[$1] }
        var roots=Set([0])
        // A fixed work budget. Changing the seed can avoid a greedy branch
        // dead end without weakening the final geometry/overlap/crossing gate.
        for root in ([n-1,n/2]+Array(hubs.prefix(8))).prefix(10) where roots.insert(root).inserted {
            if let value=attempt(graph,root: root) { return value }
        }
        return nil
    }

    private static func attempt(_ graph: MolecularDrawingGraph,root: Int,theta: Bool = false) -> Self? {
        let n = graph.atoms.count
        guard n > 0, n <= 256, graph.bonds.count <= 512,
              graph.atoms.enumerated().allSatisfy({ $0.offset == $0.element.index }) else { return nil }
        var adjacency = [[Int]](repeating: [],count: n)
        var seen = Set<String>()
        for bond in graph.bonds {
            guard (0..<n).contains(bond.a), (0..<n).contains(bond.b), bond.a != bond.b,
                  ["single","double","triple","aromatic","dative"].contains(bond.kind),
                  seen.insert("\(min(bond.a,bond.b)),\(max(bond.a,bond.b))").inserted else { return nil }
            adjacency[bond.a].append(bond.b); adjacency[bond.b].append(bond.a)
        }
        adjacency = adjacency.map { $0.sorted() }
        // Shortest cycle through each edge; unique edge sets preserve ring order.
        var rings = [[Int]](), ringKeys = Set<String>()
        for bond in graph.bonds {
            var parent = [Int](repeating: -1,count: n), queue = [bond.a], head = 0
            parent[bond.a] = bond.a
            while head < queue.count && parent[bond.b] < 0 {
                let u = queue[head]; head += 1
                for v in adjacency[u] where parent[v] < 0 {
                    if (u == bond.a && v == bond.b) || (u == bond.b && v == bond.a) { continue }
                    parent[v] = u; queue.append(v)
                }
            }
            guard parent[bond.b] >= 0 else { continue }
            var ring = [bond.b], u = bond.b
            while u != bond.a { u = parent[u]; ring.append(u) }
            let key = ring.indices.map { i in
                let a=ring[i], b=ring[(i+1)%ring.count]; return "\(min(a,b)),\(max(a,b))"
            }.sorted().joined(separator: ";")
            if ringKeys.insert(key).inserted { rings.append(ring) }
        }
        rings.sort { $0.count == $1.count ? $0.lexicographicallyPrecedes($1) : $0.count < $1.count }
        var positions = [DrawingPoint?](repeating: nil,count: n)
        var finishedRings = Set<Int>()
        if theta {
            // A theta block has three internally disjoint paths between its
            // two degree-three bridgeheads. Its short cycles need not both be
            // regular polygons. Seed the outer perimeter and the inner bridge.
            let heads=adjacency.indices.filter { adjacency[$0].count==3 }
            guard heads.count==2,adjacency.indices.allSatisfy({ heads.contains($0) || adjacency[$0].count==2 }) else { return nil }
            let a=heads[0],b=heads[1]
            var paths=[[Int]]()
            for neighbor in adjacency[a] {
                var path=[a,neighbor],previous=a,current=neighbor
                while current != b && path.count<=n {
                    guard let next=adjacency[current].first(where:{$0 != previous}) else { return nil }
                    previous=current;current=next;path.append(current)
                }
                guard current==b else { return nil };paths.append(path)
            }
            paths.sort { $0.count == $1.count ? $0.lexicographicallyPrecedes($1) : $0.count>$1.count }
            let outer=paths[0]+paths[1].reversed().dropFirst().dropLast()
            guard Set(outer).count==outer.count else { return nil }
            let radius=1/(2*sin(.pi/Double(outer.count)))
            for (i,atom) in outer.enumerated() { let angle=Double(i)*2 * .pi/Double(outer.count);positions[atom] = .init(x:radius*cos(angle),y:radius*sin(angle)) }
            let d=positions[b]!-positions[a]!,normal=DrawingPoint(x:-d.y,y:d.x)*(1/max(d.length,0.001))
            for i in 1..<(paths[2].count-1) {
                let t=Double(i)/Double(paths[2].count-1)
                positions[paths[2][i]]=positions[a]!+d*t+normal*(0.35*sin(t*2 * .pi))
            }
        }
        if rings.isEmpty && n>2 && graph.bonds.count==n-1 {
            func longestPath(from start: Int) -> [Int] {
                var parent=[Int](repeating: -1,count:n),queue=[start],head=0
                parent[start]=start
                while head<queue.count {
                    let i=queue[head];head += 1
                    for j in adjacency[i] where parent[j]<0 { parent[j]=i;queue.append(j) }
                }
                guard queue.count==n else { return [] }
                var path=[queue.last!]
                while path.last! != start { path.append(parent[path.last!]) }
                return path
            }
            if let end=longestPath(from:root).first {
                let backbone=longestPath(from:end)
                positions[backbone[0]] = .init(x:0,y:0)
                var angle=0.0
                for i in 1..<backbone.count {
                    if i>1 {
                        let central=backbone[i-1]
                        let linear=graph.bonds.contains { ($0.a==central || $0.b==central) && $0.kind=="triple" }
                        if !linear { angle=abs(angle)<0.1 ? .pi/3 : 0 }
                    }
                    positions[backbone[i]] = positions[backbone[i-1]]! + DrawingPoint(x:cos(angle),y:sin(angle))
                }
            }
        }
        func cross(_ a: DrawingPoint,_ b: DrawingPoint) -> Double { a.x*b.y-a.y*b.x }
        func clearance(_ candidate: [Int: DrawingPoint]) -> Double {
            var score = 0.0
            for i in candidate.keys.sorted() where positions[i] == nil {
                let p = candidate[i]!
                for j in 0..<n where positions[j] != nil && !candidate.keys.contains(j) {
                    let distance = p.distance(positions[j]!)
                    score += 1 / max(0.02,distance*distance*distance*distance)
                }
                score += 0.001*abs(p.y)
            }
            // Choose an uncrossed branch now instead of discovering a crossing
            // only after the complete tree has been placed.
            for edge in graph.bonds where positions[edge.a] == nil || positions[edge.b] == nil {
                guard let p=candidate[edge.a] ?? positions[edge.a],let q=candidate[edge.b] ?? positions[edge.b] else { continue }
                for old in graph.bonds {
                    if edge.a==old.a || edge.a==old.b || edge.b==old.a || edge.b==old.b { continue }
                    guard let s=positions[old.a],let t=positions[old.b] else { continue }
                    if cross(q-p,s-p)*cross(q-p,t-p) < -1e-8 && cross(t-s,p-s)*cross(t-s,q-s) < -1e-8 { score += 10_000 }
                }
            }
            return score
        }
        func ringPlacement(_ ring: [Int], seed: DrawingPoint) -> [Int: DrawingPoint]? {
            let count = ring.count, radius = 1/(2*sin(.pi/Double(count)))
            let base = (0..<count).map { DrawingPoint(x: radius*cos(Double($0)*2 * .pi/Double(count)),y: radius*sin(Double($0)*2 * .pi/Double(count))) }
            let known = ring.indices.filter { positions[ring[$0]] != nil }
            var best: [Int: DrawingPoint]?, bestScore = Double.infinity
            for sign in [-1.0,1.0] {
                let mirror = base.map { DrawingPoint(x: $0.x,y: $0.y*sign) }
                let anchor = known.first ?? 0
                let angles: [Double]
                if known.count >= 2 {
                    let second = known[1], target = positions[ring[second]]! - positions[ring[anchor]]!, original = mirror[second]-mirror[anchor]
                    guard abs(target.length-original.length) < 0.12 else { continue }
                    angles = [atan2(target.y,target.x)-atan2(original.y,original.x)]
                } else { angles = (0..<12).map { Double($0)*Double.pi/6 } }
                for angle in angles {
                    let origin = positions[ring[anchor]] ?? seed
                    var candidate = [Int: DrawingPoint]()
                    for i in ring.indices {
                        let p = mirror[i]-mirror[anchor]
                        candidate[ring[i]] = origin + DrawingPoint(x: p.x*cos(angle)-p.y*sin(angle),y: p.x*sin(angle)+p.y*cos(angle))
                    }
                    guard known.allSatisfy({ candidate[ring[$0]]!.distance(positions[ring[$0]]!) < 0.12 }) else { continue }
                    let score = clearance(candidate)
                    if score < bestScore { best = candidate; bestScore = score }
                }
            }
            return best
        }
        while positions.contains(where: { $0 == nil }) {
            var progressed = false
            for (index,ring) in rings.enumerated() where !finishedRings.contains(index) {
                if ring.allSatisfy({ positions[$0] != nil }) { finishedRings.insert(index); continue }
                guard ring.contains(where: { positions[$0] != nil }) else { continue }
                if let candidate = ringPlacement(ring, seed: .init(x: 0,y: 0)) {
                    for (i,p) in candidate where positions[i] == nil { positions[i] = p }
                    finishedRings.insert(index); progressed = true; break
                }
            }
            if progressed { continue }
            if let edge = graph.bonds.first(where: { (positions[$0.a] == nil) != (positions[$0.b] == nil) }) {
                let a = positions[edge.a] == nil ? edge.b : edge.a, b = a == edge.a ? edge.b : edge.a
                let neighbors = adjacency[a].filter { positions[$0] != nil }
                let linear = graph.bonds.contains { ($0.a == a || $0.b == a) && $0.kind == "triple" }
                var best = DrawingPoint(x: 0,y: 0), bestScore = Double.infinity
                for step in 0..<12 {
                    let angle = Double(step)*Double.pi/6
                    let direction = DrawingPoint(x: cos(angle),y: sin(angle)), p = positions[a]! + direction
                    var score = clearance([b:p])
                    for neighbor in neighbors {
                        let v = positions[neighbor]! - positions[a]!
                        let cosine = (v.x*direction.x+v.y*direction.y)/v.length
                        let target = linear ? -1.0 : -0.5
                        score += 40*pow(cosine-target,2)
                    }
                    // Compact zigzag chains; deterministic tie break.
                    score += 0.01*abs(p.y) + Double(step)*0.000001
                    if score < bestScore { bestScore=score; best=p }
                }
                positions[b] = best
            } else {
                let offset = (positions.compactMap { $0?.x }.max().map { $0+3 }) ?? 0
                let availableRings=rings.indices.filter { !finishedRings.contains($0) && rings[$0].allSatisfy { positions[$0] == nil } }
                let nextRing=positions.allSatisfy({ $0 == nil }) ? (availableRings.first { rings[$0].contains(root) } ?? availableRings.first) : availableRings.first
                if let index = nextRing,
                   let candidate = ringPlacement(rings[index],seed: .init(x: offset,y: 0)) {
                    for (i,p) in candidate { positions[i]=p }; finishedRings.insert(index)
                } else if let index = positions.allSatisfy({ $0 == nil }) ? root : positions.firstIndex(where: { $0 == nil }) { positions[index] = .init(x: offset,y: 0) }
            }
        }
        var points = positions.compactMap { $0 }
        func valid(_ points: [DrawingPoint]) -> Bool {
            guard points.allSatisfy({ $0.x.isFinite && $0.y.isFinite }),
                  graph.bonds.allSatisfy({ (0.85...1.15).contains(points[$0.a].distance(points[$0.b])) }) else { return false }
            for i in 0..<n { for j in 0..<i where points[i].distance(points[j]) < 0.42 { return false } }
            for i in graph.bonds.indices {
                let a=graph.bonds[i]
                for j in 0..<i {
                    let b=graph.bonds[j]
                    if a.a == b.a || a.a == b.b || a.b == b.a || a.b == b.b { continue }
                    let p=points[a.a],q=points[a.b],s=points[b.a],t=points[b.b]
                    if cross(q-p,s-p)*cross(q-p,t-p) < -1e-8 && cross(t-s,p-s)*cross(t-s,q-s) < -1e-8 { return false }
                }
            }
            return true
        }
        if !valid(points) && !graph.bonds.contains(where: { $0.kind=="triple" }) {
            // Local collision repair for an otherwise rejected template.
            // Linear triple bonds are never deformed. Keep the same quality gate.
            let original=points
            for _ in 0..<180 {
                var force=[DrawingPoint](repeating: .init(x:0,y:0),count:n)
                for bond in graph.bonds {
                    let d=points[bond.b]-points[bond.a],length=max(d.length,0.001),f=d*((length-1)/length)
                    force[bond.a]=force[bond.a]+f;force[bond.b]=force[bond.b]-f
                }
                for i in 0..<n {
                    for j in 0..<i where !adjacency[i].contains(j) {
                        var d=points[i]-points[j],length=d.length
                        if length<0.001 { let a=Double(i*17+j*31);d = .init(x:cos(a)*0.001,y:sin(a)*0.001);length=0.001 }
                        if length<0.7 {
                            let f=d*((0.7-length)/length)
                            force[i]=force[i]+f;force[j]=force[j]-f
                        }
                    }
                    force[i]=force[i]+(original[i]-points[i])*0.003
                }
                for i in 0..<n {
                    let f=force[i]*0.15,step=f.length>0.08 ? f*(0.08/f.length) : f
                    points[i]=points[i]+step
                }
            }
        }
        guard valid(points) else { return nil }
        var alkenes=[Alkene](),handledDirections=Set<Int>()
        func directions(at center:Int) -> [(index:Int,neighbor:Int,sign:Int)] {
            graph.bonds.indices.compactMap { i in
                let b=graph.bonds[i]
                guard b.kind=="single",let d=b.direction,d != 0 else { return nil }
                if b.a==center { return (i,b.b,d) }
                if b.b==center { return (i,b.a,-d) }
                return nil
            }
        }
        func product(_ p:[DrawingPoint],a:Int,b:Int,left:Int,right:Int) -> Double {
            let axis=p[b]-p[a]
            return cross(axis,p[left]-p[a])*cross(axis,p[right]-p[b])
        }
        func component(from start:Int,excluding other:Int) -> [Int]? {
            var queue=[start],seen=Set([start]),head=0
            while head<queue.count {
                let i=queue[head];head += 1
                for j in adjacency[i] {
                    if (i==start && j==other)||(i==other && j==start) { continue }
                    if seen.insert(j).inserted { queue.append(j) }
                }
            }
            return seen.contains(other) ? nil : queue
        }
        func fanValid(_ p:[DrawingPoint],center:Int,other:Int) -> Bool {
            let others=adjacency[center].filter{$0 != other}
            guard others.count==2 else { return others.count<2 }
            let axis=p[other]-p[center]
            return cross(axis,p[others[0]]-p[center])*cross(axis,p[others[1]]-p[center]) < -1e-8
        }
        func annotationsMatch(_ p:[DrawingPoint],a:Int,b:Int) -> Bool {
            guard fanValid(p,center:a,other:b),fanValid(p,center:b,other:a) else { return false }
            return directions(at:a).allSatisfy { l in directions(at:b).allSatisfy { r in
                let value=product(p,a:a,b:b,left:l.neighbor,right:r.neighbor)
                return abs(value)>1e-8 && (value>0)==(l.sign*r.sign>0)
            } }
        }
        for bond in graph.bonds where bond.kind=="double" {
            let left=directions(at:bond.a),right=directions(at:bond.b)
            guard let l=left.first,let r=right.first else { continue }
            // At a trigonal endpoint the two substituents must straddle the
            // double-bond axis. A deformed ring template can otherwise satisfy
            // one slash while depicting the opposite CIP-selected geometry.
            for (center,other) in [(bond.a,bond.b),(bond.b,bond.a)] where !fanValid(points,center:center,other:other) {
                let neighbors=adjacency[center].filter{$0 != other}
                guard neighbors.count==2,let side=component(from:center,excluding:other) else { continue }
                let u=points[neighbors[0]]-points[center],v=points[neighbors[1]]-points[center]
                let sum=u*(1/u.length)+v*(1/v.length),target=points[center]-points[other]
                guard sum.length>0.01 else { continue }
                let angle=atan2(target.y,target.x)-atan2(sum.y,sum.x),origin=points[center]
                var candidate=points
                for i in side {
                    let d=points[i]-origin
                    candidate[i]=origin+DrawingPoint(x:d.x*cos(angle)-d.y*sin(angle),y:d.x*sin(angle)+d.y*cos(angle))
                }
                if valid(candidate) && fanValid(candidate,center:center,other:other) && alkenes.allSatisfy({ annotationsMatch(candidate,a:$0.a,b:$0.b) }) { points=candidate }
            }
            let wantPositive=l.sign*r.sign>0
            let initial=product(points,a:bond.a,b:bond.b,left:l.neighbor,right:r.neighbor)
            if abs(initial)<1e-8 || (initial>0) != wantPositive {
                let sides=[component(from:bond.a,excluding:bond.b),component(from:bond.b,excluding:bond.a)].compactMap{$0}.sorted{$0.count<$1.count}
                let origin=points[bond.a],axis=points[bond.b]-origin
                let normal=DrawingPoint(x:-axis.y,y:axis.x)*(1/axis.length)
                for side in sides {
                    var candidate=points
                    for i in side {
                        let d=points[i]-origin,distance=d.x*normal.x+d.y*normal.y
                        candidate[i]=points[i]-normal*(2*distance)
                    }
                    let p=product(candidate,a:bond.a,b:bond.b,left:l.neighbor,right:r.neighbor)
                    if abs(p)>1e-8 && (p>0)==wantPositive && valid(candidate) { points=candidate;break }
                }
            }
            // Verify every annotated substituent, not just the first pair.
            let matches=annotationsMatch(points,a:bond.a,b:bond.b)
            if matches {
                alkenes.append(.init(a:bond.a,b:bond.b))
                handledDirections.formUnion(left.map(\.index));handledDirections.formUnion(right.map(\.index))
            }
        }
        let labels: [Label] = graph.atoms.compactMap { atom in
            let isolated = adjacency[atom.index].isEmpty
            let neighbors=adjacency[atom.index]
            var straight=false
            if neighbors.count==2 {
                let a=points[neighbors[0]]-points[atom.index],b=points[neighbors[1]]-points[atom.index]
                straight=(a.x*b.x+a.y*b.y)/(a.length*b.length) < -0.99
            }
            // A degree-two carbon on a straight segment has no visible vertex.
            // Write it explicitly so a bridge/sp chain never loses an atom.
            guard atom.symbol != "C" || isolated || straight || atom.formal_charge != 0 || atom.isotope != nil else { return nil }
            let hydrogen = (atom.bracket_h_count ?? 0)+(atom.implicit_h_count ?? 0)
            let h = atom.symbol == "H" || hydrogen == 0 ? "" : "H"+(hydrogen == 1 ? "" : String(hydrogen))
            let charge = atom.formal_charge == 0 ? "" : (abs(atom.formal_charge) == 1 ? "" : String(abs(atom.formal_charge)))+(atom.formal_charge > 0 ? "+" : "−")
            return Label(point: points[atom.index],text: (atom.isotope.map(String.init) ?? "")+atom.symbol+h+charge,element: atom.symbol)
        }
        let labelled = Set(graph.atoms.filter { a in labels.contains(where: { $0.point == points[a.index] }) }.map(\.index))
        var wedges=[Wedge](), resolvedStereo=Set<Int>(), wedgeEdges=Set<String>()
        for atom in graph.atoms where atom.chirality_token != nil {
            guard let token=atom.chirality_token,["@","@@","@TH1","@TH2"].contains(token),
                  let order=atom.stereo_neighbor_order,order.count==4,Set(order).count==4,
                  order.allSatisfy({ $0 == -1 || adjacency[atom.index].contains($0) }),
                  order.filter({ $0 >= 0 }).count==adjacency[atom.index].count else { continue }
            let candidates=adjacency[atom.index].filter { neighbor in
                let key="\(min(atom.index,neighbor)),\(max(atom.index,neighbor))"
                return !wedgeEdges.contains(key) && graph.bonds.contains { b in
                    b.kind=="single" && ((b.a==atom.index && b.b==neighbor)||(b.b==atom.index && b.a==neighbor))
                }
            }.sorted {
                let a=graph.atoms[$0].chirality_token == nil,b=graph.atoms[$1].chirality_token == nil
                return a == b ? $0 < $1 : a
            }
            for neighbor in candidates {
                // Neighbor order is the parser's SMILES order, including the
                // implicit-H slot. Lift one neighbor out of the drawing plane;
                // a signed tetrahedral volume determines wedge versus hash.
                let vectors=order.map { i -> (Double,Double,Double) in
                    let p=i == -1 ? points[atom.index] : points[i]
                    return (p.x,p.y,i == neighbor ? 1 : 0)
                }
                let u=vectors[0],v=vectors[1],w=vectors[2],t=vectors[3]
                let ax=u.0-t.0,ay=u.1-t.1,az=u.2-t.2
                let bx=v.0-t.0,by=v.1-t.1,bz=v.2-t.2
                let cx=w.0-t.0,cy=w.1-t.1,cz=w.2-t.2
                let volume=ax*(by*cz-bz*cy)-ay*(bx*cz-bz*cx)+az*(bx*cy-by*cx)
                guard abs(volume)>0.1 else { continue }
                let positive=["@","@TH1"].contains(token)
                let hashed=(volume>0) != positive
                wedges.append(.init(center: atom.index,neighbor: neighbor,a: points[atom.index],b: points[neighbor],hashed: hashed))
                wedgeEdges.insert("\(min(atom.index,neighbor)),\(max(atom.index,neighbor))")
                resolvedStereo.insert(atom.index);break
            }
        }
        var strokes = [Stroke]()
        for bond in graph.bonds {
            if wedgeEdges.contains("\(min(bond.a,bond.b)),\(max(bond.a,bond.b))") { continue }
            let a=points[bond.a],b=points[bond.b],delta=b-a, direction=delta*(1/delta.length)
            let normal=DrawingPoint(x: -direction.y,y: direction.x)
            let ring = rings.first { r in r.indices.contains { i in
                let x=r[i],y=r[(i+1)%r.count]; return (x == bond.a && y == bond.b)||(x == bond.b && y == bond.a)
            } }
            var inward = 1.0
            if let ring {
                let center = ring.map { points[$0] }.reduce(DrawingPoint(x: 0,y: 0),+) * (1/Double(ring.count))
                let d=center-(a+b)*0.5; inward = d.x*normal.x+d.y*normal.y < 0 ? -1 : 1
            }
            let lines: [(Double,Bool)]
            switch bond.kind {
            case "triple": lines=[(-0.13,false),(0,false),(0.13,false)]
            case "double": lines=ring == nil ? [(-0.075,false),(0.075,false)] : [(0,false),(0.16*inward,false)]
            case "aromatic": lines=[(0,false),(0.15*inward,true)]
            default: lines=[(0,false)]
            }
            for (offset,dashed) in lines {
                let trimA=labelled.contains(bond.a) ? 0.19 : (ring != nil && offset != 0 ? 0.16 : 0)
                let trimB=labelled.contains(bond.b) ? 0.19 : (ring != nil && offset != 0 ? 0.16 : 0)
                strokes.append(.init(a: a+direction*trimA+normal*offset,b: b-direction*trimB+normal*offset,dashed: dashed))
            }
            if bond.kind=="dative" {
                let tip=b-direction*(labelled.contains(bond.b) ? 0.19 : 0)
                let base=tip-direction*0.15
                strokes.append(.init(a:base+normal*0.075,b:tip,dashed:false))
                strokes.append(.init(a:base-normal*0.075,b:tip,dashed:false))
            }
        }
        return Self(points: points,strokes: strokes,labels: labels,wedges: wedges,alkenes: alkenes,
                    stereochemistryOmitted: graph.unrepresented_stereo==true || graph.atoms.contains { $0.chirality_token != nil && !resolvedStereo.contains($0.index) } || graph.bonds.indices.contains { (graph.bonds[$0].direction ?? 0) != 0 && !handledDirections.contains($0) })
    }

    /// Same primitives as the native canvas, useful for reproducible inspection.
    /// Molecular coordinates use upward Y. Screen/SVG Y points downward, so
    /// rendering MUST reverse Y while retaining the wedge/hash direction.
    public func screenPoint(_ p: DrawingPoint) -> DrawingPoint { .init(x:p.x,y:-p.y) }

    public func svg() -> String {
        let screen=points.map(screenPoint)
        let minX=(screen.map(\.x).min() ?? 0)-0.65, minY=(screen.map(\.y).min() ?? 0)-0.65
        let width=(screen.map(\.x).max() ?? 0)-minX+0.65, height=(screen.map(\.y).max() ?? 0)-minY+0.65
        func number(_ x: Double) -> String { String(format: "%.5f",locale: Locale(identifier: "en_US_POSIX"),x) }
        func y(_ value:Double) -> String { number(screenPoint(.init(x:0,y:value)).y) }
        func escape(_ s: String) -> String { s.replacingOccurrences(of: "&",with: "&amp;").replacingOccurrences(of: "<",with: "&lt;").replacingOccurrences(of: "\"",with: "&quot;") }
        var s="<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"\(number(minX)) \(number(minY)) \(number(width)) \(number(height))\"><rect x=\"\(number(minX))\" y=\"\(number(minY))\" width=\"\(number(width))\" height=\"\(number(height))\" fill=\"white\"/>"
        for stroke in strokes { s += "<path d=\"M\(number(stroke.a.x)),\(y(stroke.a.y)) L\(number(stroke.b.x)),\(y(stroke.b.y))\" stroke=\"#162527\" stroke-width=\".035\"\(stroke.dashed ? " stroke-dasharray=\".09 .06\"" : "")/>" }
        for wedge in wedges {
            let d=wedge.b-wedge.a,normal=DrawingPoint(x: -d.y,y: d.x)*(0.15/d.length)
            if wedge.hashed {
                for i in 1...7 {
                    let t=Double(i)/8,p=wedge.a+d*t,a=p-normal*t,b=p+normal*t
                    s += "<path d=\"M\(number(a.x)),\(y(a.y)) L\(number(b.x)),\(y(b.y))\" stroke=\"#162527\" stroke-width=\".025\"/>"
                }
            } else {
                let a=wedge.b-normal,b=wedge.b+normal
                s += "<path d=\"M\(number(wedge.a.x)),\(y(wedge.a.y)) L\(number(a.x)),\(y(a.y)) L\(number(b.x)),\(y(b.y)) Z\" fill=\"#162527\"/>"
            }
        }
        for label in labels {
            s += "<text x=\"\(number(label.point.x))\" y=\"\(y(label.point.y))\" text-anchor=\"middle\" dominant-baseline=\"central\" font-family=\"sans-serif\" font-size=\".28\" fill=\"#162527\" stroke=\"white\" stroke-width=\".1\" paint-order=\"stroke\">\(escape(label.text))</text>"
        }
        return s+"</svg>"
    }
}
