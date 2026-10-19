"""Pinned non-stereo, single-bond canonical partition experiment, not production.

Adapts RDKit new_canon.cpp/.h (Copyright2014 Greg Landrum, Roger Sayle
pseudocode), FindRings and hanoiSort. BSD3 notice: crates/smiles/RDKIT-DATA-LICENSE.txt.
Stable key sorting reproduces hanoi's stable equivalence-group merge; partition
activation/refinement and tie breaking retain the pinned algorithm's order.
RDKit supplies graph/properties/fast-ring metadata in this experiment only.
"""
from rdkit import Chem


class Ranking:
    def __init__(self, mol):
        self.n = mol.GetNumAtoms()
        assert 0 < self.n <= 64
        assert all(b.GetBondType() == Chem.BondType.SINGLE for b in mol.GetBonds())
        assert all(a.GetChiralTag() == Chem.ChiralType.CHI_UNSPECIFIED for a in mol.GetAtoms())
        self.adj = [[a.GetIdx() for a in atom.GetNeighbors()] for atom in mol.GetAtoms()]
        self.base = [(a.GetAtomMapNum(),a.GetDegree(),a.GetAtomicNum(),a.GetIsotope(),
                      a.GetTotalNumHs(),a.GetFormalCharge() & 0xffffffff)
                     for a in mol.GetAtoms()]
        Chem.FastFindRings(mol)
        rings = mol.GetRingInfo().AtomRings()
        self.ring_counts = [sum(i in ring for ring in rings) for i in range(self.n)]
        self.rank = [0]*self.n
        self.order = list(range(self.n))
        self.count = [self.n]+[0]*(self.n-1)
        self.next = [-2]*self.n
        self.active = -1
        self.profiles = [((),())]*self.n

    def key(self, atom, special):
        neighbors = tuple(sorted((self.rank[j] for j in self.adj[atom]),reverse=True))
        if special:
            return (*self.profiles[atom], neighbors)
        return (self.rank[atom],self.base[atom],neighbors)

    def activate(self):
        self.active = -1
        self.next = [-2]*self.n
        position = 0
        while position < self.n:
            head = self.order[position]
            assert self.count[head] > 0
            if self.count[head] > 1:
                self.next[head] = self.active;self.active = head
                position += self.count[head]
            else:
                position += 1

    def touch(self, vertices):
        partitions = set(self.rank[j] for v in vertices for j in self.adj[v])
        for position in sorted(partitions):
            head = self.order[position]
            if self.count[head] > 1 and self.next[head] == -2:
                self.next[head] = self.active;self.active = head

    def refine(self, special=False):
        iterations = 0
        while self.active != -1:
            iterations += 1;assert iterations < 100_000
            head = self.active;self.active = self.next[head];self.next[head] = -2
            length, offset = self.count[head],self.rank[head]
            members = self.order[offset:offset+length]
            keys = {v:self.key(v,special) for v in members}
            members.sort(key=keys.__getitem__)
            self.order[offset:offset+length] = members
            for v in members:self.count[v] = 0
            group = 0
            for i in range(1,length+1):
                if i == length or keys[members[i]] != keys[members[group]]:
                    self.count[members[group]] = i-group;group=i
            first = self.count[members[0]]
            current = offset
            for i in range(first,length):
                v = members[i]
                if self.count[v]:current = offset+i
                self.rank[v] = current
            self.touch(members[first:])

    def ring_profiles(self):
        for seed in range(self.n):
            if not self.ring_counts[seed]:continue
            visited = [False]*self.n;last = [False]*self.n
            level = [seed];numbers=[];revisited=[]
            while level:
                current = [False]*self.n;following=[];counts=[0]*self.n
                for v in level:
                    if not self.ring_counts[v]:continue
                    last[v] = True;visited[v] = True
                    for j in self.adj[v]:
                        if not visited[j]:
                            visited[j] = True;current[j] = True;following.append(j)
                for v in following:
                    for j in self.adj[v]:
                        if current[j] or last[j]:counts[j] += 1
                revisited.extend(sorted(c for c in counts if c)+[-1])
                numbers.extend([len(following),-1]);last = current;level = following
            self.profiles[seed] = (tuple(numbers),tuple(revisited))

    def run(self):
        self.activate();self.refine()
        # The special chirality phase is inert for this explicitly non-stereo subset.
        sym = sum(self.count[v] for v in self.order if self.ring_counts[v] and self.count[v]>2)
        ring_atoms = sum(c>0 for c in self.ring_counts)
        branching = any(self.ring_counts[v]>1 and self.count[v]>1 for v in self.order)
        if any(c==0 for c in self.count) and ring_atoms and 2*sym>ring_atoms and branching:
            self.ring_profiles();self.activate();self.refine(special=True)
        position = 0
        while position < self.n:
            head = self.order[position];old = self.rank[head]
            while self.count[head] > 1:
                end = self.rank[head] + self.count[head]-1
                vertex = self.order[end];self.rank[vertex] = end
                self.count[head] -= 1;self.count[vertex] = 1
                if not self.adj[vertex]:continue
                self.touch([vertex]);self.refine()
            if self.rank[head] == old:position += 1
        return self.rank


if __name__ == '__main__':
    import hashlib,json,random
    from pathlib import Path
    data = json.loads(Path('tasks/evidence/catalog289-sps/remaining-metal5-diagnostic.json').read_text())
    params = Chem.SmilesParserParams();params.sanitize=False;params.removeHs=False
    controls = []
    rng = random.Random(20261018)
    for row in data:
        if row['index'] in [18084,53180]:continue  # valid NT geometry is a separate rank-prefix proof
        mol = Chem.MolFromSmiles(row['smiles'],params);mol.UpdatePropertyCache(False)
        for trial in range(65):
            order = list(range(mol.GetNumAtoms()))
            if trial:rng.shuffle(order)
            current = Chem.RenumberAtoms(mol,order)
            actual = Ranking(current).run();expected=list(Chem.CanonicalRankAtoms(current))
            assert actual == expected,(row['index'],trial,actual,expected)
        controls.append(dict(index=row['index'],permutations=65,exact_pinned_ranks=True))
    result = dict(controls=controls,production_implementation=False,
                  source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  scope='Non-stereo single-bond graph partition subset only; RDKit graph/fast-ring metadata in prototype')
    Path('tasks/evidence/catalog289-sps/cluster-canonical-prototype.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))
