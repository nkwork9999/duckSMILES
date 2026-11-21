"""Batch descriptor reads on a host-provisioned duckSMILES connection."""
from collections.abc import Mapping
from dataclasses import dataclass
from enum import Enum
from sys import float_info
from time import sleep

from duckdb import ConstraintException, DuckDBPyConnection, TransactionException


@dataclass(frozen=True)
class DescriptorRow:
    formula: str | None
    weight: float | None


def _validate_molecules(molecules):
    if not isinstance(molecules, list) or any(value is not None and type(value) is not str for value in molecules):
        raise TypeError("molecules must be a list of SMILES strings or None")


class MoleculeReader:
    def __init__(self, connection: DuckDBPyConnection):
        self._db = connection

    def describe(self, molecules: list[str | None]) -> list[DescriptorRow]:
        """Compute each descriptor per distinct SMILES, then restore input order."""
        _validate_molecules(molecules)
        rows = self._db.execute("""
            WITH inputs AS (
                SELECT smiles, ord FROM unnest(?::VARCHAR[]) WITH ORDINALITY AS t(smiles, ord)
            ), properties AS MATERIALIZED (
                SELECT smiles, mol_formula(smiles) AS formula, mol_weight(smiles) AS weight
                FROM (SELECT DISTINCT smiles FROM inputs WHERE smiles IS NOT NULL)
            )
            SELECT formula, weight FROM inputs LEFT JOIN properties USING (smiles) ORDER BY ord
        """, [molecules]).fetchall()
        return [DescriptorRow(*row) for row in rows]


@dataclass(frozen=True)
class DescriptorRevision:
    implementation: str
    standardization: str

    def __post_init__(self):
        for value in (self.implementation, self.standardization):
            if type(value) is not str or not value or len(value) > 255:
                raise ValueError("describe each implementation and standardization policy with a version ID")


@dataclass(frozen=True)
class CachedDescriptors:
    rows: tuple[DescriptorRow, ...]
    new_entries: int


class DescriptorCache:
    """Keep native descriptor results, including invalid inputs, per host-provided revision."""

    def __init__(self, connection: DuckDBPyConnection):
        self._db = connection
        with connection.cursor() as tx:
            tx.execute("BEGIN TRANSACTION")
            tx.execute("CREATE SCHEMA IF NOT EXISTS ducksmiles_cache")
            tx.execute("""
                CREATE TABLE IF NOT EXISTS ducksmiles_cache.descriptors (
                    implementation VARCHAR, standardization VARCHAR, digest UUID, smiles VARCHAR,
                    formula VARCHAR, weight DOUBLE,
                    PRIMARY KEY (implementation, standardization, digest, smiles)
                )
            """)
            tx.execute("COMMIT")

    def describe(self, revision: DescriptorRevision, molecules: list[str | None]) -> CachedDescriptors:
        if type(revision) is not DescriptorRevision:
            raise TypeError("expected a DescriptorRevision")
        _validate_molecules(molecules)
        unique = list(dict.fromkeys(value for value in molecules if value is not None))
        if not unique:
            return CachedDescriptors(tuple(DescriptorRow(None, None) for _ in molecules), 0)
        scope = [revision.implementation, revision.standardization]
        read = """
            WITH cached AS (
                SELECT digest, smiles, formula, weight FROM ducksmiles_cache.descriptors
                WHERE implementation = ? AND standardization = ?
            ), inputs AS (
                SELECT smiles, md5(smiles)::UUID AS digest FROM unnest(?::VARCHAR[]) AS t(smiles)
            )
            SELECT smiles, formula, weight FROM inputs JOIN cached USING (digest, smiles)
        """
        write = """
            INSERT INTO ducksmiles_cache.descriptors
            WITH input AS (
                SELECT smiles, md5(smiles)::UUID AS digest FROM unnest(?::VARCHAR[]) AS t(smiles)
            )
            SELECT ?, ?, digest, smiles, mol_formula(smiles), mol_weight(smiles)
            FROM input ON CONFLICT DO NOTHING RETURNING smiles, formula, weight
        """
        descriptors = {}
        missing = unique
        new_entries = 0
        with self._db.cursor() as connection:
            for attempt in range(8):
                found = connection.execute(read, [*scope, missing]).fetchall()
                descriptors.update({smiles: DescriptorRow(formula, weight) for smiles, formula, weight in found})
                missing = [value for value in missing if value not in descriptors]
                if not missing:
                    break
                try:
                    added = connection.execute(write, [missing, *scope]).fetchall()
                except (ConstraintException, TransactionException) as error:
                    if "duplicate key" not in str(error).lower() or attempt == 7:
                        raise
                    sleep(min(0.002 * 2 ** attempt, 0.05))
                    continue
                new_entries += len(added)
                descriptors.update({smiles: DescriptorRow(formula, weight) for smiles, formula, weight in added})
                missing = [value for value in missing if value not in descriptors]
                if not missing:
                    break
            if missing:
                raise RuntimeError("descriptor cache fill did not settle")
        empty = DescriptorRow(None, None)
        return CachedDescriptors(tuple(descriptors.get(value, empty) for value in molecules), new_entries)


@dataclass(frozen=True)
class MorganFingerprint:
    radius: int = 2
    bits: int = 2048

    def __post_init__(self):
        if type(self.radius) is not int or not 0 <= self.radius <= 32:
            raise ValueError("radius must be an integer in [0, 32]")
        if type(self.bits) is not int or not 1 <= self.bits <= 1_048_576:
            raise ValueError("bits must be an integer in [1, 1048576]")


@dataclass(frozen=True)
class MassRange:
    lower: float
    upper: float

    def __post_init__(self):
        if any(type(value) not in (int, float) or not 0 <= value <= float_info.max
               for value in (self.lower, self.upper)) or self.lower > self.upper:
            raise ValueError("mass bounds must be finite, nonnegative and ordered")


@dataclass(frozen=True)
class PatternId:
    value: str

    def __post_init__(self):
        if type(self.value) is not str or not self.value:
            raise ValueError("pattern ID must be a nonempty string")


@dataclass(frozen=True)
class FingerprintFilter:
    reference: str
    threshold: float
    definition: MorganFingerprint = MorganFingerprint()

    def __post_init__(self):
        if type(self.reference) is not str or not self.reference or type(self.definition) is not MorganFingerprint:
            raise ValueError("provide a reference SMILES and MorganFingerprint")
        if type(self.threshold) not in (int, float) or not 0 <= self.threshold <= 1:
            raise ValueError("threshold must be in [0, 1]")


@dataclass(frozen=True)
class MoleculeSearchSpec:
    mass: MassRange | None = None
    fingerprint: FingerprintFilter | None = None
    pattern: PatternId | None = None
    limit: int = 100

    def __post_init__(self):
        for value, expected in ((self.mass, MassRange), (self.fingerprint, FingerprintFilter), (self.pattern, PatternId)):
            if value is not None and type(value) is not expected:
                raise TypeError(f"expected a {expected.__name__}")
        if type(self.limit) is not int or not 1 <= self.limit <= 10_000:
            raise ValueError("limit must be an integer in [1, 10000]")


@dataclass(frozen=True)
class MoleculeSearchHit:
    molecule_id: int
    mass: float | None
    score: float | None


@dataclass(frozen=True)
class MoleculeSearchResult:
    state: str
    hits: tuple[MoleculeSearchHit, ...]


class MoleculeSearch:
    """Search one host-registered relation with a closed set of predicate types."""

    def __init__(self, connection: DuckDBPyConnection, relation: str, patterns: Mapping[PatternId, str]):
        if type(relation) is not str or not relation or "\0" in relation:
            raise ValueError("relation must be a nonempty identifier without NUL")
        self._db = connection
        self._relation = '"' + relation.replace('"', '""') + '"'
        self._patterns = dict(patterns)
        if any(type(key) is not PatternId or type(pattern) is not str or not pattern
               for key, pattern in self._patterns.items()):
            raise TypeError("provide host-defined PatternId to SMARTS mappings")

    def search(self, spec: MoleculeSearchSpec) -> MoleculeSearchResult:
        if type(spec) is not MoleculeSearchSpec:
            raise TypeError("expected a MoleculeSearchSpec")
        mass_sql, structure_sql, score_sql, score_filter = "TRUE", "TRUE", "NULL::DOUBLE", "TRUE"
        parameters = []
        if spec.mass is not None:
            mass_sql = "isfinite(mass) AND mass BETWEEN ?::DOUBLE AND ?::DOUBLE"
            parameters.extend([spec.mass.lower, spec.mass.upper])
        if spec.pattern is not None:
            structure_sql = "mol_has_substructure(smiles, ?::VARCHAR)"
            parameters.append(self._patterns[spec.pattern])
        if spec.fingerprint is not None:
            fp = spec.fingerprint
            definition = [fp.definition.radius, fp.definition.bits]
            reference, = self._db.execute("SELECT morgan_fp_bits(?::VARCHAR, ?, ?)",
                                          [fp.reference, *definition]).fetchone()
            if reference is None:
                return MoleculeSearchResult("invalid_reference", ())
            score_sql = "tanimoto_bit(morgan_fp_bits(smiles, ?, ?), ?::BLOB)"
            score_filter = "score >= ?::DOUBLE"
            parameters.extend([*definition, reference, fp.threshold])
        rows = self._db.execute(f"""
            WITH eligible AS MATERIALIZED (
                SELECT molecule_id, smiles, CASE WHEN isfinite(mass) THEN mass END AS mass
                FROM main.{self._relation} WHERE {mass_sql}
            ), matched AS MATERIALIZED (
                SELECT molecule_id, mass, smiles FROM eligible WHERE {structure_sql}
            ), scored AS (
                SELECT molecule_id, mass, {score_sql} AS score FROM matched
            )
            SELECT molecule_id, mass, score FROM scored WHERE {score_filter}
            ORDER BY score DESC NULLS LAST, molecule_id LIMIT ?
        """, [*parameters, spec.limit]).fetchall()
        return MoleculeSearchResult("ok", tuple(MoleculeSearchHit(*row) for row in rows))


@dataclass(frozen=True)
class SimilarityHit:
    molecule_id: int
    score: float


@dataclass(frozen=True)
class SimilarityResult:
    hits: tuple[SimilarityHit, ...]
    candidates_scored: int


@dataclass(frozen=True)
class SimilarityEdge:
    source: int
    target: int
    score: float


@dataclass(frozen=True)
class SimilarityGraphSummary:
    vertices: int
    comparisons: int
    edges: int


@dataclass(frozen=True)
class ScaffoldSample:
    molecule_id: int
    scaffold: str
    activity: str | None


@dataclass(frozen=True)
class ScaffoldSplit:
    molecule_id: int
    scaffold: str | None
    partition: str | None


class ScaffoldSampler:
    """Sample a registered relation with unique BIGINT molecule IDs and SMILES."""

    def __init__(self, connection: DuckDBPyConnection, relation: str):
        if type(relation) is not str or not relation or "\0" in relation:
            raise ValueError("relation must be a nonempty identifier without NUL")
        self._db = connection
        self._relation = '"' + relation.replace('"', '""') + '"'
        self._scaffolds = f"""
            WITH input AS MATERIALIZED (
                SELECT molecule_id, smiles, activity FROM main.{self._relation} WHERE molecule_id IS NOT NULL
            ), scaffolds AS MATERIALIZED (
                SELECT smiles, murcko_scaffold(smiles) AS scaffold
                FROM (SELECT DISTINCT smiles FROM input WHERE smiles IS NOT NULL)
            )
        """

    @staticmethod
    def _validate_seed(seed):
        if type(seed) is not int or not -(1 << 63) <= seed < (1 << 63):
            raise ValueError("seed must be a signed 64-bit integer")

    def sample(self, *, per_stratum: int, seed: int) -> list[ScaffoldSample]:
        """Choose stable hash ranks by Murcko scaffold and activity class."""
        self._validate_seed(seed)
        if type(per_stratum) is not int or not 1 <= per_stratum <= 1_000_000:
            raise ValueError("per_stratum must be in [1, 1000000]")
        rows = self._db.execute(self._scaffolds + """
            , ranked AS (
                SELECT molecule_id, scaffold, activity, row_number() OVER (
                    PARTITION BY scaffold, activity
                    ORDER BY md5(?::VARCHAR || '/' || molecule_id::VARCHAR), molecule_id
                ) AS rank
                FROM input JOIN scaffolds USING (smiles) WHERE scaffold IS NOT NULL
            )
            SELECT molecule_id, scaffold, activity FROM ranked WHERE rank <= ?
            ORDER BY scaffold, activity NULLS LAST, rank
        """, [seed, per_stratum]).fetchall()
        return [ScaffoldSample(*row) for row in rows]

    def split(self, *, train: float, validation: float, seed: int) -> list[ScaffoldSplit]:
        """Keep each scaffold together, including the empty acyclic scaffold."""
        self._validate_seed(seed)
        if any(type(value) not in (int, float) or not 0 <= value < 1 for value in (train, validation)):
            raise ValueError("split fractions must be in [0, 1)")
        if train == 0 or train + validation >= 1:
            raise ValueError("train and test fractions must be positive")
        rows = self._db.execute(self._scaffolds + """
            , scores AS (
                SELECT DISTINCT scaffold,
                    ('0x' || substr(md5(?::VARCHAR || '/' || scaffold), 1, 8))::UBIGINT /
                      4294967296.0 AS score
                FROM scaffolds WHERE scaffold IS NOT NULL
            ), allocations AS (
                SELECT scaffold, CASE WHEN score < ? THEN 'train'
                   WHEN score < ? THEN 'validation' ELSE 'test' END AS partition
                FROM scores
            )
            SELECT molecule_id, scaffold, partition FROM input LEFT JOIN scaffolds USING (smiles)
            LEFT JOIN allocations USING (scaffold) ORDER BY molecule_id
        """, [seed, train, train + validation]).fetchall()
        return [ScaffoldSplit(*row) for row in rows]


class SparseFeatures:
    """Store nonzero and explicitly missing coordinates against an immutable feature dictionary."""

    def __init__(self, connection: DuckDBPyConnection, features: Mapping[str, int]):
        if not features or any(type(name) is not str or not name or type(identifier) is not int
                               or not 0 <= identifier < (1 << 63) for name, identifier in features.items()):
            raise ValueError("features map nonempty names to nonnegative signed 64-bit IDs")
        if len(set(features.values())) != len(features):
            raise ValueError("feature IDs must be unique")
        self._db = connection
        self._features = dict(features)
        with connection.cursor() as tx:
            tx.execute("BEGIN TRANSACTION")
            tx.execute("CREATE SCHEMA IF NOT EXISTS ducksmiles_sparse")
            tx.execute("""
                CREATE TABLE IF NOT EXISTS ducksmiles_sparse.dictionary (
                    name VARCHAR PRIMARY KEY, id BIGINT UNIQUE NOT NULL
                );
                CREATE TABLE IF NOT EXISTS ducksmiles_sparse.versions (version VARCHAR PRIMARY KEY);
                CREATE TABLE IF NOT EXISTS ducksmiles_sparse.molecules (
                    version VARCHAR, id BIGINT, PRIMARY KEY (version, id)
                );
                CREATE TABLE IF NOT EXISTS ducksmiles_sparse.coordinates (
                    version VARCHAR, id BIGINT, feature BIGINT, value DOUBLE,
                    PRIMARY KEY (version, id, feature), CHECK (value IS NULL OR isfinite(value))
                );
            """)
            existing = dict(tx.execute("SELECT name, id FROM ducksmiles_sparse.dictionary").fetchall())
            if existing and existing != self._features:
                raise ValueError("feature dictionary does not match the existing IDs")
            if not existing:
                tx.execute("""
                    INSERT INTO ducksmiles_sparse.dictionary
                    SELECT unnest(?::VARCHAR[]), unnest(?::BIGINT[])
                """, [list(self._features), list(self._features.values())])
            tx.execute("COMMIT")

    @staticmethod
    def _validate(version, molecules):
        if type(version) is not str or not version:
            raise ValueError("supply a nonempty version")
        if (type(molecules) is not list or any(type(key) is not int or not -(1 << 63) <= key < 1 << 63
                                               for key in molecules)
                or len(set(molecules)) != len(molecules)):
            raise ValueError("provide unique signed 64-bit molecule IDs")

    def store(self, version: str, molecules: list[int], terms: list[tuple[int, str, float | None]]) -> int:
        """Write one immutable matrix. A missing contribution marks its coordinate missing."""
        self._validate(version, molecules)
        valid_molecules = set(molecules)
        if type(terms) is not list:
            raise TypeError("expected coordinate terms")
        rows = []
        for term in terms:
            if type(term) is not tuple or len(term) != 3:
                raise TypeError("expected (molecule ID, feature name, value) terms")
            molecule, name, value = term
            if type(molecule) is not int or molecule not in valid_molecules:
                raise KeyError(molecule)
            feature = self._features[name]
            if value is not None and (type(value) not in (int, float) or not -float_info.max <= value <= float_info.max):
                raise ValueError("coordinate values must be finite numbers or None")
            rows.append({"molecule": molecule, "feature": feature, "value": value})
        with self._db.cursor() as tx:
            tx.execute("BEGIN TRANSACTION")
            tx.execute("INSERT INTO ducksmiles_sparse.versions VALUES (?)", [version])
            tx.execute("INSERT INTO ducksmiles_sparse.molecules SELECT ?, unnest(?::BIGINT[])", [version, molecules])
            count, = tx.execute("""
                INSERT INTO ducksmiles_sparse.coordinates
                WITH terms AS (
                    SELECT t.molecule AS id, t.feature AS feature, t.value AS value
                    FROM unnest(?::STRUCT(molecule BIGINT, feature BIGINT, value DOUBLE)[]) AS input(t)
                ), merged AS (
                    SELECT id, feature, CASE WHEN count(value) = count(*) THEN fsum(value) END AS value
                    FROM terms GROUP BY id, feature
                )
                SELECT ?, id, feature, value FROM merged WHERE value IS DISTINCT FROM 0
            """, [rows, version]).fetchone()
            tx.execute("COMMIT")
        return count

    def _require_molecules(self, connection, version, molecules):
        self._validate(version, molecules)
        exists = connection.execute("SELECT version FROM ducksmiles_sparse.versions WHERE version = ?",
                                    [version]).fetchone()
        if exists is None:
            raise KeyError(version)
        missing = connection.execute("""
            SELECT id FROM unnest(?::BIGINT[]) AS input(id)
            ANTI JOIN (SELECT id FROM ducksmiles_sparse.molecules WHERE version = ?) USING (id)
            LIMIT 1
        """, [molecules, version]).fetchone()
        if missing is not None:
            raise KeyError(missing[0])

    def coordinates(self, version: str, molecules: list[int]) -> list[tuple[int, int, float | None]]:
        """Return stored coordinates in molecule and feature ID order."""
        with self._db.cursor() as db:
            self._require_molecules(db, version, molecules)
            return db.execute("""
                SELECT id, feature, value FROM ducksmiles_sparse.coordinates
                WHERE version = ? AND id IN (SELECT unnest(?::BIGINT[])) ORDER BY id, feature
            """, [version, molecules]).fetchall()

    def dense(self, version: str, molecules: list[int], names: list[str]) -> list[tuple[int, tuple[float | None, ...]]]:
        """Read at most 262144 cells, retaining request order and explicit missing values."""
        if type(names) is not list or any(type(name) is not str for name in names) or len(set(names)) != len(names):
            raise ValueError("provide a list of unique feature names")
        features = [self._features[name] for name in names]
        with self._db.cursor() as db:
            self._require_molecules(db, version, molecules)
            if len(molecules) * len(features) > 262_144:
                raise ValueError("a dense read may contain at most 262144 cells")
            if not features:
                return [(molecule, ()) for molecule in molecules]
            rows = db.execute("""
                WITH requested_rows AS (SELECT * FROM unnest(?::BIGINT[]) WITH ORDINALITY AS t(id, position)),
                     requested_features AS (
                         SELECT * FROM unnest(?::BIGINT[]) WITH ORDINALITY AS t(feature, position)
                     ), values_at_version AS (
                         SELECT id, feature, value FROM ducksmiles_sparse.coordinates WHERE version = ?
                     )
                SELECT r.id, list(CASE WHEN v.id IS NULL THEN 0 ELSE v.value END ORDER BY f.position)
                FROM requested_rows r CROSS JOIN requested_features f LEFT JOIN values_at_version v
                  ON (r.id, f.feature) = (v.id, v.feature)
                GROUP BY r.id, r.position ORDER BY r.position
            """, [molecules, features, version]).fetchall()
        return [(molecule, tuple(values)) for molecule, values in rows]


class FingerprintIndex:
    """Store Morgan fingerprints for explicit, host-provided generations."""

    def __init__(self, connection: DuckDBPyConnection):
        self._db = connection
        connection.execute("CREATE SCHEMA IF NOT EXISTS ducksmiles_fingerprints")
        connection.execute("""
            CREATE TABLE IF NOT EXISTS ducksmiles_fingerprints.generations (
                generation VARCHAR, radius INTEGER, width INTEGER, revision BIGINT NOT NULL,
                PRIMARY KEY (generation, radius, width)
            )
        """)
        connection.execute("""
            CREATE TABLE IF NOT EXISTS ducksmiles_fingerprints.entries (
                generation VARCHAR, radius INTEGER, width INTEGER, id BIGINT,
                bits BLOB NOT NULL, popcount INTEGER NOT NULL,
                PRIMARY KEY (generation, radius, width, id),
                CHECK (octet_length(bits) = (width + 7) // 8 AND popcount = bit_count(bits::BIT))
            )
        """)
        connection.execute("""
            CREATE TABLE IF NOT EXISTS ducksmiles_fingerprints.graphs (
                version VARCHAR PRIMARY KEY, generation VARCHAR NOT NULL, radius INTEGER NOT NULL,
                width INTEGER NOT NULL, fingerprint_revision BIGINT NOT NULL, threshold DOUBLE NOT NULL,
                method VARCHAR NOT NULL CHECK (method = 'exact')
            )
        """)
        connection.execute("""
            CREATE TABLE IF NOT EXISTS ducksmiles_fingerprints.vertices (
                version VARCHAR, id BIGINT, PRIMARY KEY (version, id)
            )
        """)
        connection.execute("""
            CREATE TABLE IF NOT EXISTS ducksmiles_fingerprints.edges (
                version VARCHAR, src BIGINT, dst BIGINT, score DOUBLE NOT NULL,
                PRIMARY KEY (version, src, dst), CHECK (src < dst)
            )
        """)

    @staticmethod
    def _definition(generation, definition):
        if type(generation) is not str or not generation:
            raise ValueError("generation must be a nonempty string")
        if type(definition) is not MorganFingerprint:
            raise TypeError("expected a MorganFingerprint")
        return [generation, definition.radius, definition.bits]

    def replace(self, generation: str, definition: MorganFingerprint,
                molecules: list[tuple[int, str | None]]) -> int:
        """Replace one fingerprint generation atomically and return the valid row count."""
        key = self._definition(generation, definition)
        if type(molecules) is not list:
            raise TypeError("expected a list of (id, SMILES) tuples")
        rows = []
        seen = set()
        for record in molecules:
            if type(record) is not tuple or len(record) != 2:
                raise TypeError("expected (id, SMILES) tuples")
            molecule_id, smiles = record
            if type(molecule_id) is not int or not -(1 << 63) <= molecule_id < (1 << 63):
                raise ValueError("molecule IDs must be signed 64-bit integers")
            if molecule_id in seen:
                raise ValueError("duplicate molecule ID")
            if smiles is not None and type(smiles) is not str:
                raise TypeError("SMILES must be a string or None")
            seen.add(molecule_id)
            rows.append({"id": molecule_id, "smiles": smiles})

        with self._db.cursor() as tx:
            tx.execute("BEGIN TRANSACTION")
            tx.execute("""
                INSERT INTO ducksmiles_fingerprints.generations VALUES (?, ?, ?, 1)
                ON CONFLICT (generation, radius, width) DO UPDATE SET revision = revision + 1
            """, key)
            tx.execute(
                "DELETE FROM ducksmiles_fingerprints.entries WHERE (generation, radius, width) = (?, ?, ?)", key
            )
            inserted, = tx.execute("""
                INSERT INTO ducksmiles_fingerprints.entries
                WITH parsed AS MATERIALIZED (
                    SELECT m.id AS id, morgan_fp_bits(m.smiles, ?, ?) AS bits
                    FROM unnest(?::STRUCT(id BIGINT, smiles VARCHAR)[]) AS t(m)
                )
                SELECT ?, ?, ?, id, bits, bit_count(bits::BIT) FROM parsed WHERE bits IS NOT NULL
            """, [definition.radius, definition.bits, rows, *key]).fetchone()
            tx.execute("COMMIT")
        return inserted

    def search(self, generation: str, definition: MorganFingerprint, smiles: str, *,
               threshold: float = 0.0, limit: int = 10) -> SimilarityResult:
        self._definition(generation, definition)
        self._validate_search(threshold, limit)
        if type(smiles) is not str:
            raise TypeError("SMILES must be a string")
        bits, = self._db.execute("SELECT morgan_fp_bits(?, ?, ?)", [smiles, definition.radius, definition.bits]).fetchone()
        if bits is None:
            return SimilarityResult((), 0)
        return self.search_bits(generation, definition, bits, threshold=threshold, limit=limit)

    @staticmethod
    def _validate_search(threshold, limit):
        if type(threshold) not in (int, float) or not 0 <= threshold <= 1:
            raise ValueError("threshold must be in [0, 1]")
        if type(limit) is not int or not 1 <= limit <= 1_000_000:
            raise ValueError("limit must be an integer in [1, 1000000]")

    def search_bits(self, generation: str, definition: MorganFingerprint, bits: bytes, *,
                    threshold: float = 0.0, limit: int = 10) -> SimilarityResult:
        key = self._definition(generation, definition)
        self._validate_search(threshold, limit)
        if type(bits) is not bytes or len(bits) != (definition.bits + 7) // 8:
            raise ValueError("query fingerprint byte length does not match its definition")
        if definition.bits % 8 and bits[-1] >> (definition.bits % 8):
            raise ValueError("query fingerprint has bits beyond its declared width")
        popcount = int.from_bytes(bits, "little").bit_count()
        hits, candidates = self._db.execute("""
            WITH candidates AS MATERIALIZED (
                SELECT id, bits FROM ducksmiles_fingerprints.entries
                CROSS JOIN (SELECT ?::INTEGER AS query_count, ?::DOUBLE AS threshold)
                WHERE (generation, radius, width) = (?, ?, ?)
                  AND coalesce(least(popcount, query_count)::DOUBLE /
                               nullif(greatest(popcount, query_count), 0), 0) >= threshold
            ), top_hits AS (
                SELECT id, tanimoto_bit(bits, ?::BLOB) AS score FROM candidates
                WHERE score >= ? ORDER BY score DESC, id LIMIT ?
            )
            SELECT list(struct_pack(id := id, score := score) ORDER BY score DESC, id),
                   (SELECT count(*) FROM candidates)
            FROM top_hits
        """, [popcount, threshold, *key, bits, threshold, limit]).fetchone()
        return SimilarityResult(tuple(SimilarityHit(h["id"], h["score"]) for h in hits or ()), candidates)

    def build_graph(self, generation: str, definition: MorganFingerprint, *, version: str,
                    threshold: float, block_size: int = 256) -> SimilarityGraphSummary:
        """Materialize a new immutable exact graph from a single fingerprint snapshot."""
        key = self._definition(generation, definition)
        self._validate_search(threshold, 1)
        if type(version) is not str or not version:
            raise ValueError("graph version must be a nonempty string")
        if type(block_size) is not int or not 1 <= block_size <= 1024:
            raise ValueError("block_size must be an integer in [1, 1024]")
        with self._db.cursor() as tx:
            tx.execute("BEGIN TRANSACTION")
            if tx.execute("SELECT 1 FROM ducksmiles_fingerprints.graphs WHERE version = ?", [version]).fetchone():
                raise ValueError("graph version already exists")
            metadata = tx.execute("""
                SELECT revision FROM ducksmiles_fingerprints.generations
                WHERE (generation, radius, width) = (?, ?, ?)
            """, key).fetchone()
            if metadata is None:
                raise KeyError(generation)
            tx.execute("INSERT INTO ducksmiles_fingerprints.graphs VALUES (?, ?, ?, ?, ?, ?, 'exact')",
                       [version, *key, metadata[0], threshold])
            vertex_count, = tx.execute("""
                INSERT INTO ducksmiles_fingerprints.vertices
                SELECT ?, id FROM ducksmiles_fingerprints.entries WHERE (generation, radius, width) = (?, ?, ?)
            """, [version, *key]).fetchone()
            blocks = tx.execute("""
                WITH ranked AS (
                    SELECT id, (row_number() OVER (ORDER BY id) - 1) // ? AS block
                    FROM ducksmiles_fingerprints.vertices WHERE version = ?
                )
                SELECT min(id), max(id) FROM ranked GROUP BY block ORDER BY block
            """, [block_size, version]).fetchall()
            edge_count = 0
            for left_index, left in enumerate(blocks):
                for right in blocks[left_index:]:
                    count, = tx.execute("""
                        INSERT INTO ducksmiles_fingerprints.edges
                        WITH fingerprints AS NOT MATERIALIZED (
                            SELECT id, bits FROM ducksmiles_fingerprints.entries
                            WHERE (generation, radius, width) = (?, ?, ?)
                        ), scored AS MATERIALIZED (
                            SELECT a.id AS src, b.id AS dst, tanimoto_bit(a.bits, b.bits) AS score
                            FROM fingerprints AS a JOIN fingerprints AS b ON a.id < b.id
                            WHERE a.id BETWEEN ? AND ? AND b.id BETWEEN ? AND ?
                        )
                        SELECT ?, src, dst, score FROM scored WHERE score >= ?
                    """, [*key, *left, *right, version, threshold]).fetchone()
                    edge_count += count
            tx.execute("COMMIT")
        return SimilarityGraphSummary(vertex_count, vertex_count * (vertex_count - 1) // 2, edge_count)

    def graph_vertices(self, version: str) -> list[int]:
        self._graph_threshold(version)
        rows = self._db.execute("""
            SELECT id FROM ducksmiles_fingerprints.vertices WHERE version = ? ORDER BY id
        """, [version]).fetchall()
        return [row[0] for row in rows]

    def graph_edges(self, version: str, *, threshold: float) -> list[SimilarityEdge]:
        self._validate_search(threshold, 1)
        if threshold < self._graph_threshold(version):
            raise ValueError("requested threshold is below the stored graph threshold")
        rows = self._db.execute("""
            SELECT src, dst, score FROM ducksmiles_fingerprints.edges
            WHERE version = ? AND score >= ? ORDER BY src, dst
        """, [version, threshold]).fetchall()
        return [SimilarityEdge(*row) for row in rows]

    def _graph_threshold(self, version):
        if type(version) is not str or not version:
            raise ValueError("graph version must be a nonempty string")
        row = self._db.execute("SELECT threshold FROM ducksmiles_fingerprints.graphs WHERE version = ?",
                               [version]).fetchone()
        if row is None:
            raise KeyError(version)
        return row[0]


@dataclass(frozen=True)
class DaltonMass:
    value: float

    def __post_init__(self):
        if type(self.value) not in (float, int) or not 0 < self.value <= float_info.max:
            raise ValueError("mass must be a finite positive number of daltons")


@dataclass(frozen=True)
class MassMatch:
    reference_index: int
    peak_index: int
    error_ppm: float


class MassWindowSearch:
    """Trusted embedded operation for neutral masses expressed in daltons."""

    def __init__(self, connection: DuckDBPyConnection):
        self._db = connection

    def match(self, references: tuple[DaltonMass, ...], peaks: tuple[DaltonMass, ...],
              tolerance_ppm: float) -> list[MassMatch]:
        """Return inclusive matches in input order, with zero-based indices."""
        if type(tolerance_ppm) not in (int, float) or not 0 <= tolerance_ppm <= 1_000_000:
            raise ValueError("tolerance must be between 0 and 1,000,000 ppm")
        for masses in (references, peaks):
            if type(masses) is not tuple or any(type(mass) is not DaltonMass for mass in masses):
                raise TypeError("expected a tuple of DaltonMass values")
        rows = self._db.execute("""
            WITH reference_masses AS (
                SELECT mass, ord, mass * (?::DOUBLE / 1000000) AS delta
                FROM unnest(?::DOUBLE[]) WITH ORDINALITY AS t(mass, ord)
            ), windows AS (
                SELECT mass, ord, mass - delta AS lo, mass + delta AS hi FROM reference_masses
            ), peaks AS (
                SELECT mass, ord FROM unnest(?::DOUBLE[]) WITH ORDINALITY AS t(mass, ord)
            )
            SELECT w.ord - 1, p.ord - 1, (p.mass - w.mass) / w.mass * 1000000
            FROM windows w JOIN peaks p ON p.mass BETWEEN w.lo AND w.hi
            ORDER BY w.ord, p.ord
        """, [tolerance_ppm, [m.value for m in references], [m.value for m in peaks]]).fetchall()
        return [MassMatch(*row) for row in rows]


class ConcentrationBasis(Enum):
    MASS = "g/L"
    MOLAR = "mol/L"


@dataclass(frozen=True)
class NormalizedConcentration:
    measurement_id: int
    source_value: float | None
    source_unit: str | None
    value: float | None
    unit: str
    rule_version: str
    reason: str | None


class ConcentrationTable:
    """Convert a registered measurement relation to a compatible concentration basis."""

    RULE_VERSION = "concentration-1"

    def __init__(self, connection: DuckDBPyConnection, relation: str):
        self._db = connection
        if type(relation) is not str or not relation or "\0" in relation:
            raise ValueError("relation must be a nonempty identifier without NUL")
        quoted = '"' + relation.replace('"', '""') + '"'
        self._normalized = f"""
            WITH factors(unit_name, base, scale) AS (
                VALUES ('g/L', 'g/L', 1.0), ('mg/L', 'g/L', 0.001), ('ug/mL', 'g/L', 0.001),
                       ('µg/mL', 'g/L', 0.001), ('µg/L', 'g/L', 0.000001),
                       ('mol/L', 'mol/L', 1.0), ('mmol/L', 'mol/L', 0.001),
                       ('umol/L', 'mol/L', 0.000001), ('µmol/L', 'mol/L', 0.000001)
            ), desired AS (SELECT ?::VARCHAR AS target),
            readings AS (
                SELECT id, value AS source_value, unit, target,
                       CASE WHEN value IS NULL THEN 'missing_value'
                            WHEN NOT isfinite(value) THEN 'nonfinite_value'
                            WHEN scale IS NULL THEN 'unknown_unit'
                            WHEN base <> target AND
                                (grams_per_mole IS NULL OR NOT isfinite(grams_per_mole) OR grams_per_mole <= 0)
                                THEN 'molar_mass_required' END AS problem,
                       CASE WHEN base = target THEN value * scale
                            WHEN target = 'g/L' THEN value * scale * grams_per_mole
                            ELSE value * scale / nullif(grams_per_mole, 0) END AS converted
                FROM main.{quoted} LEFT JOIN factors ON unit = unit_name CROSS JOIN desired
            ), normalized AS (
                SELECT id, source_value, unit, target,
                       CASE WHEN problem IS NULL AND isfinite(converted) THEN converted END AS result,
                       CASE WHEN problem IS NOT NULL THEN problem
                            WHEN NOT isfinite(converted) THEN 'nonfinite_result' END AS problem
                FROM readings
            )
        """

    @staticmethod
    def _basis(basis):
        if type(basis) is not ConcentrationBasis:
            raise TypeError("expected a ConcentrationBasis")
        return basis.value

    def read(self, basis: ConcentrationBasis) -> list[NormalizedConcentration]:
        rows = self._db.execute(self._normalized + """
            SELECT id, source_value, unit, result, target, problem FROM normalized ORDER BY id
        """, [self._basis(basis)]).fetchall()
        return [NormalizedConcentration(i, value, unit, result, target, self.RULE_VERSION, problem)
                for i, value, unit, result, target, problem in rows]

    def summarize(self, basis: ConcentrationBasis) -> tuple[int, int, float | None]:
        """Return valid count, excluded count and total without exporting raw rows."""
        return self._db.execute(self._normalized + """
            SELECT count(result), count(*) - count(result), fsum(result) FROM normalized
        """, [self._basis(basis)]).fetchone()
