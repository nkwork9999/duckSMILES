import os
import tempfile
import unittest
from pathlib import Path

import duckdb

from ducksmiles_embedded import FingerprintIndex, MorganFingerprint, ScaffoldSampler, SimilarityEdge, SimilarityHit


class FingerprintTests(unittest.TestCase):
    def setUp(self):
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.path = str(Path(directory.name) / "molecules.duckdb")
        self.db = self.connect()
        self.addCleanup(self.db.close)
        self.index = FingerprintIndex(self.db)
        self.definition = MorganFingerprint(radius=2, bits=128)

    def connect(self):
        extension = os.getenv("DUCKSMILES_TEST_EXTENSION")
        db = duckdb.connect(self.path, config={"allow_unsigned_extensions": bool(extension)})
        db.load_extension(extension or "ducksmiles")
        return db

    def test_results_match_complete_sort_and_include_exact_threshold(self):
        molecules = [(4, "CCO"), (1, "CCO"), (3, "CCN"), (2, "O"), (5, None), (6, "invalid")]
        self.assertEqual(self.index.replace("batch-1", self.definition, molecules), 4)
        for other in (MorganFingerprint(3, 128), MorganFingerprint(2, 256)):
            self.index.replace("batch-1", other, [(99, "CCO")])
        self.index.replace("batch-2", self.definition, [(98, "CCO")])
        expected = self.db.execute("""
            SELECT id, tanimoto_bit(morgan_fp_bits(smiles, 2, 128),
                                   morgan_fp_bits('CCO', 2, 128)) score
            FROM unnest(?::STRUCT(id BIGINT, smiles VARCHAR)[]) AS t(m)
            CROSS JOIN LATERAL (SELECT m.id, m.smiles)
            WHERE score IS NOT NULL ORDER BY score DESC, id
        """, [[{"id": i, "smiles": s} for i, s in molecules]]).fetchall()
        for threshold in (0.0, expected[-2][1], 0.5, 0.8, 0.95, 1.0):
            for limit in (1, 20):
                result = self.index.search("batch-1", self.definition, "CCO", threshold=threshold, limit=limit)
                self.assertEqual(result.hits, tuple(SimilarityHit(i, score)
                    for i, score in expected if score >= threshold)[:limit])
                self.assertLessEqual(result.candidates_scored, len(expected))
        self.assertEqual(self.index.search("batch-1", self.definition, "invalid").hits, ())

    def test_replacement_is_atomic_and_visible_after_reopen(self):
        self.index.replace("batch", self.definition, [(1, "CCO")])
        with self.assertRaises(ValueError):
            self.index.replace("batch", self.definition, [(2, "CCO"), (2, "O")])
        with self.connect() as reader:
            self.assertEqual(FingerprintIndex(reader).search("batch", self.definition, "CCO").hits,
                             (SimilarityHit(1, 1.0),))
        self.assertEqual(self.index.replace("batch", self.definition, []), 0)
        self.assertEqual(self.index.search("batch", self.definition, "CCO").hits, ())

    def test_empty_bitsets_follow_native_zero_similarity(self):
        self.index.replace("batch", self.definition, [(1, "CCO")])
        self.db.execute("""
            UPDATE ducksmiles_fingerprints.entries SET bits = ?, popcount = 0
            WHERE generation = 'batch'
        """, [bytes(16)])
        all_hits = self.index.search_bits("batch", self.definition, bytes(16), threshold=0)
        self.assertEqual(all_hits.hits, (SimilarityHit(1, 0.0),))
        above_zero = self.index.search_bits("batch", self.definition, bytes(16), threshold=0.001)
        self.assertEqual(above_zero.hits, ())
        self.assertEqual(above_zero.candidates_scored, 0)

    def test_graph_compares_upper_triangle_and_retains_isolated_vertices(self):
        self.index.replace("sample", self.definition, [(9, "CCO"), (17, "CCO"), (3, "CCO"),
                                                        (-5, "O"), (4, "CCN"), (8, None)])
        built = self.index.build_graph("sample", self.definition, version="g1", threshold=0.9, block_size=2)
        self.assertEqual((built.vertices, built.comparisons, built.edges), (5, 10, 3))
        self.assertEqual(self.index.graph_vertices("g1"), [-5, 3, 4, 9, 17])
        self.assertEqual(self.index.graph_edges("g1", threshold=1.0),
                         [SimilarityEdge(3, 9, 1.0), SimilarityEdge(3, 17, 1.0), SimilarityEdge(9, 17, 1.0)])
        with self.assertRaises(ValueError):
            self.index.graph_edges("g1", threshold=0.5)
        with self.assertRaises(ValueError):
            self.index.build_graph("sample", self.definition, version="g1", threshold=0.9)
        self.index.replace("sample", self.definition, [(1, "CCO")])
        self.assertEqual(self.index.graph_vertices("g1"), [-5, 3, 4, 9, 17])
        self.assertEqual(len(self.index.graph_edges("g1", threshold=0.9)), 3)

    def test_scaffold_sampler_is_stable_across_thread_counts(self):
        self.db.execute("""
            CREATE TABLE library AS SELECT i::BIGINT AS molecule_id,
                CASE i % 2 WHEN 0 THEN 'CCO' ELSE 'c1ccccc1' END AS smiles,
                (i % 7)::VARCHAR AS activity
            FROM range(10000) t(i)
        """)
        sampler = ScaffoldSampler(self.db, "library")
        self.db.execute("SET threads = 1")
        first = sampler.sample(per_stratum=5, seed=87)
        self.db.execute("SET threads = 4")
        repeated = sampler.sample(per_stratum=5, seed=87)
        self.assertEqual(first, repeated)
        self.assertEqual(len(first), 70)
        self.assertEqual(len({(row.scaffold, row.activity) for row in first}), 14)
        self.assertNotEqual(first, sampler.sample(per_stratum=5, seed=88))
