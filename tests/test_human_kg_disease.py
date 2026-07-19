import io
import json
import tempfile
import unittest
import zipfile
from pathlib import Path
from unittest import mock

import numpy as np
import pandas as pd

from scripts import human_kg_disease as kg
from scripts import public_evidence


class HumanKgDiseaseTest(unittest.TestCase):
    genes = ["A", "B", "C", "D"]

    @staticmethod
    def _write_csv(path: Path, rows: list[dict]):
        pd.DataFrame(rows).to_csv(path, index=False)

    @staticmethod
    def _write_chea(path: Path, down: bool = False):
        relation = "downregulates" if down else "upregulates"
        frame = pd.DataFrame(
            {
                "source": ["source"],
                "relation": [relation],
                "target": ["target"],
                "z_score": [2.0],
                "p_value": [0.01],
                "source_label": ["A"],
                "target_label": ["B"],
            }
        )
        with zipfile.ZipFile(path, "w") as archive:
            archive.writestr(
                f"Transcription Factor.{relation}.Transcription Factor.edges.csv",
                frame.to_csv(index=False),
            )

    def _write_cache(self, root: Path):
        assets = root / "assets"
        assets.mkdir()
        term_rows = [
            {
                "source": "genetic",
                "relation": "association",
                "target": "A",
                "source_label": "Synthetic genetic disease",
                "target_label": "A",
            },
            {
                "source": "genetic",
                "relation": "association",
                "target": "B",
                "source_label": "Synthetic genetic disease",
                "target_label": "B",
            },
            {
                "source": "sepsis",
                "relation": "association",
                "target": "C",
                "source_label": "Synthetic sepsis",
                "target_label": "C",
            },
            {
                "source": "sepsis",
                "relation": "association",
                "target": "D",
                "source_label": "Synthetic sepsis",
                "target_label": "D",
            },
        ]
        manifest_assets = {}
        for asset_id in ["disgenet", "jensen_diseases", "gwas_catalog"]:
            path = assets / f"{asset_id}.csv"
            self._write_csv(path, term_rows)
            manifest_assets[asset_id] = {"path": str(path.relative_to(root))}

        archs = assets / "archs4_tf.csv"
        self._write_csv(
            archs,
            [
                {
                    "source": "archs-human",
                    "relation": "ARCHS4_TFs_Coexpression",
                    "target": "B",
                    "source_label": "A human tf ARCHS4 coexpression",
                    "target_label": "B",
                },
                {
                    "source": "archs-mouse",
                    "relation": "ARCHS4_TFs_Coexpression",
                    "target": "D",
                    "source_label": "A mouse tf ARCHS4 coexpression",
                    "target_label": "D",
                },
            ],
        )
        manifest_assets["archs4_tf"] = {"path": str(archs.relative_to(root))}

        ppi = assets / "string_ppi.csv"
        self._write_csv(
            ppi,
            [
                {
                    "source": "A",
                    "relation": "STRING-db PPI",
                    "target": "B",
                    "source_label": "A",
                    "target_label": "B",
                    "combined_score": 900,
                }
            ],
        )
        manifest_assets["string_ppi"] = {"path": str(ppi.relative_to(root))}

        chea = assets / "chea_node_draw.zip"
        self._write_chea(chea)
        manifest_assets["chea_node_draw"] = {"path": str(chea.relative_to(root))}
        (root / kg.MONDO_TERMS_FILE).write_text(
            json.dumps(
                {
                    "terms": [
                        {
                            "id": "MONDO:GENETIC",
                            "name": "genetic disease",
                            "synonyms": [],
                            "parent_ids": [],
                            "family": "genetic disease",
                        },
                        {
                            "id": "MONDO:0001",
                            "name": "Synthetic genetic disease",
                            "synonyms": [],
                            "parent_ids": ["MONDO:GENETIC"],
                            "family": "genetic disease",
                        },
                        {
                            "id": "MONDO:0002",
                            "name": "Synthetic sepsis",
                            "synonyms": [],
                            "parent_ids": [],
                            "family": "infectious disease",
                        },
                    ]
                }
            )
        )
        (root / kg.MANIFEST_FILE).write_text(
            json.dumps(
                {
                    "schema_version": 1,
                    "species": "human",
                    "assets": manifest_assets,
                    "mondo": {"terms_path": kg.MONDO_TERMS_FILE},
                    "network_primary_validation": {"status": "not_run", "approved": False},
                }
            )
        )

    def _write_evidence_cache(self, root: Path):
        public_evidence.write_evidence_cache(
            root,
            [
                {
                    "record_id": "human-perturbation",
                    "resource": "perturbseqr",
                    "species": "human",
                    "species_provenance": "library_name",
                    "kind": "gene_perturbation",
                    "up_genes": ["A", "B"],
                    "down_genes": [],
                    "perturbation_gene": "GENEZ",
                    "perturbation_mode": "KO",
                    "library_name": "Perturb Atlas Human",
                    "source_label": "GENEZ KO up",
                }
            ],
        )

    def _build(self, kg_root: Path, evidence_root: Path):
        signed_genes = pd.Series({"A": 1.0, "B": 0.6, "C": 0.2, "D": -0.1})
        signed_programs = pd.Series({"P0": 1.0})
        dspin = pd.DataFrame(
            [[0.0, 0.0, 0.0, 0.0], [1.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0]],
            index=self.genes,
            columns=self.genes,
        )
        gene_to_program = pd.DataFrame({"P0": [1.0, 0.5, 0.0, 0.0]}, index=self.genes)
        return kg.build_human_kg_hypotheses(
            kg_root,
            evidence_root,
            signed_genes,
            signed_programs,
            dspin,
            gene_to_program,
            pd.DataFrame({"cell_type": ["Monocyte"]}),
            permutations=20,
        )

    def test_human_only_cache_separates_disease_mimic_and_mechanism_outputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            kg_root = root / "kg"
            evidence_root = root / "evidence"
            kg_root.mkdir()
            self._write_cache(kg_root)
            self._write_evidence_cache(evidence_root)
            result = self._build(kg_root, evidence_root)

        self.assertEqual(result["human_diseases"].iloc[0]["disease_name"], "Synthetic genetic disease")
        self.assertEqual(result["human_diseases"].iloc[0]["hypothesis_class"], "Monogenic/genetic disease")
        self.assertEqual(result["human_mimics"].iloc[0]["disease_name"], "Synthetic sepsis")
        self.assertEqual(result["human_mimics"].iloc[0]["hypothesis_class"], "Infection/sepsis-like mimic")
        self.assertIn("GENEZ", set(result["human_mechanisms"]["perturbation_gene"]))
        self.assertFalse(result["network_concordance"].empty)
        self.assertFalse(result["manifold_modules"].empty)

    def test_explicit_human_filter_rejects_mouse_and_unknown_terms(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "trrust.csv"
            self._write_csv(
                path,
                [
                    {"source": "h", "source_label": "A human", "target_label": "B"},
                    {"source": "m", "source_label": "A mouse", "target_label": "B"},
                    {"source": "u", "source_label": "A regulator", "target_label": "B"},
                ],
            )
            terms, stats = kg.load_term_sets(path, "trrust", "regulatory", pd.Index(self.genes), "label_human")

        self.assertEqual([term["term_label"] for term in terms], ["A human"])
        self.assertEqual(stats["rejected_species_edges"], 2)

    def test_chea_species_labels_reject_mouse_mixed_and_unknown_edges(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "chea.zip"
            frame = pd.DataFrame(
                {
                    "source": ["a", "c", "d", "b"],
                    "relation": ["upregulates"] * 4,
                    "target": ["b", "d", "a", "c"],
                    "z_score": [1.0] * 4,
                    "source_label": ["A", "C", "D", "B"],
                    "target_label": ["B", "D", "A", "C"],
                    "species": ["human", "mouse", "mixed", "unknown"],
                }
            )
            with zipfile.ZipFile(path, "w") as archive:
                archive.writestr(
                    "Transcription Factor.upregulates.Transcription Factor.edges.csv",
                    frame.to_csv(index=False),
                )
            matrix = kg._load_chea_matrix(path, self.genes)

        self.assertGreater(matrix[1, 0], 0.0)
        self.assertEqual(matrix[3, 2], 0.0)
        self.assertEqual(matrix[0, 3], 0.0)
        self.assertEqual(matrix[2, 1], 0.0)

    def test_directed_sign_and_orientation_are_not_treated_as_unsigned(self):
        seed = np.array([1.0, 0.0])
        forward = np.array([[0.0, 0.0], [1.0, 0.0]])  # A -> B
        reversed_edge = np.array([[0.0, 1.0], [0.0, 0.0]])  # B -> A
        negative = np.array([[0.0, 0.0], [-1.0, 0.0]])

        self.assertEqual(kg.propagate_signed(seed, forward)[1][1], 1.0)
        self.assertEqual(kg.propagate_signed(seed, reversed_edge)[1][1], 0.0)
        self.assertEqual(kg.propagate_signed(seed, negative)[1][1], -1.0)

    def test_ppi_is_symmetric_and_manifold_is_deterministic(self):
        with tempfile.TemporaryDirectory() as tmp:
            ppi_path = Path(tmp) / "ppi.csv"
            self._write_csv(
                ppi_path,
                [
                    {
                        "source": "A",
                        "relation": "STRING-db PPI",
                        "target": "B",
                        "source_label": "A",
                        "target_label": "B",
                        "combined_score": 900,
                    }
                ],
            )
            ppi = kg._load_ppi_matrix(ppi_path, self.genes)

        self.assertEqual(ppi[0, 1], ppi[1, 0])
        coexpression = kg._archs4_tf_coexpression_matrix(
            [{"term_label": "A human tf ARCHS4 coexpression", "genes": {"B"}}],
            self.genes,
        )
        self.assertEqual(coexpression[0, 1], coexpression[1, 0])
        ppi_before = ppi.copy()
        forward = np.zeros((4, 4))
        forward[1, 0] = 1.0
        reversed_edge = np.zeros((4, 4))
        reversed_edge[0, 1] = 1.0
        self.assertNotEqual(kg.propagate_signed(np.array([1.0, 0.0, 0.0, 0.0]), forward)[1][1], kg.propagate_signed(np.array([1.0, 0.0, 0.0, 0.0]), reversed_edge)[1][1])
        np.testing.assert_allclose(ppi, ppi_before)
        first = kg._multi_view_embeddings(forward, forward, ppi, ppi)
        second = kg._multi_view_embeddings(forward, forward, ppi, ppi)
        np.testing.assert_allclose(first[0], second[0])
        np.testing.assert_allclose(first[1], second[1])

    def test_cache_only_execution_does_not_read_the_network(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            kg_root = root / "kg"
            evidence_root = root / "evidence"
            kg_root.mkdir()
            self._write_cache(kg_root)
            self._write_evidence_cache(evidence_root)
            with mock.patch("scripts.human_kg_disease.urllib.request.urlopen", side_effect=AssertionError("network read")):
                result = self._build(kg_root, evidence_root)

        self.assertFalse(result["human_diseases"].empty)

    def test_empty_manifest_recovers_known_local_assets_without_network(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self._write_cache(root)
            manifest_path = root / kg.MANIFEST_FILE
            manifest = json.loads(manifest_path.read_text())
            manifest["assets"] = {}
            manifest_path.write_text(json.dumps(manifest))
            with mock.patch("scripts.human_kg_disease.urllib.request.urlopen", side_effect=AssertionError("network read")):
                recovered, warnings = kg.read_kg_cache(root)

        self.assertIn("disgenet", recovered["assets"])
        self.assertTrue(any("recovered" in warning.lower() for warning in warnings))

    def test_bootstrap_gate_requires_positive_confidence_interval(self):
        approved = kg.bootstrap_mrr_improvement([4, 4, 4], [1, 1, 1], bootstrap_iterations=100)
        rejected = kg.bootstrap_mrr_improvement([1, 1, 1], [4, 4, 4], bootstrap_iterations=100)

        self.assertTrue(approved["approved"])
        self.assertFalse(rejected["approved"])

    def test_adaptive_empirical_pvalue_extends_only_floor_hits(self):
        metric = np.array([1.0, 0.0, 0.0, 0.0])
        members = np.array([0])
        bins = np.zeros(4, dtype=int)
        floor_hit = kg._empirical_pvalue(
            2.0,
            metric,
            members,
            bins,
            permutations=2,
            rng=np.random.default_rng(17),
            adaptive=True,
            max_permutations=10,
            batch_permutations=2,
        )
        resolved = kg._empirical_pvalue(
            -1.0,
            metric,
            members,
            bins,
            permutations=2,
            rng=np.random.default_rng(17),
            adaptive=True,
            max_permutations=10,
            batch_permutations=2,
        )

        self.assertEqual(floor_hit["permutations"], 10)
        self.assertEqual(floor_hit["exceedances"], 0)
        self.assertTrue(floor_hit["pvalue_is_upper_bound"])
        self.assertAlmostEqual(floor_hit["pvalue"], 1 / 11)
        self.assertEqual(resolved["permutations"], 2)
        self.assertEqual(resolved["exceedances"], 2)
        self.assertFalse(resolved["pvalue_is_upper_bound"])

    def test_leave_one_resource_out_is_source_only_and_reports_all_retrieval_modes(self):
        rows = []
        for library in ["disgenet", "jensen_diseases", "gwas_catalog"]:
            for index in range(3):
                rows.append(
                    {
                        "role": "direct_disease",
                        "library": library,
                        "mondo_id": f"MONDO:{index}",
                        "_member_indices": [index],
                    }
                )
        result = kg.validate_leave_one_resource_out(
            pd.DataFrame(rows),
            np.zeros((4, 4)),
            np.zeros((4, 4)),
            np.eye(4),
            minimum_cases=2,
            bootstrap_iterations=50,
        )

        self.assertEqual(result["status"], "complete")
        self.assertEqual(result["case_count"], 9)
        self.assertIn("direct_mrr", result["retrieval_modes"])
        self.assertIn("directed_mrr", result["retrieval_modes"])
        self.assertIn("unsigned_manifold_mrr", result["retrieval_modes"])
        self.assertIn("hybrid_mrr", result["retrieval_modes"])

    def test_human_perturbseqr_manifest_drops_legacy_mouse_metadata(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            kg_root = root / "kg"
            evidence_root = root / "evidence"
            kg_root.mkdir()
            self._write_cache(kg_root)
            public_evidence.write_evidence_cache(
                evidence_root,
                [
                    {
                        "record_id": "human-perturbation",
                        "resource": "perturbseqr",
                        "species": "human",
                        "species_provenance": "library_name",
                        "kind": "gene_perturbation",
                        "up_genes": ["A"],
                        "down_genes": [],
                        "perturbation_gene": "GENEZ",
                        "perturbation_mode": "KO",
                        "library_name": "Perturb Atlas Human",
                        "source_label": "GENEZ KO up",
                    }
                ],
                manifest={
                    "schema_version": 1,
                    "sources": {
                        "perturbseqr": {"source_url": "https://example.test/perturb", "species": ["human"]},
                        "rummageo": {"species": ["human", "mouse"]},
                    },
                    "ortholog_map": {"species": ["human", "mouse"]},
                },
            )
            result = self._build(kg_root, evidence_root)

        rendered = json.dumps(result["manifest"]).casefold()
        self.assertNotIn("mouse", rendered)
        self.assertEqual(
            result["manifest"]["human_perturbseqr"]["human_perturbseqr_source"]["species"],
            ["human"],
        )


if __name__ == "__main__":
    unittest.main()
