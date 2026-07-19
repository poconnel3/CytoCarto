import tempfile
import unittest
from pathlib import Path
from unittest import mock

import pandas as pd

from scripts import public_evidence


class PublicEvidenceTest(unittest.TestCase):
    def setUp(self):
        self.query_genes = pd.Series({"A": 1.0, "B": -0.8, "C": 0.25})
        self.query_programs = pd.Series({"P0": 1.0, "P1": -0.25})
        self.gene_network = pd.DataFrame(
            [[0.0, 0.2, 0.0], [0.2, 0.0, 0.0], [0.0, 0.0, 0.0]],
            index=["A", "B", "C"],
            columns=["A", "B", "C"],
        )
        self.gene_to_program = pd.DataFrame(
            {"P0": [1.0, -1.0, 0.0], "P1": [0.0, 0.0, 1.0]},
            index=["A", "B", "C"],
        )
        self.cell_types = pd.DataFrame({"cell_type": ["Monocyte"]})

    @staticmethod
    def disease_record(record_id, species, study_id, disease_name, up_genes, down_genes):
        return {
            "record_id": record_id,
            "resource": "rummageo",
            "species": species,
            "species_provenance": "source_background",
            "kind": "disease_signature",
            "up_genes": up_genes,
            "down_genes": down_genes,
            "study_id": study_id,
            "source_label": record_id,
            "source_url": f"https://example.test/{record_id}",
            "disease_terms": [
                {
                    "id": f"MONDO:{record_id}",
                    "name": disease_name,
                    "family": "Evidence-backed immune conditions",
                }
            ],
        }

    @staticmethod
    def write_cache(cache_root, records, include_orthologs=True):
        public_evidence.write_evidence_cache(
            cache_root,
            records,
            {
                "sources": {
                    "fixture": {
                        "source_url": "https://example.test",
                        "species": ["human", "mouse"],
                        "source_version": "fixture",
                        "retrieved_at": "2026-07-13T00:00:00",
                    }
                }
            },
        )
        if include_orthologs:
            pd.DataFrame(
                {
                    "human_gene": ["A", "B", "C"],
                    "mouse_gene": ["A", "B", "C"],
                    "one_to_one": [1, 1, 1],
                }
            ).to_csv(cache_root / public_evidence.ORTHOLOG_FILE, sep="\t", index=False)

    def score(self, cache_root):
        return public_evidence.build_public_evidence_hypotheses(
            cache_root,
            self.query_genes,
            self.query_programs,
            self.gene_network,
            self.gene_to_program,
            self.cell_types,
        )

    def test_species_are_separate_directional_and_duplicate_studies_are_collapsed(self):
        records = [
            self.disease_record(
                "human-concordant", "human", "GSE1", "Novel Evidence Disorder", ["A"], ["B"]
            ),
            self.disease_record(
                "human-duplicate", "human", "GSE1", "Novel Evidence Disorder", ["A", "C"], ["B"]
            ),
            {
                "record_id": "rummagene-corroboration",
                "resource": "rummagene",
                "species": "human",
                "species_provenance": "source_metadata",
                "kind": "disease_signature",
                "up_genes": ["A"],
                "down_genes": ["B"],
                "study_id": "PMC99",
                "source_label": "Explicitly human published signature",
                "disease_terms": [
                    {
                        "id": "MONDO:human-concordant",
                        "name": "Novel Evidence Disorder",
                        "family": "Evidence-backed immune conditions",
                    }
                ],
            },
            self.disease_record(
                "human-reversed", "human", "GSE2", "Opposed Evidence Disorder", ["B"], ["A"]
            ),
            self.disease_record(
                "mouse-concordant", "mouse", "GSEM1", "Mouse Model Disorder", ["A"], ["B"]
            ),
            {
                "record_id": "rummagene-unverified",
                "resource": "rummagene",
                "species": "human",
                "kind": "disease_signature",
                "up_genes": ["A"],
                "down_genes": ["B"],
                "disease_terms": [{"id": "MONDO:reject", "name": "Rejected", "family": "Rejected"}],
            },
            {
                "record_id": "mixed-species",
                "resource": "rummageo",
                "species": "mixed",
                "kind": "disease_signature",
                "up_genes": ["A"],
                "down_genes": ["B"],
            },
            {
                "record_id": "rat-species",
                "resource": "rummageo",
                "species": "rat",
                "kind": "disease_signature",
                "up_genes": ["A"],
                "down_genes": ["B"],
            },
            {
                "record_id": "perturb-genez",
                "resource": "perturbseqr",
                "species": "human",
                "species_provenance": "library_name",
                "kind": "gene_perturbation",
                "up_genes": ["A"],
                "down_genes": ["B"],
                "perturbation_gene": "GENEZ",
                "library_name": "Perturb Atlas Human",
                "source_label": "GENEZ knockout up",
            },
            {
                "record_id": "chemical-perturbation",
                "resource": "perturbseqr",
                "species": "human",
                "kind": "gene_perturbation",
                "up_genes": ["A"],
                "perturbation_gene": "CHEM1",
                "library_name": "CREEDS chemical perturbations",
            },
        ]
        with tempfile.TemporaryDirectory() as tmp:
            cache_root = Path(tmp)
            self.write_cache(cache_root, records)
            result = self.score(cache_root)

        human = result["human_diseases"]
        mouse = result["mouse_diseases"]
        self.assertEqual(human.iloc[0]["disease_name"], "Novel Evidence Disorder")
        self.assertEqual(human.iloc[0]["independent_study_count"], 2)
        self.assertEqual(human.iloc[0]["resource_count"], 2)
        self.assertNotIn("Mouse Model Disorder", set(human["disease_name"]))
        self.assertEqual(set(mouse["disease_name"]), {"Mouse Model Disorder"})
        self.assertGreater(
            float(human.loc[human["disease_name"] == "Novel Evidence Disorder", "evidence_score"].iloc[0]),
            float(human.loc[human["disease_name"] == "Opposed Evidence Disorder", "evidence_score"].iloc[0]),
        )
        self.assertGreater(
            float(human.loc[human["disease_name"] == "Opposed Evidence Disorder", "counterevidence_score"].iloc[0]),
            0.0,
        )
        self.assertIn("GENEZ", set(result["human_mechanisms"]["perturbation_gene"]))
        self.assertNotIn("GENEZ", set(human["disease_name"]))
        self.assertGreaterEqual(result["manifest"]["rejected_record_count"], 4)

    def test_cache_only_execution_does_not_read_network_and_missing_orthologs_disable_mouse(self):
        records = [
            self.disease_record("human-only", "human", "GSE1", "Cache-only Disorder", ["A"], ["B"])
        ]
        with tempfile.TemporaryDirectory() as tmp:
            cache_root = Path(tmp)
            self.write_cache(cache_root, records, include_orthologs=False)
            with mock.patch("scripts.public_evidence.urllib.request.urlopen", side_effect=AssertionError("network read")):
                result = self.score(cache_root)

        self.assertFalse(result["human_diseases"].empty)
        self.assertTrue(result["mouse_diseases"].empty)
        self.assertTrue(any("ortholog" in warning.lower() for warning in result["warnings"]))

    def test_refresh_outage_preserves_cached_records(self):
        records = [
            self.disease_record("cached", "human", "GSE1", "Cached Disorder", ["A"], ["B"])
        ]
        with tempfile.TemporaryDirectory() as tmp:
            cache_root = Path(tmp)
            self.write_cache(cache_root, records)
            with (
                mock.patch.object(public_evidence, "refresh_mondo_ontology", return_value=([], "MONDO refresh failed")),
                mock.patch.object(public_evidence, "_rummageo_records", side_effect=OSError("outage")),
                mock.patch.object(public_evidence, "_perturbseqr_records", side_effect=OSError("outage")),
            ):
                manifest, warnings = public_evidence.refresh_public_evidence_cache(
                    cache_root,
                    self.query_genes,
                    {"A": "A", "B": "B", "C": "C"},
                )
            cached, _, _ = public_evidence.read_evidence_cache(cache_root)

        self.assertEqual(manifest["record_count"], 1)
        self.assertEqual(cached[0]["record_id"], "cached")
        self.assertTrue(any("retained cached evidence" in warning for warning in warnings))

    def test_mesh_terms_are_accepted_as_disease_concepts(self):
        record, warning = public_evidence.normalize_evidence_record(
            {
                "record_id": "mesh-record",
                "resource": "rummageo",
                "species": "human",
                "up_genes": ["A"],
                "mesh_terms": [{"id": "D012345", "name": "MeSH Evidence Disorder", "family": "MeSH family"}],
            }
        )

        self.assertIsNone(warning)
        self.assertEqual(record["disease_terms"][0]["id"], "MESH:D012345")

    def test_mouse_symbol_case_and_network_translation_are_preserved(self):
        with tempfile.TemporaryDirectory() as tmp:
            cache_root = Path(tmp)
            pd.DataFrame(
                {
                    "human_gene": ["IL18", "IL1B"],
                    "mouse_gene": ["Il18", "Il1b"],
                    "one_to_one": [1, 1],
                }
            ).to_csv(cache_root / public_evidence.ORTHOLOG_FILE, sep="\t", index=False)
            mapping, warnings = public_evidence.load_ortholog_map(cache_root)

        network = pd.DataFrame([[0.0, 1.0], [1.0, 0.0]], index=["IL18", "IL1B"], columns=["IL18", "IL1B"])
        translated = public_evidence.translate_gene_network_to_mouse(network, mapping)
        vector = public_evidence.record_gene_vector({"up_genes": ["IL18"]}, pd.Index(["Il18", "Il1b"]))

        self.assertFalse(warnings)
        self.assertEqual(mapping, {"IL18": "Il18", "IL1B": "Il1b"})
        self.assertEqual(list(translated.index), ["Il18", "Il1b"])
        self.assertEqual(vector["Il18"], 1.0)

    def test_generic_ontology_terms_are_not_disease_hypotheses(self):
        terms = public_evidence.match_disease_terms(
            "Cancer cohort with juvenile idiopathic arthritis",
            [
                {"id": "MONDO:0000001", "name": "disease", "family": "disease"},
                {"id": "MONDO:0004992", "name": "cancer", "family": "disease"},
                {
                    "id": "MONDO:0005515",
                    "name": "juvenile idiopathic arthritis",
                    "family": "arthritis",
                    "synonyms": [],
                },
            ],
        )

        self.assertEqual(terms, [{"id": "MONDO:0005515", "name": "juvenile idiopathic arthritis", "family": "arthritis"}])

    def test_perturbseqr_live_refresh_queries_species_specific_libraries(self):
        search_calls = []

        def graphql(endpoint, query, variables, timeout):
            search_calls.append((variables["genes"], variables["libraries"]))
            return {"geneSetGeneSearch": {"nodes": []}}

        with mock.patch.object(public_evidence, "_graphql_request", side_effect=graphql):
            records, source_info = public_evidence._perturbseqr_records(
                {"human": pd.Series({"IL18": 1.0}), "mouse": pd.Series({"Il18": 1.0})},
                limit=5,
                timeout=1,
        )

        self.assertEqual(records, [])
        self.assertEqual(
            search_calls,
            [(["IL18"], ["Perturb Atlas Human"]), (["Il18"], ["Perturb Atlas Mouse"])],
        )
        self.assertEqual(source_info["query_mode"], "library_restricted_gene_set_search")

    def test_ortholog_refresh_keeps_only_one_to_one_targets(self):
        def orthology(gene, timeout):
            return gene, {"A": "ENSMUSG1", "B": None}.get(gene)

        with tempfile.TemporaryDirectory() as tmp:
            cache_root = Path(tmp)
            with (
                mock.patch.object(public_evidence, "_one_to_one_target_id", side_effect=orthology),
                mock.patch.object(public_evidence, "_lookup_mouse_symbols", return_value={"ENSMUSG1": "A"}),
            ):
                mapping, warnings = public_evidence.refresh_ortholog_map(
                    cache_root,
                    pd.Series({"A": 1.0, "B": -0.8}),
                )
            manifest = public_evidence.read_cache_manifest(cache_root)

        self.assertEqual(mapping, {"A": "A"})
        self.assertFalse(warnings)
        self.assertEqual(manifest["ortholog_map"]["one_to_one_count"], 1)
        self.assertIn("sha256", manifest["ortholog_map"])

    def test_human_perturbseqr_refresh_never_queries_mouse(self):
        calls = []

        def perturb_records(queries, limit, timeout):
            calls.append(queries)
            return [
                {
                    "record_id": "human-only",
                    "resource": "perturbseqr",
                    "species": "human",
                    "species_provenance": "library_name",
                    "kind": "gene_perturbation",
                    "up_genes": ["A"],
                    "down_genes": [],
                    "perturbation_gene": "GENEZ",
                    "library_name": "Perturb Atlas Human",
                }
            ], {"source_url": "https://example.test", "species": ["human"]}

        with tempfile.TemporaryDirectory() as tmp:
            cache_root = Path(tmp)
            with mock.patch.object(public_evidence, "_perturbseqr_records", side_effect=perturb_records):
                manifest, warnings = public_evidence.refresh_human_perturbseqr_cache(
                    cache_root,
                    pd.Series({"A": 1.0, "B": -0.5}),
                    limit=5,
                )
            records, _, _ = public_evidence.read_evidence_cache(cache_root)

        self.assertFalse(warnings)
        self.assertEqual(list(calls[0]), ["human"])
        self.assertEqual(manifest["sources"]["perturbseqr"]["species"], ["human"])
        self.assertEqual([record["species"] for record in records], ["human"])


if __name__ == "__main__":
    unittest.main()
