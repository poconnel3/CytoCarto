import argparse
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from scripts import cytokine_dspin_interpreter as interp


class CytokineDspinInterpreterTest(unittest.TestCase):
    def test_il18ra_is_normalized_to_il18(self):
        panel = interp.load_json(Path("config/cytokine_panels/mayo_cypan.json"))
        self.assertEqual(interp.normalize_stimulus_name("IL-18Ra", panel), "IL-18")

    def test_default_panel_calculates_optional_il5_il8_il13(self):
        panel = interp.load_json(Path("config/cytokine_panels/mayo_cypan.json"))
        profile = interp.PatientProfile(
            age=None,
            sex=None,
            race=None,
            ethnicity=None,
            cytokines={"Interleukin 5": 6.8, "Interleukin 8": 6.9, "Interleukin 13": 3.1},
        )
        canonical, warnings = interp.canonicalize_cytokines(profile, panel)
        self.assertEqual(warnings, [])
        self.assertEqual(set(canonical), {"IL-5", "IL-8", "IL-13"})
        results = interp.build_analyte_results(canonical, panel, {"IL-5", "IL-8", "IL-13"})
        self.assertEqual({row.analyte_id for row in results}, {"IL-5", "IL-8", "IL-13"})
        self.assertTrue(all(row.stimulus_status == "matched" for row in results))
        self.assertTrue(all(row.elevation_score > 0 for row in results))

    def test_smoothing_uses_directed_signed_edges(self):
        seed = pd.Series({"A": 1.0, "B": 0.0, "C": 0.0})
        network = pd.DataFrame(
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [-1.0, 0.0, 0.0],
            ],
            index=["A", "B", "C"],
            columns=["A", "B", "C"],
        )

        scores = interp.smooth_scores(seed, network, restart=0.5)

        self.assertGreater(scores["B"], 0.0)
        self.assertEqual(scores["C"], 0.0)

        signed_scores = interp.smooth_signed_scores(seed, network, restart=0.5)
        self.assertGreater(signed_scores["B"], 0.0)
        self.assertLess(signed_scores["C"], 0.0)

        reverse_only = pd.DataFrame(
            [
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0],
            ],
            index=["A", "B", "C"],
            columns=["A", "B", "C"],
        )
        reverse_scores = interp.smooth_scores(seed, reverse_only, restart=0.5)
        self.assertEqual(reverse_scores["B"], 0.0)

    def test_gene_to_program_regulators_preserve_sign(self):
        program_scores = pd.Series({"P0": 1.0})
        gene_to_program = pd.DataFrame({"P0": [1.0, -1.0]}, index=["POS", "NEG"])

        _, regulator_score, final_score = interp.build_gene_scores(
            program_scores,
            [],
            gene_to_program,
            None,
        )

        self.assertGreater(regulator_score["POS"], regulator_score["NEG"])
        self.assertGreater(final_score["POS"], final_score["NEG"])
        self.assertEqual(regulator_score["NEG"], 0.0)

    def test_exact_cell_type_enrichment_uses_per_cell_program_matrix(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            programs = ["P0-A", "P1-B"]
            np.save(
                root / "program_program_representation.npy",
                np.array([[1.0, 0.0], [0.9, 0.1], [0.0, 1.0], [0.1, 0.9]]),
            )
            pd.DataFrame({"program": programs}).to_csv(root / "program_program_names.tsv", sep="\t", index=False)
            pd.DataFrame({"obs_name": ["C1", "C2", "C3", "C4"]}).to_csv(
                root / "program_obs_names.tsv", sep="\t", index=False
            )
            metadata = root / "cell_metadata.csv"
            pd.DataFrame(
                {
                    "Cell_ID": ["C1", "C2", "C3", "C4"],
                    "cell_type": ["Mono", "Mono", "B cell", "B cell"],
                }
            ).to_csv(metadata, index=False)

            table = interp.exact_cell_type_enrichment(
                pd.Series({"P0-A": 1.0, "P1-B": 0.0}),
                root,
                metadata,
                top_fraction=0.5,
                chunk_size=2,
            )

            top = table.iloc[0]
            self.assertEqual(top["cell_type"], "Mono")
            self.assertEqual(top["top_cells"], 2)
            self.assertEqual(top["mode"], "exact_per_cell")
            contributions = json.loads(top["program_contributions"])
            self.assertEqual(top["top_contributing_program"], "P0-A")
            self.assertGreater(top["top_program_contribution"], 0.9)
            self.assertAlmostEqual(sum(contributions.values()), top["mean_similarity"], places=7)

    def test_exact_cell_type_enrichment_rejects_row_count_mismatch(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            np.save(root / "program_program_representation.npy", np.array([[1.0], [0.0]]))
            pd.DataFrame({"program": ["P0-A"]}).to_csv(root / "program_program_names.tsv", sep="\t", index=False)
            pd.DataFrame({"obs_name": ["C1", "C2"]}).to_csv(root / "program_obs_names.tsv", sep="\t", index=False)
            metadata = root / "cell_metadata.csv"
            pd.DataFrame({"Cell_ID": ["C1"], "cell_type": ["Mono"]}).to_csv(metadata, index=False)

            with self.assertRaisesRegex(ValueError, "fewer rows|row count"):
                interp.exact_cell_type_enrichment(pd.Series({"P0-A": 1.0}), root, metadata, chunk_size=1)

    def test_exact_cell_type_enrichment_rejects_cell_id_mismatch(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            np.save(root / "program_program_representation.npy", np.array([[1.0], [0.0]]))
            pd.DataFrame({"program": ["P0-A"]}).to_csv(root / "program_program_names.tsv", sep="\t", index=False)
            pd.DataFrame({"obs_name": ["C1", "C2"]}).to_csv(root / "program_obs_names.tsv", sep="\t", index=False)
            metadata = root / "cell_metadata.csv"
            pd.DataFrame({"Cell_ID": ["C1", "WRONG"], "cell_type": ["Mono", "B cell"]}).to_csv(
                metadata, index=False
            )

            with self.assertRaisesRegex(ValueError, "Cell_ID order"):
                interp.exact_cell_type_enrichment(pd.Series({"P0-A": 1.0}), root, metadata, chunk_size=1)

    def test_cell_type_enrichment_falls_back_without_metadata(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            annotations = pd.DataFrame(
                {
                    "program_id": [0, 1],
                    "immune_cell_type": ["Mono", "B cell"],
                    "final_annotation": ["P0", "P1"],
                    "primary_annotation": ["P0", "P1"],
                    "biological_response_or_function": ["", ""],
                }
            )

            table, mode, warnings = interp.build_cell_type_enrichment(
                pd.Series({"P0-A": 1.0, "P1-B": 0.2}),
                annotations,
                root,
                metadata_path=root / "missing_cell_metadata.csv",
            )

            self.assertEqual(mode, "program_annotation_inferred")
            self.assertTrue(warnings)
            self.assertEqual(table.iloc[0]["cell_type"], "Mono")

    def test_covariate_weights_prefer_matching_donor(self):
        donors = {
            "Donor1": {
                "age": 30,
                "sex": "Female",
                "race": "Caucasian",
                "ethnicity": "Not Hispanic/Latino",
            },
            "Donor2": {
                "age": 70,
                "sex": "Male",
                "race": "Caucasian",
                "ethnicity": "Not Hispanic/Latino",
            },
        }
        profile = interp.PatientProfile(
            age=32,
            sex="Female",
            race="Caucasian",
            ethnicity="Not Hispanic/Latino",
            cytokines={},
        )
        weights, notes = interp.covariate_similarity_weights(donors, profile)
        self.assertGreater(weights["Donor1"], weights["Donor2"])
        self.assertTrue(any("covariate-weighted" in note for note in notes))

    def test_run_interpreter_writes_outputs_with_il18_match(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "global"
            root.mkdir()
            output = Path(tmp) / "out"

            programs = ["P0-RSAD2,TNFSF10,IFIT2", "P1-CXCL8,IL1B,NAMPT"]
            pd.DataFrame(
                {
                    "Donor1|IL-18Ra": [1.0, 0.2],
                    "Donor2|IL-18Ra": [0.8, 0.3],
                    "Donor1|IFN-alpha1": [0.9, 0.1],
                    "Donor2|IFN-alpha1": [0.7, 0.1],
                },
                index=programs,
            ).to_csv(root / "program_relative_responses.tsv", sep="\t")
            pd.DataFrame(
                [[0.0, 0.2], [0.2, 0.0]],
                index=programs,
                columns=programs,
            ).to_csv(root / "program_network.tsv", sep="\t")
            pd.DataFrame(
                {
                    programs[0]: [1.0, 0.5, 0.0],
                    programs[1]: [0.0, 0.2, 1.0],
                },
                index=["IL18", "IL18R1", "IL1B"],
            ).to_csv(root / "gene_to_program_regulators.tsv", sep="\t")
            pd.DataFrame(
                [[0.0, 0.4, 0.0], [0.4, 0.0, 0.1], [0.0, 0.1, 0.0]],
                index=["IL18", "IL18R1", "IL1B"],
                columns=["IL18", "IL18R1", "IL1B"],
            ).to_csv(root / "gene_network.tsv", sep="\t")
            pd.DataFrame(
                {
                    "module": [0, 1],
                    "node": programs,
                }
            ).to_csv(root / "program_modules.tsv", sep="\t", index=False)
            pd.DataFrame(
                {
                    "module": [0, 0, 1],
                    "node": ["IL18", "IL18R1", "IL1B"],
                }
            ).to_csv(root / "gene_modules.tsv", sep="\t", index=False)
            pd.DataFrame(
                {
                    "program_id": [0, 1],
                    "program_label": ["P0", "P1"],
                    "final_annotation": ["P0-Type I IFN response", "P1-IL1B inflammatory myeloid"],
                    "primary_annotation": ["Type I IFN", "IL1B inflammation"],
                    "immune_cell_type": ["Broad immune response", "Inflammatory monocytes"],
                    "biological_response_or_function": ["interferon antiviral defense", "acute inflammation"],
                    "confidence_tier": ["high", "high"],
                }
            ).to_csv(root / "dspin_gene_program_annotations.csv", index=False)

            args = argparse.Namespace(
                input="",
                cytokine=["IL-18=936", "IFN-alpha=40"],
                age=30,
                sex="Female",
                race=None,
                ethnicity=None,
                dspin_root=str(root),
                panel_config="config/cytokine_panels/mayo_cypan.json",
                evidence_cache=str(Path(tmp) / "evidence_cache"),
                refresh_evidence=False,
                evidence_refresh_limit=5,
                kg_cache=str(Path(tmp) / "kg_cache"),
                refresh_kg=False,
                kg_permutations=20,
                donor_workbook=str(Path(tmp) / "missing.xlsx"),
                output_dir=str(output),
                top_genes=10,
            )
            written = interp.run_interpreter(args)
            self.assertEqual(written, output)
            self.assertTrue((output / "program_scores.tsv").exists())
            self.assertTrue((output / "gene_subnetwork.tsv").exists())
            self.assertTrue((output / "cell_type_enrichment.tsv").exists())
            self.assertTrue((output / "human_disease_hypotheses.tsv").exists())
            self.assertTrue((output / "human_mimic_hypotheses.tsv").exists())
            self.assertTrue((output / "human_mechanistic_perturbations.tsv").exists())
            self.assertTrue((output / "disease_evidence_edges.tsv").exists())
            self.assertTrue((output / "disease_evidence_manifest.json").exists())
            self.assertTrue((output / "disease_network_concordance.tsv").exists())
            self.assertTrue((output / "disease_manifold_modules.tsv").exists())
            self.assertTrue((output / "disease_supporting_terms.tsv").exists())
            self.assertFalse((output / "mouse_model_hypotheses.tsv").exists())
            self.assertFalse((output / "mouse_mechanistic_perturbations.tsv").exists())
            report = json.loads((output / "patient_report.json").read_text())
            self.assertEqual(report["cell_type_enrichment_mode"], "program_annotation_inferred")
            self.assertTrue(report["top_cell_types"])
            self.assertIn("top_human_disease_hypotheses", report)
            self.assertIn("top_human_mimic_hypotheses", report)
            self.assertNotIn("top_mouse_model_hypotheses", report)
            self.assertEqual(report["disease_engine"], "human_kg_multiview")
            il18 = [row for row in report["analytes"] if row["analyte_id"] == "IL-18"][0]
            self.assertEqual(il18["stimulus_status"], "matched")
            top_program = pd.read_csv(output / "program_scores.tsv", sep="\t").iloc[0]
            self.assertEqual(top_program["program_id"], 0)

            args.refresh_evidence = True
            args.allow_external_patient_signature = False
            with self.assertRaisesRegex(ValueError, "allow-external-patient-signature"):
                interp.run_interpreter(args)


if __name__ == "__main__":
    unittest.main()
