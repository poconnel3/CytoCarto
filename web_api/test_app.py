import unittest
import json
from unittest.mock import patch

import pandas as pd

from web_api import app as api_module
from web_api.app import AnalysisRequest, AnalyteInput, build_graph, get_job, network_cytokine_catalog, request_panel_and_profile, submit_job


class CytoCartoApiTests(unittest.TestCase):
    def tearDown(self):
        api_module.configure_job_backend()

    def test_report_range_overrides_selected_panel_without_mutating_source(self):
        request = AnalysisRequest(
            analytes=[
                AnalyteInput(raw_name="Interleukin 6", value=199.8, units="pg/mL", reference_high=2.0),
                AnalyteInput(raw_name="Interleukin 2", value=2.1, units="pg/mL", reference_high=2.1),
            ],
            age=15,
            sex="male",
        )
        panel, profile, warnings = request_panel_and_profile(request)
        by_id = {row["id"]: row for row in panel["analytes"]}
        self.assertEqual(panel["panel_id"], "example_mito_cytokines")
        self.assertEqual(by_id["IL-6"]["reference_upper_pg_ml"], 2.0)
        self.assertEqual(by_id["IL-2"]["reference_upper_pg_ml"], 2.1)
        self.assertEqual(profile["cytokines"]["Interleukin 6"], 199.8)
        self.assertFalse(warnings)

    def test_ng_ml_is_converted_before_interpreter_input(self):
        request = AnalysisRequest(analytes=[AnalyteInput(raw_name="IL-8", value=0.0069, units="ng/mL", reference_high=0.003)])
        panel, profile, _ = request_panel_and_profile(request)
        by_id = {row["id"]: row for row in panel["analytes"]}
        self.assertAlmostEqual(profile["cytokines"]["IL-8"], 6.9)
        self.assertEqual(by_id["IL-8"]["reference_upper_pg_ml"], 3.0)

    @patch("web_api.app.response_genes_for_stimulus")
    @patch("web_api.app.program_network")
    def test_graph_uses_gene_intermediates_and_never_direct_analyte_program_edges(self, network_mock, response_mock):
        program = "P20-TEST"
        response_mock.return_value = [("STAT3", 0.8, program, -0.4)]
        network_mock.return_value = pd.DataFrame([[0.0]], index=[program], columns=[program])
        graph = build_graph(
            {"analytes": [
                {"analyte_id": "sIL-2R alpha", "display_name": "soluble IL-2 receptor", "elevation_score": 8.0, "stimulus": None, "gene_symbols": ["IL2RA"]},
                {"analyte_id": "IL-6", "display_name": "IL-6", "elevation_score": 2.0, "stimulus": "IL-6", "gene_symbols": ["IL6"]},
            ]},
            {
                "program_scores.tsv": [{"program": program, "program_id": 20, "score": 1.0, "immune_cell_type": "monocytes"}],
                "cell_type_enrichment.tsv": [],
                "gene_subnetwork.tsv": [],
            },
        )
        edge_kinds = [edge["kind"] for edge in graph["edges"]]
        self.assertEqual(edge_kinds, ["analyte_gene", "gene_program"])
        self.assertEqual(graph["edges"][1]["weight"], -0.4)
        self.assertFalse(any(edge["source"].startswith("analyte:") and edge["target"].startswith("program:") for edge in graph["edges"]))
        self.assertFalse(any(node["id"] == "analyte:sIL-2R alpha" for node in graph["nodes"]))

    @patch("web_api.app.program_network")
    def test_context_uses_strongest_program_cosine_contribution(self, network_mock):
        programs = ["P0-A", "P1-B"]
        network_mock.return_value = pd.DataFrame(0.0, index=programs, columns=programs)
        graph = build_graph(
            {"analytes": []},
            {
                "program_scores.tsv": [
                    {"program": "P0-A", "program_id": 0, "score": 1.0, "immune_cell_type": "B cells"},
                    {"program": "P1-B", "program_id": 1, "score": 0.8, "immune_cell_type": "T cells"},
                ],
                "cell_type_enrichment.tsv": [
                    {
                        "cell_type": "Monocyte",
                        "score": 1.0,
                        "program_contributions": '{"P0-A":0.1,"P1-B":0.4}',
                    },
                    {"cell_type": "cDC", "score": 0.5},
                ],
                "gene_subnetwork.tsv": [],
            },
        )
        context_edges = [edge for edge in graph["edges"] if edge["kind"] == "cell_context"]
        self.assertEqual(len(context_edges), 1)
        self.assertEqual(context_edges[0]["source"], "program:P1-B")
        self.assertEqual(context_edges[0]["target"], "cell:Monocyte")
        self.assertEqual(context_edges[0]["weight"], 0.4)

    @patch("web_api.app.network_cytokine_catalog")
    def test_network_only_manual_cytokine_requires_and_uses_report_reference(self, catalog_mock):
        catalog_mock.return_value = ({
            "display_name": "4-1BBL",
            "stimulus": "4-1BBL",
            "source_stimulus": "4-1BBL",
            "reference_upper_pg_ml": None,
        },)
        request = AnalysisRequest(analytes=[AnalyteInput(raw_name="4-1BBL", value=12.5, reference_high=3.0)])
        panel, profile, warnings = request_panel_and_profile(request)
        analyte = next(row for row in panel["analytes"] if row["id"] == "4-1BBL")
        self.assertEqual(analyte["stimulus"], "4-1BBL")
        self.assertEqual(analyte["reference_upper_pg_ml"], 3.0)
        self.assertEqual(profile["cytokines"]["4-1BBL"], 12.5)
        self.assertFalse(warnings)

    @patch("web_api.app.gene_response_matrix")
    @patch("web_api.app.panel_configs")
    def test_network_catalog_normalizes_il18_and_excludes_pbs(self, panels_mock, responses_mock):
        panels_mock.return_value = ((None, {"stimulus_renames": {"IL-18Ra": "IL-18"}, "analytes": []}),)
        responses_mock.return_value = pd.DataFrame(columns=["donor|IL-18Ra", "donor|IL-6", "donor|PBS"])
        network_cytokine_catalog.cache_clear()
        try:
            rows = network_cytokine_catalog()
        finally:
            network_cytokine_catalog.cache_clear()
        self.assertEqual([row["display_name"] for row in rows], ["IL-18", "IL-6"])

    def test_async_job_contract_returns_queued_then_complete(self):
        request = AnalysisRequest(analytes=[AnalyteInput(raw_name="IL-6", value=10.0)])
        api_module.configure_job_backend(
            lambda payload: "remote-job-1",
            lambda job_id: {"job_id": job_id, "status": "running"} if job_id == "remote-job-1" else {"job_id": job_id, "status": "failed"},
        )
        self.assertEqual(submit_job(request), {"job_id": "remote-job-1", "status": "queued"})
        pending = get_job("remote-job-1")
        self.assertEqual(pending.status_code, 202)
        self.assertEqual(json.loads(pending.body)["status"], "running")

        api_module.configure_job_backend(
            lambda payload: "remote-job-2",
            lambda job_id: {"job_id": job_id, "status": "complete", "result": {"ok": True}},
        )
        complete = get_job("remote-job-2")
        self.assertEqual(complete["status"], "complete")
        self.assertEqual(complete["result"], {"ok": True})


if __name__ == "__main__":
    unittest.main()
