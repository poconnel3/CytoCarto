"""Stateless local API adapter for the existing DSPIN cytokine interpreter."""

from __future__ import annotations

import argparse
import concurrent.futures
import io
import json
import os
import re
import tempfile
import time
import uuid
import zipfile
from functools import lru_cache
from pathlib import Path
from typing import Any, Callable, Literal

import numpy as np
import pandas as pd
from fastapi import FastAPI, HTTPException, Response
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse
from pydantic import BaseModel, ConfigDict, Field, field_validator

from scripts import cytokine_dspin_interpreter as interpreter


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DSPIN_ROOT = Path(
    os.environ.get(
        "DSPIN_GLOBAL_ROOT",
        "/Users/patrickoconnell/Documents/Research Projects/10M_PBMC_project/global",
    )
)
DEFAULT_DEMO_OUTPUT = REPO_ROOT / "artifacts/cytokine_dspin_interpreter/melas_human_kg_adaptive_permutations"
ARTIFACT_FILES = (
    "program_scores.tsv",
    "gene_subnetwork.tsv",
    "cell_type_enrichment.tsv",
    "human_disease_hypotheses.tsv",
    "human_mimic_hypotheses.tsv",
    "human_mechanistic_perturbations.tsv",
    "disease_evidence_edges.tsv",
    "disease_network_concordance.tsv",
    "disease_manifold_modules.tsv",
    "disease_supporting_terms.tsv",
    "disease_evidence_manifest.json",
    "patient_report.json",
    "report.md",
)
TABLE_LIMITS = {
    "disease_evidence_edges.tsv": 300,
    "disease_supporting_terms.tsv": 300,
}
BUNDLE_TTL_SECONDS = 15 * 60
BUNDLES: dict[str, tuple[float, bytes]] = {}
JOB_TTL_SECONDS = 30 * 60
LOCAL_JOB_EXECUTOR = concurrent.futures.ThreadPoolExecutor(max_workers=2, thread_name_prefix="cytocarto-job")
LOCAL_JOBS: dict[str, tuple[float, concurrent.futures.Future[dict[str, Any]]]] = {}
JOB_SUBMITTER: Callable[[dict[str, Any]], str] | None = None
JOB_STATUS_PROVIDER: Callable[[str], dict[str, Any]] | None = None


class AnalyteInput(BaseModel):
    model_config = ConfigDict(extra="forbid", str_strip_whitespace=True)

    raw_name: str = Field(min_length=1, max_length=160)
    value: float = Field(ge=0, le=10_000_000)
    qualifier: Literal["<", ">", "="] | None = None
    units: str = Field(default="pg/mL", max_length=32)
    flag: Literal["H", "L"] | None = None
    reference_low: float | None = Field(default=None, ge=0)
    reference_high: float | None = Field(default=None, gt=0)

    @field_validator("units")
    @classmethod
    def normalize_units(cls, value: str) -> str:
        normalized = value.replace(" ", "").casefold()
        if normalized in {"pg/ml", "pg/mL".casefold(), "ng/l"}:
            return "pg/mL"
        if normalized == "ng/ml":
            return "ng/mL"
        raise ValueError("Only pg/mL, ng/L, and ng/mL cytokine units are supported.")


class AnalysisRequest(BaseModel):
    model_config = ConfigDict(extra="forbid", str_strip_whitespace=True)

    analytes: list[AnalyteInput] = Field(min_length=1, max_length=100)
    age: float | None = Field(default=None, ge=0, le=120)
    sex: str | None = Field(default=None, max_length=80)
    race: str | None = Field(default=None, max_length=120)
    ethnicity: str | None = Field(default=None, max_length=120)


def dataframe_rows(path: Path, limit: int | None = None) -> list[dict[str, Any]]:
    if not path.exists() or path.stat().st_size == 0:
        return []
    frame = pd.read_csv(path, sep="\t", nrows=limit, low_memory=False)
    return json.loads(frame.replace({np.nan: None}).to_json(orient="records"))


def canonical_text(value: str) -> str:
    return interpreter.canonical_text(value)


@lru_cache(maxsize=1)
def panel_configs() -> tuple[tuple[Path, dict[str, Any]], ...]:
    configs = []
    configured_dir = os.environ.get("CYTOCARTO_PANEL_CONFIG_DIR")
    config_dirs = [
        Path(configured_dir) if configured_dir else None,
        Path(__file__).resolve().parent / "panel_configs",
        REPO_ROOT / "config/cytokine_panels",
        Path("/root/config/cytokine_panels"),
    ]
    for config_dir in dict.fromkeys(path for path in config_dirs if path is not None):
        candidates = [
            config_dir / "mayo_cypan.json",
            config_dir / "example_mito_cytokines.json",
        ]
        for path in candidates:
            if path.is_file():
                configs.append((path, json.loads(path.read_text())))
        if configs:
            break
    if not configs:
        raise RuntimeError("No cytokine panel configurations are available.")
    return tuple(configs)


def select_panel(analytes: list[AnalyteInput]) -> tuple[Path, dict[str, Any]]:
    incoming = {canonical_text(item.raw_name) for item in analytes}
    ranked: list[tuple[int, str, Path, dict[str, Any]]] = []
    for path, panel in panel_configs():
        aliases = set(interpreter.build_alias_map(panel))
        ranked.append((len(incoming & aliases), panel.get("panel_id", path.stem), path, panel))
    _, _, path, panel = max(ranked, key=lambda row: (row[0], row[1]))
    return path, json.loads(json.dumps(panel))


def pg_ml_value(item: AnalyteInput) -> float:
    return item.value * (1000.0 if item.units == "ng/mL" else 1.0)


def request_panel_and_profile(request: AnalysisRequest) -> tuple[dict[str, Any], dict[str, Any], list[str]]:
    _, panel = select_panel(request.analytes)
    alias_map = interpreter.build_alias_map(panel)
    unknown_items = [item for item in request.analytes if canonical_text(item.raw_name) not in alias_map]
    catalog_by_name = {
        canonical_text(name): row
        for row in (network_cytokine_catalog() if unknown_items else ())
        for name in (str(row["display_name"]), str(row["stimulus"]), str(row["source_stimulus"]))
    }
    for item in request.analytes:
        key = canonical_text(item.raw_name)
        if key in alias_map or key not in catalog_by_name:
            continue
        if item.reference_high is None:
            raise HTTPException(
                status_code=422,
                detail=f"Enter the laboratory upper reference limit for {item.raw_name}.",
            )
        catalog_row = catalog_by_name[key]
        reference = pg_ml_value(AnalyteInput(raw_name=item.raw_name, value=item.reference_high, units=item.units))
        analyte = {
            "id": str(catalog_row["display_name"]),
            "display_name": str(catalog_row["display_name"]),
            "aliases": [str(catalog_row["display_name"]), str(catalog_row["stimulus"]), str(catalog_row["source_stimulus"])],
            "reference_upper_pg_ml": reference,
            "stimulus": str(catalog_row["stimulus"]),
            "gene_symbols": [],
            "note": "Manually selected DSPIN network stimulus; report reference range used for this request.",
        }
        panel["analytes"].append(analyte)
        alias_map = interpreter.build_alias_map(panel)
    warnings: list[str] = []
    cytokines: dict[str, float] = {}
    overrides: dict[str, float] = {}
    for item in request.analytes:
        analyte = alias_map.get(canonical_text(item.raw_name))
        if analyte is None:
            warnings.append(f"Unrecognized cytokine analyte ignored: {item.raw_name}")
            continue
        analyte_id = str(analyte["id"])
        cytokines[item.raw_name] = pg_ml_value(item)
        if item.reference_high is not None:
            overrides[analyte_id] = pg_ml_value(
                AnalyteInput(raw_name=item.raw_name, value=item.reference_high, units=item.units)
            )
    for analyte in panel["analytes"]:
        if analyte["id"] in overrides:
            analyte["reference_upper_pg_ml"] = overrides[analyte["id"]]
            analyte["note"] = (analyte.get("note", "") + " Report reference range used for this request.").strip()
    profile = {
        "age": request.age,
        "sex": request.sex,
        "race": request.race,
        "ethnicity": request.ethnicity,
        "cytokines": cytokines,
    }
    return panel, profile, warnings


@lru_cache(maxsize=1)
def gene_to_program() -> pd.DataFrame:
    path = DEFAULT_DSPIN_ROOT / "gene_to_program_regulators.tsv"
    if not path.exists():
        raise RuntimeError(f"Missing DSPIN gene-to-program matrix: {path}")
    return pd.read_csv(path, sep="\t", index_col=0)


@lru_cache(maxsize=1)
def program_network() -> pd.DataFrame:
    path = DEFAULT_DSPIN_ROOT / "program_network.tsv"
    if not path.exists():
        raise RuntimeError(f"Missing DSPIN program network: {path}")
    return pd.read_csv(path, sep="\t", index_col=0)


@lru_cache(maxsize=1)
def gene_response_matrix() -> pd.DataFrame:
    path = DEFAULT_DSPIN_ROOT / "gene_relative_responses.tsv"
    if not path.exists():
        raise RuntimeError(f"Missing DSPIN gene-response matrix: {path}")
    return pd.read_csv(path, sep="\t", index_col=0).apply(pd.to_numeric, errors="coerce").fillna(0.0)


@lru_cache(maxsize=1)
def network_cytokine_catalog() -> tuple[dict[str, Any], ...]:
    renames: dict[str, str] = {}
    configured_references: dict[str, set[float]] = {}
    for _, panel in panel_configs():
        renames.update({str(source): str(target) for source, target in panel.get("stimulus_renames", {}).items()})
        for analyte in panel.get("analytes", []):
            stimulus = analyte.get("stimulus")
            if not stimulus:
                continue
            normalized = renames.get(str(stimulus), str(stimulus))
            configured_references.setdefault(canonical_text(normalized), set()).add(float(analyte["reference_upper_pg_ml"]))

    source_stimuli = {
        str(column).split("|", 1)[1]
        for column in gene_response_matrix().columns
        if "|" in str(column)
    }
    rows = []
    for source_stimulus in sorted(source_stimuli, key=str.casefold):
        stimulus = renames.get(source_stimulus, source_stimulus)
        if canonical_text(stimulus) == canonical_text("PBS"):
            continue
        references = configured_references.get(canonical_text(stimulus), set())
        rows.append(
            {
                "display_name": stimulus,
                "stimulus": stimulus,
                "source_stimulus": source_stimulus,
                "reference_upper_pg_ml": next(iter(references)) if len(references) == 1 else None,
            }
        )
    return tuple(rows)


def program_gene_members(program: str) -> list[dict[str, Any]]:
    matrix = gene_to_program().apply(pd.to_numeric, errors="coerce").fillna(0.0)
    if program not in matrix.columns:
        raise KeyError(program)
    rows = [
        {
            "gene": str(gene),
            "program_weight": float(weight),
            "direction": "positive" if weight >= 0 else "negative",
        }
        for gene, weight in matrix[program].items()
    ]
    return sorted(rows, key=lambda row: abs(row["program_weight"]), reverse=True)


def excluded_network_analyte(analyte: dict[str, Any]) -> bool:
    text = " ".join(str(analyte.get(key, "")) for key in ("analyte_id", "display_name")).casefold()
    compact = re.sub(r"[^a-z0-9]+", "", text)
    return "sil2ralpha" in compact or "solubleil2receptor" in compact


def program_contribution_map(cell_type_row: dict[str, Any]) -> dict[str, float]:
    value = cell_type_row.get("program_contributions")
    if isinstance(value, str):
        try:
            value = json.loads(value)
        except json.JSONDecodeError:
            return {}
    if not isinstance(value, dict):
        return {}
    contributions = {}
    for program, contribution in value.items():
        try:
            numeric = float(contribution)
        except (TypeError, ValueError):
            continue
        if np.isfinite(numeric):
            contributions[str(program)] = numeric
    return contributions


def response_genes_for_stimulus(stimulus: str, active_programs: set[str], limit: int = 2) -> list[tuple[str, float, str, float]]:
    responses = gene_response_matrix()
    matching_columns = [
        column
        for column in responses.columns
        if "|" in column and canonical_text(column.split("|", 1)[1]) == canonical_text(stimulus)
    ]
    regulator = gene_to_program().apply(pd.to_numeric, errors="coerce").fillna(0.0)
    programs = [program for program in active_programs if program in regulator.columns]
    common = responses.index.intersection(regulator.index)
    if not matching_columns or not programs or common.empty:
        return []
    response = responses.loc[common, matching_columns].mean(axis=1)
    program_weights = regulator.loc[common, programs]
    strongest_program = program_weights.abs().idxmax(axis=1)
    strongest_weight = pd.Series(
        [float(program_weights.at[gene, strongest_program.at[gene]]) for gene in common],
        index=common,
    )
    relevance = response.abs() * strongest_weight.abs()
    return [
        (str(gene), float(response.at[gene]), str(strongest_program.at[gene]), float(strongest_weight.at[gene]))
        for gene in relevance.nlargest(limit).index
    ]


def build_graph(report: dict[str, Any], tables: dict[str, list[dict[str, Any]]]) -> dict[str, list[dict[str, Any]]]:
    programs = tables["program_scores.tsv"][:5]
    cell_types = tables["cell_type_enrichment.tsv"]
    active_programs = {str(row["program"]) for row in programs}
    nodes: list[dict[str, Any]] = []
    edges: list[dict[str, Any]] = []
    seen_nodes: set[str] = set()

    def node(node_id: str, label: str, kind: str, score: float = 0.0) -> None:
        if node_id not in seen_nodes:
            seen_nodes.add(node_id)
            nodes.append({"id": node_id, "label": label, "kind": kind, "score": score})

    elevated_analytes = sorted(
        (
            row
            for row in report.get("analytes", [])
            if float(row.get("elevation_score", 0)) > 0 and not excluded_network_analyte(row)
        ),
        key=lambda row: float(row.get("elevation_score", 0)),
        reverse=True,
    )[:6]
    for analyte in elevated_analytes:
        analyte_id = f"analyte:{analyte['analyte_id']}"
        node(analyte_id, str(analyte["display_name"]), "analyte", float(analyte["elevation_score"]))
        response_genes = response_genes_for_stimulus(str(analyte.get("stimulus")), active_programs) if analyte.get("stimulus") else []
        if not response_genes:
            response_genes = [(str(gene), 1.0, "", 0.0) for gene in analyte.get("gene_symbols", [])[:1]]
        for gene, response_weight, program, program_weight in response_genes:
            gene_id = f"gene:{gene}"
            node(gene_id, gene, "gene", response_weight)
            edges.append({"id": f"{analyte_id}->{gene_id}", "source": analyte_id, "target": gene_id, "kind": "analyte_gene", "weight": response_weight})
            if program:
                edges.append({"id": f"{gene_id}->program:{program}", "source": gene_id, "target": f"program:{program}", "kind": "gene_program", "weight": program_weight})

    for row in programs:
        node(f"program:{row['program']}", f"P{row['program_id']}", "program", float(row["score"]))
    for row in cell_types:
        node(f"cell:{row['cell_type']}", str(row["cell_type"]), "cell", float(row["score"]))

    network = program_network().apply(pd.to_numeric, errors="coerce").fillna(0.0)
    program_edges = []
    for target in active_programs:
        for source in active_programs:
            if target == source or target not in network.index or source not in network.columns:
                continue
            # DSPIN convention: J[target, source], therefore source -> target.
            weight = float(network.at[target, source])
            program_edges.append((abs(weight), source, target, weight))
    seen_program_pairs: set[frozenset[str]] = set()
    selected_program_edges = []
    for edge in sorted(program_edges, reverse=True):
        _, source, target, _ = edge
        pair = frozenset((source, target))
        if pair in seen_program_pairs:
            continue
        seen_program_pairs.add(pair)
        selected_program_edges.append(edge)
        if len(selected_program_edges) == 6:
            break
    for _, source, target, weight in selected_program_edges:
        edges.append({"id": f"program:{source}->program:{target}", "source": f"program:{source}", "target": f"program:{target}", "kind": "program_program", "weight": weight})

    for cell in cell_types:
        contributions = program_contribution_map(cell)
        candidates = [
            (float(contributions.get(str(program["program"]), 0.0)), program)
            for program in programs
        ]
        positive_candidates = [candidate for candidate in candidates if candidate[0] > 0]
        if not positive_candidates:
            continue
        contribution, program = max(positive_candidates, key=lambda candidate: candidate[0])
        edges.append(
            {
                "id": f"program:{program['program']}--cell:{cell['cell_type']}",
                "source": f"program:{program['program']}",
                "target": f"cell:{cell['cell_type']}",
                "kind": "cell_context",
                "weight": contribution,
            }
        )
    return {"nodes": nodes, "edges": edges}


def create_bundle(output_dir: Path) -> str:
    stream = io.BytesIO()
    with zipfile.ZipFile(stream, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for filename in ARTIFACT_FILES:
            path = output_dir / filename
            if path.exists():
                archive.write(path, arcname=filename)
    now = time.time()
    expired = [token for token, (expires_at, _) in BUNDLES.items() if expires_at <= now]
    for token in expired:
        BUNDLES.pop(token, None)
    token = uuid.uuid4().hex
    BUNDLES[token] = (now + BUNDLE_TTL_SECONDS, stream.getvalue())
    return f"/v1/bundles/{token}"


def build_response(output_dir: Path, selected_panel: str | None = None, input_warnings: list[str] | None = None) -> dict[str, Any]:
    if not (output_dir / "patient_report.json").exists():
        raise RuntimeError("Interpreter did not write patient_report.json.")
    report = json.loads((output_dir / "patient_report.json").read_text())
    tables = {
        name: dataframe_rows(output_dir / name, TABLE_LIMITS.get(name))
        for name in ARTIFACT_FILES
        if name.endswith(".tsv")
    }
    if input_warnings:
        report["warnings"] = list(dict.fromkeys([*report.get("warnings", []), *input_warnings]))
    return {
        "report": report,
        "tables": tables,
        "download_url": create_bundle(output_dir),
        "graph": build_graph(report, tables),
        "selected_panel": selected_panel,
    }


def run_request(request: AnalysisRequest) -> dict[str, Any]:
    panel, profile, input_warnings = request_panel_and_profile(request)
    if not profile["cytokines"]:
        raise HTTPException(status_code=422, detail="No recognized cytokine values were supplied.")
    with tempfile.TemporaryDirectory(prefix="cytocarto-") as temp_dir:
        temp_root = Path(temp_dir)
        panel_path = temp_root / "panel.json"
        input_path = temp_root / "patient.json"
        output_path = temp_root / "output"
        panel_path.write_text(json.dumps(panel))
        input_path.write_text(json.dumps(profile))
        args = argparse.Namespace(
            input=str(input_path),
            cytokine=[],
            age=None,
            sex=None,
            race=None,
            ethnicity=None,
            dspin_root=str(DEFAULT_DSPIN_ROOT),
            panel_config=str(panel_path),
            evidence_cache=os.environ.get("DSPIN_EVIDENCE_CACHE", str(interpreter.DEFAULT_EVIDENCE_CACHE)),
            refresh_evidence=False,
            allow_external_patient_signature=False,
            evidence_refresh_limit=25,
            kg_cache=os.environ.get("DSPIN_KG_CACHE", str(interpreter.DEFAULT_KG_CACHE)),
            refresh_kg=False,
            kg_permutations=int(os.environ.get("CYTOCARTO_KG_PERMUTATIONS", "2000")),
            donor_workbook=os.environ.get("DSPIN_DONOR_WORKBOOK", str(interpreter.DEFAULT_DONOR_WORKBOOK)),
            output_dir=str(output_path),
            top_genes=100,
        )
        interpreter.run_interpreter(args)
        return build_response(output_path, panel.get("panel_id"), input_warnings)


def configure_job_backend(
    submitter: Callable[[dict[str, Any]], str] | None = None,
    status_provider: Callable[[str], dict[str, Any]] | None = None,
) -> None:
    """Install a deployment-specific async backend while retaining local execution."""
    global JOB_SUBMITTER, JOB_STATUS_PROVIDER
    JOB_SUBMITTER = submitter
    JOB_STATUS_PROVIDER = status_provider


def prune_local_jobs() -> None:
    cutoff = time.time() - JOB_TTL_SECONDS
    for job_id, (created_at, _) in list(LOCAL_JOBS.items()):
        if created_at < cutoff:
            LOCAL_JOBS.pop(job_id, None)


def local_job_status(job_id: str) -> dict[str, Any] | JSONResponse:
    job = LOCAL_JOBS.get(job_id)
    if job is None:
        raise HTTPException(status_code=404, detail="Analysis job not found or expired.")
    future = job[1]
    if not future.done():
        return JSONResponse(status_code=202, content={"job_id": job_id, "status": "running"})
    try:
        return {"job_id": job_id, "status": "complete", "result": future.result()}
    except HTTPException as error:
        return JSONResponse(status_code=error.status_code, content={"job_id": job_id, "status": "failed", "detail": error.detail})
    except Exception as error:
        return JSONResponse(status_code=500, content={"job_id": job_id, "status": "failed", "detail": f"CytoCarto analysis failed: {error}"})


app = FastAPI(title="CytoCarto API", version="0.1.0")
origins = [
    origin
    for origin in os.environ.get(
        "CYTOCARTO_ALLOWED_ORIGINS",
        "http://localhost:3000,http://127.0.0.1:3000",
    ).split(",")
    if origin
]
app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=False,
    allow_methods=["GET", "POST"],
    allow_headers=["Content-Type", "Authorization"],
)


@app.get("/health")
def health() -> dict[str, Any]:
    try:
        panel_config_count = len(panel_configs())
    except RuntimeError:
        panel_config_count = 0
    return {
        "status": "ok",
        "dspin_root": str(DEFAULT_DSPIN_ROOT),
        "panel_config_dir": os.environ.get("CYTOCARTO_PANEL_CONFIG_DIR", ""),
        "panel_config_count": panel_config_count,
        "demo_available": DEFAULT_DEMO_OUTPUT.exists(),
    }


@app.get("/v1/demo")
def demo() -> dict[str, Any]:
    if not DEFAULT_DEMO_OUTPUT.exists():
        raise HTTPException(status_code=404, detail="The local MELAS demo output is not available.")
    return build_response(DEFAULT_DEMO_OUTPUT, "example_mito_cytokines")


@app.get("/v1/bundles/{token}")
def download_bundle(token: str) -> Response:
    bundle = BUNDLES.pop(token, None)
    if bundle is None or bundle[0] <= time.time():
        raise HTTPException(status_code=404, detail="This export has expired. Re-run the analysis to generate a new bundle.")
    return Response(
        content=bundle[1],
        media_type="application/zip",
        headers={"Content-Disposition": "attachment; filename=cytocarto-report.zip"},
    )


@app.get("/v1/program-genes")
def get_program_genes(program: str) -> dict[str, Any]:
    try:
        return {"program": program, "genes": program_gene_members(program)}
    except KeyError as error:
        raise HTTPException(status_code=404, detail=f"Unknown DSPIN program: {program}") from error


@app.get("/v1/network-cytokines")
def get_network_cytokines() -> dict[str, Any]:
    return {
        "cytokines": [
            {key: value for key, value in row.items() if key != "source_stimulus"}
            for row in network_cytokine_catalog()
        ]
    }


@app.post("/v1/analyze")
def analyze(request: AnalysisRequest) -> dict[str, Any]:
    try:
        return run_request(request)
    except HTTPException:
        raise
    except Exception as error:
        raise HTTPException(status_code=500, detail=f"CytoCarto analysis failed: {error}") from error


@app.post("/v1/jobs")
def submit_job(request: AnalysisRequest) -> dict[str, Any]:
    prune_local_jobs()
    payload = request.model_dump()
    if JOB_SUBMITTER is not None:
        job_id = JOB_SUBMITTER(payload)
    else:
        job_id = uuid.uuid4().hex
        LOCAL_JOBS[job_id] = (time.time(), LOCAL_JOB_EXECUTOR.submit(run_request, request))
    return {"job_id": job_id, "status": "queued"}


@app.get("/v1/jobs/{job_id}", response_model=None)
def get_job(job_id: str) -> dict[str, Any] | JSONResponse:
    if JOB_STATUS_PROVIDER is not None:
        payload = JOB_STATUS_PROVIDER(job_id)
        if payload.get("status") in {"queued", "running"}:
            return JSONResponse(status_code=202, content=payload)
        return payload
    return local_job_status(job_id)
