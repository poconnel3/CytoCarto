"""Species-separated public-signature retrieval and disease-evidence scoring."""

from __future__ import annotations

import hashlib
import json
import math
import re
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


CACHE_SCHEMA_VERSION = 1
RECORDS_FILE = "signatures.jsonl"
MANIFEST_FILE = "manifest.json"
ORTHOLOG_FILE = "human_mouse_orthologs.tsv"
ONTOLOGY_FILE = "disease_ontology.json"

RUMMAGEO_GRAPHQL = "https://rummageo.com/graphql"
RUMMAGENE_GRAPHQL = "https://rummagene.com/graphql"
PERTURBSEQR_GRAPHQL = "https://perturbseqr.maayanlab.cloud/graphql"
MONDO_OBO_URL = "https://purl.obolibrary.org/obo/mondo.obo"
GEO_SOFT_URL = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?targ=self&form=text&view=quick&acc={gse}"
ENSEMBL_HOMOLOGY_URL = (
    "https://rest.ensembl.org/homology/symbol/human/{gene}?target_species=mus_musculus;type=orthologues"
)
ENSEMBL_LOOKUP_URL = "https://rest.ensembl.org/lookup/id"

SOURCE_URLS = {
    "rummageo": RUMMAGEO_GRAPHQL,
    "rummagene": RUMMAGENE_GRAPHQL,
    "perturbseqr": PERTURBSEQR_GRAPHQL,
}

ALLOWED_RESOURCES = {"rummageo", "rummagene", "perturbseqr"}
ALLOWED_SPECIES = {"human", "mouse"}
EXPLICIT_SPECIES_PROVENANCE = {"explicit", "source_metadata", "source_background", "library_name"}
BROAD_DISEASE_TERM_IDS = {"MONDO:0000001", "MONDO:0004992"}
BROAD_DISEASE_LABELS = {
    "disease",
    "cancer",
    "neoplasm",
    "malignancy",
    "tumor",
    "infectious disease",
    "genetic disease",
    "inherited disease",
    "restricted to specific location",
    "locational disease characteristic",
}

DISEASE_COLUMNS = [
    "rank",
    "disease_family",
    "disease_id",
    "disease_name",
    "evidence_score",
    "counterevidence_score",
    "direct_signature_score",
    "network_proximity_score",
    "program_coherence_score",
    "independent_study_count",
    "resource_count",
    "confidence",
    "matched_genes",
    "supporting_sources",
]

MECHANISM_COLUMNS = [
    "rank",
    "species",
    "perturbation_gene",
    "perturbation_mode",
    "phenocopy_score",
    "reversal_score",
    "network_proximity_score",
    "program_coherence_score",
    "independent_signature_count",
    "supporting_sources",
]

EDGE_COLUMNS = [
    "species",
    "role",
    "resource",
    "record_id",
    "study_id",
    "source_label",
    "disease_id",
    "disease_name",
    "disease_family",
    "perturbation_gene",
    "perturbation_mode",
    "directional_gene_similarity",
    "network_proximity_score",
    "program_coherence_score",
    "context_score",
    "evidence_score",
    "counterevidence_score",
    "matched_genes",
    "opposed_genes",
    "source_url",
]


def canonical_gene(value: Any) -> str:
    return str(value or "").strip().upper()


def display_gene(value: Any) -> str:
    return "" if value is None or pd.isna(value) else str(value).strip()


def canonical_species(value: Any) -> str:
    text = str(value or "").strip().casefold()
    if text in {"human", "homo sapiens", "homo_sapiens"}:
        return "human"
    if text in {"mouse", "mus musculus", "mus_musculus"}:
        return "mouse"
    return ""


def empty_table(columns: list[str]) -> pd.DataFrame:
    return pd.DataFrame(columns=columns)


def cache_path(cache_root: Path, filename: str) -> Path:
    return cache_root / filename


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def records_checksum(records: list[dict[str, Any]]) -> str:
    digest = hashlib.sha256()
    for record in sorted(records, key=lambda item: str(item.get("record_id") or item.get("id") or "")):
        digest.update(json.dumps(record, sort_keys=True, separators=(",", ":")).encode())
        digest.update(b"\n")
    return digest.hexdigest()


def default_manifest() -> dict[str, Any]:
    return {
        "schema_version": CACHE_SCHEMA_VERSION,
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "sources": {},
    }


def read_cache_manifest(cache_root: Path) -> dict[str, Any]:
    path = cache_path(cache_root, MANIFEST_FILE)
    if not path.exists():
        return default_manifest()
    try:
        payload = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError):
        return default_manifest()
    return payload if isinstance(payload, dict) else default_manifest()


def write_evidence_cache(
    cache_root: Path,
    records: list[dict[str, Any]],
    manifest: dict[str, Any] | None = None,
) -> dict[str, Any]:
    cache_root.mkdir(parents=True, exist_ok=True)
    record_path = cache_path(cache_root, RECORDS_FILE)
    with record_path.open("w") as handle:
        for record in records:
            handle.write(json.dumps(record, sort_keys=True) + "\n")

    payload = dict(manifest or default_manifest())
    sources = dict(payload.get("sources", {}))
    for resource in ALLOWED_RESOURCES:
        source_records = [record for record in records if str(record.get("resource", "")).casefold() == resource]
        if not source_records and resource not in sources:
            continue
        source = dict(sources.get(resource, {}))
        source.setdefault("source_url", SOURCE_URLS[resource])
        source.setdefault("source_version", "cached")
        source.setdefault("retrieved_at", "unknown")
        source["species"] = sorted(
            {
                species
                for record in source_records
                if (species := canonical_species(record.get("species"))) in ALLOWED_SPECIES
            }
        )
        source["record_count"] = len(source_records)
        source["records_sha256"] = records_checksum(source_records)
        sources[resource] = source
    payload["sources"] = sources
    payload["schema_version"] = CACHE_SCHEMA_VERSION
    payload["updated_at"] = datetime.now().isoformat(timespec="seconds")
    payload["record_count"] = len(records)
    payload["records_sha256"] = sha256_file(record_path)
    cache_path(cache_root, MANIFEST_FILE).write_text(json.dumps(payload, indent=2) + "\n")
    return payload


def read_evidence_cache(cache_root: Path) -> tuple[list[dict[str, Any]], dict[str, Any], list[str]]:
    manifest = read_cache_manifest(cache_root)
    path = cache_path(cache_root, RECORDS_FILE)
    if not path.exists():
        return [], manifest, [
            f"Public evidence cache is missing at {path}; run with --refresh-evidence or add normalized records."
        ]

    records: list[dict[str, Any]] = []
    warnings: list[str] = []
    with path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            text = line.strip()
            if not text:
                continue
            try:
                record = json.loads(text)
            except json.JSONDecodeError:
                warnings.append(f"Ignoring invalid evidence-cache JSON at {path}:{line_number}.")
                continue
            if isinstance(record, dict):
                records.append(record)
            else:
                warnings.append(f"Ignoring non-object evidence-cache entry at {path}:{line_number}.")
    return records, manifest, warnings


def _as_gene_list(value: Any) -> list[str]:
    if not isinstance(value, list):
        return []
    return sorted({gene for raw in value if (gene := canonical_gene(raw))})


def is_specific_disease_term(term_id: str, term_name: str) -> bool:
    return term_id not in BROAD_DISEASE_TERM_IDS and term_name.strip().casefold() not in BROAD_DISEASE_LABELS


def normalize_evidence_record(record: dict[str, Any]) -> tuple[dict[str, Any] | None, str | None]:
    resource = str(record.get("resource", "")).strip().casefold()
    species = canonical_species(record.get("species"))
    if resource not in ALLOWED_RESOURCES:
        return None, f"unsupported resource {resource or '<missing>'}"
    if species not in ALLOWED_SPECIES:
        return None, "missing explicit human/mouse species"

    normalized = dict(record)
    normalized["resource"] = resource
    normalized["species"] = species
    normalized["record_id"] = str(record.get("record_id") or record.get("id") or "").strip()
    if not normalized["record_id"]:
        return None, "missing record_id"
    normalized["up_genes"] = _as_gene_list(record.get("up_genes"))
    normalized["down_genes"] = _as_gene_list(record.get("down_genes"))
    if not normalized["up_genes"] and not normalized["down_genes"]:
        return None, "missing signed gene sets"

    normalized["kind"] = str(record.get("kind", "disease_signature")).strip().casefold()
    normalized["study_id"] = str(record.get("study_id", "")).strip()
    normalized["source_label"] = str(
        record.get("source_label") or record.get("condition") or record.get("title") or normalized["record_id"]
    ).strip()
    normalized["source_url"] = str(record.get("source_url", "")).strip()
    normalized["cell_type"] = str(record.get("cell_type", "")).strip()
    normalized["perturbation_gene"] = canonical_gene(record.get("perturbation_gene"))
    normalized["perturbation_mode"] = str(record.get("perturbation_mode", "")).strip()
    normalized["library_name"] = str(record.get("library_name", "")).strip()
    normalized["species_provenance"] = str(record.get("species_provenance", "")).strip().casefold()
    try:
        normalized["source_confidence"] = min(max(float(record.get("source_confidence", 1.0)), 0.0), 1.0)
    except (TypeError, ValueError):
        normalized["source_confidence"] = 1.0

    if resource == "rummagene" and normalized["species_provenance"] not in EXPLICIT_SPECIES_PROVENANCE:
        return None, "Rummagene record lacks explicit source species provenance"
    if resource == "perturbseqr":
        library = normalized["library_name"].casefold()
        if normalized["kind"] != "gene_perturbation" or not normalized["perturbation_gene"]:
            return None, "Perturb-Seqr record is not a gene perturbation"
        if "chem" in library or "drug" in library or "cmap" in library:
            return None, "Perturb-Seqr chemical signature excluded"

    source_terms = record.get("disease_terms", [])
    if not isinstance(source_terms, list):
        source_terms = []
    terms = list(source_terms)
    mesh_terms = record.get("mesh_terms", [])
    if isinstance(mesh_terms, list):
        terms = [*terms, *mesh_terms]
    normalized_terms = []
    for term_index, term in enumerate(terms):
        if not isinstance(term, dict):
            continue
        term_id = str(term.get("id", "")).strip()
        term_name = str(term.get("name", "")).strip()
        if term_index >= len(source_terms) and term_id and not term_id.casefold().startswith("mesh:"):
            term_id = f"MESH:{term_id}"
        if term_id and term_name and is_specific_disease_term(term_id, term_name):
            normalized_terms.append(
                {
                    "id": term_id,
                    "name": term_name,
                    "family": str(term.get("family", "")).strip() or "Unclassified evidence-backed conditions",
                }
            )
    normalized["disease_terms"] = normalized_terms
    return normalized, None


def load_ortholog_map(cache_root: Path) -> tuple[dict[str, str], list[str]]:
    path = cache_path(cache_root, ORTHOLOG_FILE)
    if not path.exists():
        return {}, [
            f"Mouse model evidence unavailable because one-to-one ortholog map is missing at {path}."
        ]
    frame = pd.read_csv(path, sep="\t")
    required = {"human_gene", "mouse_gene", "one_to_one"}
    missing = required - set(frame.columns)
    if missing:
        return {}, [f"Ignoring ortholog map missing required columns: {', '.join(sorted(missing))}."]
    allowed = frame["one_to_one"].astype(str).str.casefold().isin({"1", "true", "yes"})
    candidate_pairs = [
        (canonical_gene(row.human_gene), display_gene(row.mouse_gene))
        for row in frame.loc[allowed, ["human_gene", "mouse_gene"]].itertuples(index=False)
        if canonical_gene(row.human_gene) and display_gene(row.mouse_gene)
    ]
    mouse_counts = pd.Series([canonical_gene(mouse_gene) for _, mouse_gene in candidate_pairs]).value_counts()
    mapping = {
        human_gene: mouse_gene
        for human_gene, mouse_gene in candidate_pairs
        if mouse_counts[canonical_gene(mouse_gene)] == 1
    }
    if not mapping:
        return {}, ["No high-confidence one-to-one human-mouse orthologs were available."]
    return mapping, []


def map_human_query_to_mouse(query: pd.Series, orthologs: dict[str, str]) -> pd.Series:
    values: dict[str, float] = {}
    for gene, value in query.items():
        mouse_gene = orthologs.get(canonical_gene(gene))
        if not mouse_gene:
            continue
        numeric = float(value)
        if abs(numeric) > abs(values.get(mouse_gene, 0.0)):
            values[mouse_gene] = numeric
    return normalize_signed_series(pd.Series(values, dtype=float))


def _one_to_one_target_id(gene: str, timeout: int) -> tuple[str, str | None]:
    url = ENSEMBL_HOMOLOGY_URL.format(gene=urllib.parse.quote(gene, safe=""))
    payload = _json_request(url, timeout=timeout)
    data = payload.get("data", []) if isinstance(payload, dict) else []
    homologies = data[0].get("homologies", []) if data and isinstance(data[0], dict) else []
    target_ids = {
        str(item.get("target", {}).get("id", "")).strip()
        for item in homologies
        if isinstance(item, dict)
        and item.get("type") == "ortholog_one2one"
        and item.get("target", {}).get("species") == "mus_musculus"
    }
    target_ids.discard("")
    return gene, next(iter(target_ids)) if len(target_ids) == 1 else None


def _lookup_mouse_symbols(target_ids: list[str], timeout: int) -> dict[str, str]:
    symbols: dict[str, str] = {}
    for start in range(0, len(target_ids), 1000):
        payload = _json_request(ENSEMBL_LOOKUP_URL, {"ids": target_ids[start : start + 1000]}, timeout=timeout)
        if not isinstance(payload, dict):
            continue
        for target_id, item in payload.items():
            if not isinstance(item, dict) or item.get("species") != "mus_musculus":
                continue
            symbol = display_gene(item.get("display_name"))
            if symbol:
                symbols[str(target_id)] = symbol
    return symbols


def refresh_ortholog_map(
    cache_root: Path,
    signed_gene_scores: pd.Series,
    limit: int = 75,
    timeout: int = 45,
) -> tuple[dict[str, str], list[str]]:
    existing, warnings = load_ortholog_map(cache_root)
    warnings = [warning for warning in warnings if "is missing" not in warning]
    stale_generated_genes: set[str] = set()
    ortholog_path = cache_path(cache_root, ORTHOLOG_FILE)
    if ortholog_path.exists():
        raw_orthologs = pd.read_csv(ortholog_path, sep="\t")
        if {"human_gene", "source_version"}.issubset(raw_orthologs.columns) and "symbol_case_preserved" not in raw_orthologs:
            stale_generated_genes = {
                canonical_gene(gene)
                for gene in raw_orthologs.loc[
                    raw_orthologs["source_version"].astype(str).eq("Ensembl REST"), "human_gene"
                ]
                if canonical_gene(gene)
            }
    if stale_generated_genes:
        existing = {gene: mouse_gene for gene, mouse_gene in existing.items() if gene not in stale_generated_genes}
    ranked = normalize_signed_series(signed_gene_scores).abs().sort_values(ascending=False)
    query_genes = [canonical_gene(gene) for gene in ranked.index if float(ranked.loc[gene]) > 0][:limit]
    pending = [gene for gene in query_genes if gene and (gene not in existing or gene in stale_generated_genes)]
    target_ids: dict[str, str] = {}
    failed = 0
    if pending:
        with ThreadPoolExecutor(max_workers=6) as executor:
            futures = {executor.submit(_one_to_one_target_id, gene, timeout): gene for gene in pending}
            for future in as_completed(futures):
                try:
                    gene, target_id = future.result()
                except (urllib.error.URLError, TimeoutError, OSError, RuntimeError, json.JSONDecodeError):
                    failed += 1
                    continue
                if target_id:
                    target_ids[gene] = target_id
    try:
        mouse_symbols = _lookup_mouse_symbols(sorted(set(target_ids.values())), timeout) if target_ids else {}
    except (urllib.error.URLError, TimeoutError, OSError, RuntimeError, json.JSONDecodeError):
        mouse_symbols = {}
        failed += len(target_ids)

    additions = {
        gene: mouse_symbols[target_id]
        for gene, target_id in target_ids.items()
        if target_id in mouse_symbols
    }
    candidate_mapping = {**existing, **additions}
    target_counts = pd.Series([canonical_gene(mouse_gene) for mouse_gene in candidate_mapping.values()]).value_counts()
    mapping = {
        human_gene: mouse_gene
        for human_gene, mouse_gene in candidate_mapping.items()
        if target_counts[canonical_gene(mouse_gene)] == 1
    }
    if mapping:
        cache_root.mkdir(parents=True, exist_ok=True)
        retrieved_at = datetime.now().isoformat(timespec="seconds")
        pd.DataFrame(
            [
                {
                    "human_gene": gene,
                    "mouse_gene": mouse_gene,
                    "one_to_one": True,
                    "source_url": ENSEMBL_HOMOLOGY_URL,
                    "source_version": "Ensembl REST",
                    "symbol_case_preserved": True,
                    "retrieved_at": retrieved_at,
                }
                for gene, mouse_gene in sorted(mapping.items())
            ]
        ).to_csv(cache_path(cache_root, ORTHOLOG_FILE), sep="\t", index=False)
        ortholog_path = cache_path(cache_root, ORTHOLOG_FILE)
        manifest = read_cache_manifest(cache_root)
        manifest["ortholog_map"] = {
            "source_url": ENSEMBL_HOMOLOGY_URL,
            "species": ["human", "mouse"],
            "source_version": "Ensembl REST",
            "retrieved_at": retrieved_at,
            "input_gene_count": len(query_genes),
            "one_to_one_count": len(mapping),
            "sha256": sha256_file(ortholog_path),
        }
        cache_path(cache_root, MANIFEST_FILE).write_text(json.dumps(manifest, indent=2) + "\n")
    if failed:
        warnings.append(f"Ensembl one-to-one ortholog refresh failed for {failed} query gene(s); cached mappings were retained.")
    if not mapping:
        warnings.append("No high-confidence one-to-one human-mouse orthologs were available for the signed query.")
    return mapping, warnings


def normalize_signed_series(series: pd.Series) -> pd.Series:
    values = series.astype(float).replace([np.inf, -np.inf], np.nan).fillna(0.0)
    maximum = float(values.abs().max()) if not values.empty else 0.0
    return values * 0.0 if maximum <= 0 else values / maximum


def smooth_signed_scores(seed: pd.Series, network: pd.DataFrame | None, restart: float = 0.8) -> pd.Series:
    seed = seed.astype(float).fillna(0.0)
    if network is None or network.empty:
        return normalize_signed_series(seed)
    common = [gene for gene in seed.index if gene in network.index and gene in network.columns]
    if not common:
        return normalize_signed_series(seed)
    matrix = network.loc[common, common].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    # pandas may expose a read-only NumPy view under recent releases.
    values = matrix.to_numpy(dtype=float, copy=True)
    np.fill_diagonal(values, 0.0)
    matrix = pd.DataFrame(values, index=matrix.index, columns=matrix.columns)
    row_sum = matrix.abs().sum(axis=1).replace(0, np.nan)
    propagated = matrix.div(row_sum, axis=0).fillna(0.0).dot(seed.loc[common])
    result = seed.copy()
    result.loc[common] = restart * seed.loc[common] + (1.0 - restart) * propagated
    return normalize_signed_series(result)


def signed_cosine(left: pd.Series, right: pd.Series) -> float:
    common = left.index.intersection(right.index)
    if common.empty:
        return 0.0
    left_values = left.loc[common].astype(float).to_numpy()
    right_values = right.loc[common].astype(float).to_numpy()
    denominator = float(np.linalg.norm(left_values) * np.linalg.norm(right_values))
    return 0.0 if denominator <= 0 else float(np.dot(left_values, right_values) / denominator)


def record_gene_vector(record: dict[str, Any], index: pd.Index) -> pd.Series:
    vector = pd.Series(0.0, index=index.astype(str), dtype=float)
    index_by_gene = {canonical_gene(gene): gene for gene in vector.index}
    for gene in record.get("up_genes", []):
        if (target := index_by_gene.get(canonical_gene(gene))) is not None:
            vector.loc[target] = 1.0
    for gene in record.get("down_genes", []):
        if (target := index_by_gene.get(canonical_gene(gene))) is not None:
            vector.loc[target] = -1.0
    return vector


def translate_gene_network_to_mouse(
    gene_network: pd.DataFrame | None,
    orthologs: dict[str, str],
) -> pd.DataFrame | None:
    if gene_network is None or gene_network.empty or not orthologs:
        return None
    common = [
        str(gene)
        for gene in gene_network.index
        if str(gene) in gene_network.columns and canonical_gene(gene) in orthologs
    ]
    if not common:
        return None
    translated = gene_network.loc[common, common].copy()
    translated.index = [orthologs[canonical_gene(gene)] for gene in common]
    translated.columns = translated.index
    return translated


def translate_gene_to_program_to_mouse(
    gene_to_program: pd.DataFrame | None,
    orthologs: dict[str, str],
) -> pd.DataFrame | None:
    if gene_to_program is None or gene_to_program.empty or not orthologs:
        return None
    common = [str(gene) for gene in gene_to_program.index if canonical_gene(gene) in orthologs]
    if not common:
        return None
    translated = gene_to_program.loc[common].copy()
    translated.index = [orthologs[canonical_gene(gene)] for gene in common]
    return translated


def directional_gene_similarity(query: pd.Series, source: pd.Series) -> tuple[float, list[str], list[str]]:
    common = query.index.intersection(source.index)
    matched: list[str] = []
    opposed: list[str] = []
    match_weight = 0.0
    oppose_weight = 0.0
    query_weight = float(query.abs().sum())
    for gene in common:
        query_value = float(query.loc[gene])
        source_value = float(source.loc[gene])
        if query_value == 0 or source_value == 0:
            continue
        if math.copysign(1.0, query_value) == math.copysign(1.0, source_value):
            matched.append(str(gene))
            match_weight += abs(query_value)
        else:
            opposed.append(str(gene))
            oppose_weight += abs(query_value)
    overlap = match_weight + oppose_weight
    if overlap <= 0 or query_weight <= 0:
        return 0.0, matched, opposed
    sign_agreement = (match_weight - oppose_weight) / overlap
    coverage = min(1.0, overlap / query_weight)
    return float(sign_agreement * math.sqrt(coverage)), matched, opposed


def project_signature_to_programs(source_genes: pd.Series, gene_to_program: pd.DataFrame | None) -> pd.Series:
    if gene_to_program is None or gene_to_program.empty:
        return pd.Series(dtype=float)
    common = source_genes.index.intersection(gene_to_program.index.astype(str))
    if common.empty:
        return pd.Series(0.0, index=gene_to_program.columns.astype(str), dtype=float)
    matrix = gene_to_program.copy()
    matrix.index = matrix.index.astype(str)
    matrix = matrix.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    return normalize_signed_series(matrix.loc[common].T.dot(source_genes.loc[common]))


def context_score(record: dict[str, Any], patient_cell_types: set[str]) -> float:
    source_cell_type = str(record.get("cell_type", "")).strip().casefold()
    if not source_cell_type:
        return 0.0
    return 1.0 if source_cell_type in patient_cell_types else -1.0


def score_record(
    record: dict[str, Any],
    query_genes: pd.Series,
    query_programs: pd.Series,
    gene_network: pd.DataFrame | None,
    gene_to_program: pd.DataFrame | None,
    patient_cell_types: set[str],
) -> dict[str, Any]:
    source_genes = record_gene_vector(record, query_genes.index)
    direct_score, matched, opposed = directional_gene_similarity(query_genes, source_genes)
    propagated = smooth_signed_scores(source_genes, gene_network)
    network_score = signed_cosine(query_genes, propagated)
    source_programs = project_signature_to_programs(source_genes, gene_to_program)
    program_score = signed_cosine(query_programs, source_programs)
    source_context = context_score(record, patient_cell_types)
    signed_score = 0.60 * direct_score + 0.25 * network_score + 0.10 * program_score + 0.05 * source_context
    signed_score *= float(record.get("source_confidence", 1.0))
    return {
        "directional_gene_similarity": direct_score,
        "network_proximity_score": network_score,
        "program_coherence_score": program_score,
        "context_score": source_context,
        "evidence_score": max(0.0, signed_score),
        "counterevidence_score": max(0.0, -signed_score),
        "matched_genes": matched,
        "opposed_genes": opposed,
    }


def disease_terms(record: dict[str, Any]) -> list[dict[str, str]]:
    return [term for term in record.get("disease_terms", []) if term.get("id") and term.get("name")]


def _edge_row(
    record: dict[str, Any],
    score: dict[str, Any],
    role: str,
    term: dict[str, str] | None = None,
) -> dict[str, Any]:
    return {
        "species": record["species"],
        "role": role,
        "resource": record["resource"],
        "record_id": record["record_id"],
        "study_id": record.get("study_id", ""),
        "source_label": record.get("source_label", ""),
        "disease_id": term.get("id", "") if term else "",
        "disease_name": term.get("name", "") if term else "",
        "disease_family": term.get("family", "") if term else "",
        "perturbation_gene": record.get("perturbation_gene", ""),
        "perturbation_mode": record.get("perturbation_mode", ""),
        "directional_gene_similarity": score["directional_gene_similarity"],
        "network_proximity_score": score["network_proximity_score"],
        "program_coherence_score": score["program_coherence_score"],
        "context_score": score["context_score"],
        "evidence_score": score["evidence_score"],
        "counterevidence_score": score["counterevidence_score"],
        "matched_genes": ";".join(score["matched_genes"]),
        "opposed_genes": ";".join(score["opposed_genes"]),
        "source_url": record.get("source_url", ""),
    }


def aggregate_diseases(edges: pd.DataFrame, species: str) -> pd.DataFrame:
    relevant = edges.loc[
        (edges["species"] == species) & (edges["role"].isin(["disease_evidence", "disease_counterevidence"]))
    ].copy()
    if relevant.empty:
        return empty_table(DISEASE_COLUMNS)

    rows: list[dict[str, Any]] = []
    group_columns = ["disease_family", "disease_id", "disease_name"]
    for keys, disease_edges in relevant.groupby(group_columns, dropna=False):
        disease_family, disease_id, disease_name = keys
        study_key = disease_edges["study_id"].astype(str).where(
            disease_edges["study_id"].astype(str).str.len() > 0,
            disease_edges["record_id"].astype(str),
        )
        collapsed = disease_edges.assign(_study_key=disease_edges["resource"].astype(str) + ":" + study_key)
        selected = []
        for _, study_edges in collapsed.groupby("_study_key", sort=False):
            strongest = study_edges.assign(
                _strength=study_edges[["evidence_score", "counterevidence_score"]].max(axis=1)
            ).sort_values("_strength", ascending=False).iloc[0]
            selected.append(strongest)
        independent = pd.DataFrame(selected)
        positive = independent["evidence_score"].clip(lower=0.0, upper=1.0).to_numpy(dtype=float)
        negative = independent["counterevidence_score"].clip(lower=0.0, upper=1.0).to_numpy(dtype=float)
        evidence_score = float(1.0 - np.prod(1.0 - positive)) if positive.size else 0.0
        counter_score = float(1.0 - np.prod(1.0 - negative)) if negative.size else 0.0
        resource_count = int(independent.loc[independent["evidence_score"] > 0, "resource"].nunique())
        independent_count = int((independent["evidence_score"] > 0).sum())
        if independent_count >= 2 and resource_count >= 2:
            confidence = "multi-source"
        elif independent_count >= 2:
            confidence = "multi-study"
        elif independent_count == 1:
            confidence = "single-study"
        else:
            confidence = "counterevidence-only"
        positive_rows = independent.loc[independent["evidence_score"] > 0]
        sources = positive_rows["source_label"].astype(str).drop_duplicates().head(5).tolist()
        genes = []
        for value in positive_rows["matched_genes"].astype(str):
            genes.extend(gene for gene in value.split(";") if gene)
        rows.append(
            {
                "disease_family": disease_family or "Unclassified evidence-backed conditions",
                "disease_id": disease_id,
                "disease_name": disease_name,
                "evidence_score": evidence_score,
                "counterevidence_score": counter_score,
                "direct_signature_score": float(positive_rows["directional_gene_similarity"].mean())
                if not positive_rows.empty
                else 0.0,
                "network_proximity_score": float(positive_rows["network_proximity_score"].mean())
                if not positive_rows.empty
                else 0.0,
                "program_coherence_score": float(positive_rows["program_coherence_score"].mean())
                if not positive_rows.empty
                else 0.0,
                "independent_study_count": independent_count,
                "resource_count": resource_count,
                "confidence": confidence,
                "matched_genes": ";".join(sorted(set(genes))[:25]),
                "supporting_sources": "; ".join(sources),
            }
        )
    table = pd.DataFrame(rows, columns=DISEASE_COLUMNS[1:])
    table = table.sort_values(
        ["evidence_score", "counterevidence_score", "independent_study_count"],
        ascending=[False, True, False],
    ).reset_index(drop=True)
    table.insert(0, "rank", np.arange(1, len(table) + 1))
    return table[DISEASE_COLUMNS]


def aggregate_mechanisms(edges: pd.DataFrame, species: str) -> pd.DataFrame:
    relevant = edges.loc[
        (edges["species"] == species)
        & (edges["resource"] == "perturbseqr")
        & edges["perturbation_gene"].astype(str).str.len().gt(0)
    ].copy()
    if relevant.empty:
        return empty_table(MECHANISM_COLUMNS)

    rows: list[dict[str, Any]] = []
    for keys, group in relevant.groupby(["perturbation_gene", "perturbation_mode"], dropna=False):
        gene, mode = keys
        signature_key = group["study_id"].astype(str).where(
            group["study_id"].astype(str).str.len() > 0, group["record_id"].astype(str)
        )
        group = group.assign(_signature_key=signature_key)
        selected = []
        for _, values in group.groupby("_signature_key", sort=False):
            selected.append(
                values.assign(_strength=values[["evidence_score", "counterevidence_score"]].max(axis=1))
                .sort_values("_strength", ascending=False)
                .iloc[0]
            )
        independent = pd.DataFrame(selected)
        rows.append(
            {
                "species": species,
                "perturbation_gene": gene,
                "perturbation_mode": mode,
                "phenocopy_score": float(independent["evidence_score"].max()),
                "reversal_score": float(independent["counterevidence_score"].max()),
                "network_proximity_score": float(independent["network_proximity_score"].mean()),
                "program_coherence_score": float(independent["program_coherence_score"].mean()),
                "independent_signature_count": int(len(independent)),
                "supporting_sources": "; ".join(independent["source_label"].astype(str).drop_duplicates().head(5)),
            }
        )
    table = pd.DataFrame(rows, columns=MECHANISM_COLUMNS[1:])
    table = table.sort_values(["phenocopy_score", "reversal_score"], ascending=[False, False]).reset_index(drop=True)
    table.insert(0, "rank", np.arange(1, len(table) + 1))
    return table[MECHANISM_COLUMNS]


def build_public_evidence_hypotheses(
    cache_root: Path,
    signed_gene_scores: pd.Series,
    signed_program_scores: pd.Series,
    gene_network: pd.DataFrame | None,
    gene_to_program: pd.DataFrame | None,
    cell_type_table: pd.DataFrame,
) -> dict[str, Any]:
    records, cache_manifest, warnings = read_evidence_cache(cache_root)
    orthologs, ortholog_warnings = load_ortholog_map(cache_root)
    warnings.extend(ortholog_warnings)
    human_query = normalize_signed_series(signed_gene_scores.copy())
    mouse_query = map_human_query_to_mouse(human_query, orthologs) if orthologs else pd.Series(dtype=float)
    queries = {"human": human_query, "mouse": mouse_query}
    networks = {
        "human": gene_network,
        "mouse": translate_gene_network_to_mouse(gene_network, orthologs),
    }
    regulator_matrices = {
        "human": gene_to_program,
        "mouse": translate_gene_to_program_to_mouse(gene_to_program, orthologs),
    }
    patient_cell_types = {
        str(value).strip().casefold()
        for value in cell_type_table.get("cell_type", pd.Series(dtype=str)).head(10)
        if str(value).strip()
    }

    edges: list[dict[str, Any]] = []
    rejected: list[str] = []
    accepted = 0
    for raw_record in records:
        record, reason = normalize_evidence_record(raw_record)
        if record is None:
            rejected.append(reason or "invalid record")
            continue
        query = queries[record["species"]]
        if query.empty or float(query.abs().sum()) <= 0:
            continue
        accepted += 1
        score = score_record(
            record,
            query,
            signed_program_scores,
            networks[record["species"]],
            regulator_matrices[record["species"]],
            patient_cell_types,
        )
        if record["resource"] == "perturbseqr":
            role = "mechanistic_phenocopy" if score["evidence_score"] >= score["counterevidence_score"] else "mechanistic_reversal"
            edges.append(_edge_row(record, score, role))
            continue
        terms = disease_terms(record)
        if not terms:
            edges.append(_edge_row(record, score, "unlabeled_signature"))
            continue
        role = "disease_evidence" if score["evidence_score"] >= score["counterevidence_score"] else "disease_counterevidence"
        for term in terms:
            edges.append(_edge_row(record, score, role, term))

    edge_table = pd.DataFrame(edges, columns=EDGE_COLUMNS)
    if edge_table.empty:
        edge_table = empty_table(EDGE_COLUMNS)
    human_diseases = aggregate_diseases(edge_table, "human")
    mouse_diseases = aggregate_diseases(edge_table, "mouse")
    human_mechanisms = aggregate_mechanisms(edge_table, "human")
    mouse_mechanisms = aggregate_mechanisms(edge_table, "mouse")
    if rejected:
        warnings.append(f"Excluded {len(rejected)} public evidence record(s): {sorted(set(rejected))[0]}.")
    if human_diseases.empty:
        warnings.append("No ontology-normalized human disease evidence matched the current cache.")
    if not orthologs:
        warnings.append("Mouse model hypotheses were not scored because no one-to-one ortholog map was available.")

    manifest = {
        "cache_root": str(cache_root),
        "cache_manifest": cache_manifest,
        "accepted_record_count": accepted,
        "rejected_record_count": len(rejected),
        "human_query_gene_count": int((human_query != 0).sum()),
        "mouse_query_gene_count": int((mouse_query != 0).sum()),
        "species_are_separate": True,
    }
    return {
        "human_diseases": human_diseases,
        "mouse_diseases": mouse_diseases,
        "human_mechanisms": human_mechanisms,
        "mouse_mechanisms": mouse_mechanisms,
        "edges": edge_table,
        "manifest": manifest,
        "warnings": warnings,
    }


def _graphql_request(endpoint: str, query: str, variables: dict[str, Any], timeout: int = 45) -> dict[str, Any]:
    payload = json.dumps({"query": query, "variables": variables}).encode()
    request = urllib.request.Request(
        endpoint,
        data=payload,
        headers={"Accept": "application/json", "Content-Type": "application/json"},
        method="POST",
    )
    with urllib.request.urlopen(request, timeout=timeout) as response:
        result = json.loads(response.read().decode())
    if result.get("errors"):
        raise RuntimeError(result["errors"][0].get("message", "GraphQL query failed."))
    return result.get("data", {})


def _json_request(url: str, payload: dict[str, Any] | None = None, timeout: int = 45) -> dict[str, Any]:
    data = json.dumps(payload).encode() if payload is not None else None
    request = urllib.request.Request(
        url,
        data=data,
        headers={"Accept": "application/json", "Content-Type": "application/json"},
        method="POST" if data is not None else "GET",
    )
    with urllib.request.urlopen(request, timeout=timeout) as response:
        result = json.loads(response.read().decode())
    if not isinstance(result, dict):
        raise RuntimeError("Expected a JSON object from the orthology service.")
    return result


def _top_signed_genes(query: pd.Series, direction: int, limit: int) -> list[str]:
    values = query.loc[query * direction > 0].sort_values(ascending=direction < 0)
    return [display_gene(gene) for gene in values.head(limit).index if display_gene(gene)]


def _direction_from_label(label: str) -> str:
    match = re.search(r"(?:^|[\\s:_-])(up|down|dn)(?:$|[\\s:_-])", label.casefold())
    if not match:
        return "up"
    return "down" if match.group(1) in {"down", "dn"} else "up"


def _genes_from_nodes(value: Any) -> list[str]:
    if not isinstance(value, dict):
        return []
    nodes = value.get("nodes", [])
    if not isinstance(nodes, list):
        return []
    return _as_gene_list([node.get("symbol") for node in nodes if isinstance(node, dict)])


def _gse_ids(value: Any) -> list[str]:
    if not isinstance(value, dict):
        return []
    ids: list[str] = []
    for node in value.get("nodes", []):
        if not isinstance(node, dict):
            continue
        ids.extend(part.strip() for part in str(node.get("gse", "")).split(",") if part.strip())
    return sorted(set(ids))


def _fetch_text(url: str, timeout: int = 45) -> str:
    request = urllib.request.Request(url, headers={"User-Agent": "DSPIN-public-evidence/1.0"})
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return response.read().decode(errors="replace")


def _parse_mondo_obo(text: str) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    current: dict[str, Any] = {}
    for line in text.splitlines() + ["[Term]"]:
        if line == "[Term]":
            if current.get("id") and current.get("name") and not current.get("obsolete"):
                records.append(current)
            current = {"synonyms": [], "parents": [], "obsolete": False}
            continue
        if line.startswith("id: MONDO:"):
            current["id"] = line.split("id:", 1)[1].strip()
        elif line.startswith("name: "):
            current["name"] = line.split("name:", 1)[1].strip()
        elif line.startswith("synonym: "):
            match = re.match(r'synonym: "([^"]+)"', line)
            if match:
                current["synonyms"].append(match.group(1))
        elif line.startswith("is_a: MONDO:"):
            current["parents"].append(line.split("!", 1)[0].split("is_a:", 1)[1].strip())
        elif line == "is_obsolete: true":
            current["obsolete"] = True
    by_id = {record["id"]: record for record in records}
    for record in records:
        parents = [by_id[parent]["name"] for parent in record["parents"] if parent in by_id]
        record["family"] = parents[0] if parents else record["name"]
        record.pop("parents", None)
        record.pop("obsolete", None)
    return records


def load_disease_ontology(cache_root: Path) -> list[dict[str, Any]]:
    path = cache_path(cache_root, ONTOLOGY_FILE)
    if not path.exists():
        return []
    try:
        data = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError):
        return []
    terms = data.get("terms", []) if isinstance(data, dict) else []
    return [term for term in terms if isinstance(term, dict) and term.get("id") and term.get("name")]


def refresh_mondo_ontology(cache_root: Path, timeout: int = 45) -> tuple[list[dict[str, Any]], str | None]:
    try:
        terms = _parse_mondo_obo(_fetch_text(MONDO_OBO_URL, timeout=timeout))
    except (urllib.error.URLError, TimeoutError, OSError, RuntimeError) as exc:
        return load_disease_ontology(cache_root), f"MONDO refresh failed: {exc}"
    cache_root.mkdir(parents=True, exist_ok=True)
    cache_path(cache_root, ONTOLOGY_FILE).write_text(
        json.dumps({"source_url": MONDO_OBO_URL, "retrieved_at": datetime.now().isoformat(timespec="seconds"), "terms": terms})
    )
    return terms, None


def build_disease_term_index(ontology: list[dict[str, Any]]) -> dict[str, list[tuple[dict[str, Any], str]]]:
    labels: list[tuple[dict[str, Any], str, list[str]]] = []
    token_counts: Counter[str] = Counter()
    for term in ontology:
        term_id = str(term.get("id", ""))
        term_name = str(term.get("name", ""))
        if not is_specific_disease_term(term_id, term_name):
            continue
        for label in [term_name, *[str(value) for value in term.get("synonyms", [])]]:
            label = label.strip()
            tokens = re.findall(r"[a-z0-9]+", label.casefold())
            if len(label) < 6 or not tokens:
                continue
            labels.append((term, label, tokens))
            token_counts.update(set(tokens))

    index: dict[str, list[tuple[dict[str, Any], str]]] = {}
    for term, label, tokens in labels:
        anchor = min(tokens, key=lambda token: (token_counts[token], token))
        index.setdefault(anchor, []).append((term, label))
    return index


def match_disease_terms(
    text: str,
    ontology: list[dict[str, Any]],
    limit: int = 3,
    term_index: dict[str, list[tuple[dict[str, Any], str]]] | None = None,
) -> list[dict[str, str]]:
    if not text or not ontology:
        return []
    haystack = text.casefold()
    index = term_index if term_index is not None else build_disease_term_index(ontology)
    candidates = [
        candidate
        for token in set(re.findall(r"[a-z0-9]+", haystack))
        for candidate in index.get(token, [])
    ]
    matched: dict[str, dict[str, str]] = {}
    for term, label in candidates:
        term_id = str(term.get("id", ""))
        term_name = str(term.get("name", ""))
        if term_id in matched or not is_specific_disease_term(term_id, term_name):
            continue
        normalized_label = label.casefold()
        if normalized_label in haystack and re.search(
            rf"(?<![a-z0-9]){re.escape(normalized_label)}(?![a-z0-9])", haystack
        ):
            matched[term_id] = {
                "id": term_id,
                "name": term_name,
                "family": str(term.get("family", "")) or term_name,
            }
    return sorted(matched.values(), key=lambda term: (-len(term["name"]), term["name"]))[:limit]


def _geo_metadata(gse: str, timeout: int) -> dict[str, str]:
    try:
        text = _fetch_text(GEO_SOFT_URL.format(gse=urllib.parse.quote(gse)), timeout=timeout)
    except (urllib.error.URLError, TimeoutError, OSError):
        return {"title": "", "summary": ""}
    values: dict[str, list[str]] = {"title": [], "summary": []}
    for line in text.splitlines():
        if line.startswith("!Series_title ="):
            values["title"].append(line.split("=", 1)[1].strip())
        elif line.startswith("!Series_summary =") or line.startswith("!Series_overall_design ="):
            values["summary"].append(line.split("=", 1)[1].strip())
    return {key: " ".join(value) for key, value in values.items()}


def _rummageo_records(
    queries: dict[str, pd.Series],
    ontology: list[dict[str, Any]],
    limit: int,
    timeout: int,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    backgrounds = _graphql_request(
        RUMMAGEO_GRAPHQL,
        "{ backgrounds(first: 10) { nodes { id species } } }",
        {},
        timeout,
    ).get("backgrounds", {}).get("nodes", [])
    background_ids = {canonical_species(node.get("species")): node.get("id") for node in backgrounds if isinstance(node, dict)}
    query = """
    query($id: UUID!, $genes: [String]!, $first: Int!) {
      background(id: $id) {
        enrich(genes: $genes, first: $first) {
          nodes {
            pvalue adjPvalue oddsRatio nOverlap
            geneSet {
              id term species
              genes(first: 1000) { nodes { symbol } }
              geneSetGsesById { nodes { gse species } }
            }
          }
        }
      }
    }
    """
    records: list[dict[str, Any]] = []
    metadata_cache: dict[str, dict[str, str]] = {}
    term_index = build_disease_term_index(ontology)
    for species, patient_query in queries.items():
        background_id = background_ids.get(species)
        if not background_id or patient_query.empty:
            continue
        for direction in (1, -1):
            genes = _top_signed_genes(patient_query, direction, 75)
            if not genes:
                continue
            data = _graphql_request(
                RUMMAGEO_GRAPHQL,
                query,
                {"id": background_id, "genes": genes, "first": limit},
                timeout,
            )
            nodes = data.get("background", {}).get("enrich", {}).get("nodes", [])
            for result in nodes:
                gene_set = result.get("geneSet", {}) if isinstance(result, dict) else {}
                if canonical_species(gene_set.get("species")) != species:
                    continue
                symbols = _genes_from_nodes(gene_set.get("genes"))
                source_direction = _direction_from_label(str(gene_set.get("term", "")))
                gse_ids = _gse_ids(gene_set.get("geneSetGsesById"))
                metadata = []
                for gse in gse_ids[:3]:
                    if gse not in metadata_cache:
                        metadata_cache[gse] = _geo_metadata(gse, timeout)
                    metadata.append(metadata_cache[gse])
                title = " ".join(item.get("title", "") for item in metadata).strip()
                summary = " ".join(item.get("summary", "") for item in metadata).strip()
                source_text = " ".join([str(gene_set.get("term", "")), title, summary])
                record = {
                    "record_id": f"rummageo:{gene_set.get('id')}",
                    "resource": "rummageo",
                    "species": species,
                    "species_provenance": "source_background",
                    "kind": "disease_signature",
                    "up_genes": symbols if source_direction == "up" else [],
                    "down_genes": symbols if source_direction == "down" else [],
                    "study_id": ",".join(gse_ids),
                    "source_label": str(gene_set.get("term", "")),
                    "source_url": f"https://rummageo.com/gene-set/{gene_set.get('id')}",
                    "source_confidence": 1.0,
                    "disease_terms": match_disease_terms(source_text, ontology, term_index=term_index),
                    "metadata_text": source_text,
                    "source_stats": {
                        key: result.get(key) for key in ("pvalue", "adjPvalue", "oddsRatio", "nOverlap")
                    },
                }
                records.append(record)
    return records, {"source_url": RUMMAGEO_GRAPHQL, "species": ["human", "mouse"]}


def _perturbseqr_records(
    queries: dict[str, pd.Series],
    limit: int,
    timeout: int,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    query = """
    query($genes: [String], $libraries: [String], $first: Int) {
      geneSetGeneSearch(genes: $genes, libraryNames: $libraries, first: $first) {
        nodes {
          id term perturbation
          library { name }
          genes(first: 1000) { nodes { symbol } }
        }
      }
    }
    """
    allowed_libraries = {"Perturb Atlas Human": "human", "Perturb Atlas Mouse": "mouse"}
    records: list[dict[str, Any]] = []
    for query_species, patient_query in queries.items():
        if patient_query.empty:
            continue
        library_name = f"Perturb Atlas {query_species.title()}"
        if allowed_libraries.get(library_name) != query_species:
            continue
        for direction in (1, -1):
            genes = _top_signed_genes(patient_query, direction, 75)
            if not genes:
                continue
            data = _graphql_request(
                PERTURBSEQR_GRAPHQL,
                query,
                {"genes": genes, "libraries": [library_name], "first": limit},
                timeout,
            )
            search = data.get("geneSetGeneSearch") if isinstance(data, dict) else None
            nodes = search.get("nodes") if isinstance(search, dict) else []
            for gene_set in nodes if isinstance(nodes, list) else []:
                if not isinstance(gene_set, dict):
                    continue
                library = gene_set.get("library")
                result_library = str(library.get("name", "")) if isinstance(library, dict) else ""
                species = allowed_libraries.get(result_library)
                if species != query_species:
                    continue
                symbols = _genes_from_nodes(gene_set.get("genes"))
                source_direction = _direction_from_label(str(gene_set.get("term", "")))
                records.append(
                    {
                        "record_id": f"perturbseqr:{gene_set.get('id')}",
                        "resource": "perturbseqr",
                        "species": species,
                        "species_provenance": "library_name",
                        "kind": "gene_perturbation",
                        "up_genes": symbols if source_direction == "up" else [],
                        "down_genes": symbols if source_direction == "down" else [],
                        "source_label": str(gene_set.get("term", "")),
                        "source_url": f"https://perturbseqr.maayanlab.cloud/gene-set/{gene_set.get('id')}",
                        "perturbation_gene": str(gene_set.get("perturbation", "")),
                        "perturbation_mode": "",
                        "library_name": result_library,
                        "source_confidence": 1.0,
                    }
                )
    return records, {
        "source_url": PERTURBSEQR_GRAPHQL,
        "species": ["human", "mouse"],
        "chemical_data": "excluded",
        "query_mode": "library_restricted_gene_set_search",
    }


def refresh_public_evidence_cache(
    cache_root: Path,
    signed_gene_scores: pd.Series,
    orthologs: dict[str, str],
    limit: int = 25,
    timeout: int = 45,
) -> tuple[dict[str, Any], list[str]]:
    existing, manifest, warnings = read_evidence_cache(cache_root)
    warnings = [warning for warning in warnings if "cache is missing" not in warning]
    ontology, ontology_warning = refresh_mondo_ontology(cache_root, timeout=timeout)
    if ontology_warning:
        warnings.append(ontology_warning)
    human_query = normalize_signed_series(signed_gene_scores)
    queries = {
        "human": human_query,
        "mouse": map_human_query_to_mouse(human_query, orthologs) if orthologs else pd.Series(dtype=float),
    }
    refresh_records: list[dict[str, Any]] = []
    sources = dict(manifest.get("sources", {}))
    refreshed_at = datetime.now().isoformat(timespec="seconds")
    try:
        records, source_info = _rummageo_records(queries, ontology, limit, timeout)
        refresh_records.extend(records)
        sources["rummageo"] = {
            **source_info,
            "source_version": "live_graphql",
            "retrieved_at": refreshed_at,
            "status": "refreshed",
            "record_count": len(records),
            "records_sha256": records_checksum(records),
        }
    except (urllib.error.URLError, TimeoutError, OSError, RuntimeError) as exc:
        warnings.append(f"RummaGEO refresh failed; retained cached evidence: {exc}")
        sources["rummageo"] = {
            "source_url": RUMMAGEO_GRAPHQL,
            "species": ["human", "mouse"],
            "source_version": "live_graphql",
            "retrieved_at": refreshed_at,
            "status": "refresh_failed",
            "record_count": 0,
        }

    sources["rummagene"] = {
        "source_url": RUMMAGENE_GRAPHQL,
        "source_version": "live_graphql",
        "retrieved_at": refreshed_at,
        "status": "skipped_unverified_species_provenance",
        "record_count": 0,
        "species": [],
        "records_sha256": records_checksum([]),
    }
    warnings.append(
        "Live Rummagene refresh skipped because its GraphQL response does not expose explicit source species provenance; "
        "explicitly tagged cached Rummagene records remain eligible."
    )

    try:
        records, source_info = _perturbseqr_records(queries, limit, timeout)
        refresh_records.extend(records)
        sources["perturbseqr"] = {
            **source_info,
            "source_version": "live_graphql",
            "retrieved_at": refreshed_at,
            "status": "refreshed",
            "record_count": len(records),
            "records_sha256": records_checksum(records),
        }
    except (urllib.error.URLError, TimeoutError, OSError, RuntimeError) as exc:
        warnings.append(f"Perturb-Seqr refresh failed; retained cached evidence: {exc}")
        sources["perturbseqr"] = {
            "source_url": PERTURBSEQR_GRAPHQL,
            "species": ["human", "mouse"],
            "source_version": "live_graphql",
            "retrieved_at": refreshed_at,
            "status": "refresh_failed",
            "record_count": 0,
        }

    merged = {str(record.get("record_id") or record.get("id")): record for record in existing if isinstance(record, dict)}
    for record in refresh_records:
        merged[str(record["record_id"])] = record
    updated_manifest = write_evidence_cache(
        cache_root,
        list(merged.values()),
        {**manifest, "sources": sources},
    )
    return updated_manifest, warnings


def refresh_human_perturbseqr_cache(
    cache_root: Path,
    signed_gene_scores: pd.Series,
    limit: int = 25,
    timeout: int = 45,
) -> tuple[dict[str, Any], list[str]]:
    """Refresh only explicitly human Perturb-Seqr gene perturbations.

    The retired disease scorer remains available above for historical reference,
    but the active interpreter no longer needs mouse orthologs or disease
    signatures when refreshing mechanistic Perturb-Seqr evidence.
    """
    existing, manifest, warnings = read_evidence_cache(cache_root)
    warnings = [warning for warning in warnings if "cache is missing" not in warning]
    refreshed_at = datetime.now().isoformat(timespec="seconds")
    try:
        records, source_info = _perturbseqr_records(
            {"human": normalize_signed_series(signed_gene_scores)},
            limit,
            timeout,
        )
    except (urllib.error.URLError, TimeoutError, OSError, RuntimeError) as exc:
        warnings.append(f"Human Perturb-Seqr refresh failed; retained cached evidence: {exc}")
        records = []
        source_info = {
            "source_url": PERTURBSEQR_GRAPHQL,
            "species": ["human"],
            "chemical_data": "excluded",
            "query_mode": "library_restricted_gene_set_search",
        }
        status = "refresh_failed"
    else:
        status = "refreshed"

    merged = {str(record.get("record_id") or record.get("id")): record for record in existing if isinstance(record, dict)}
    for record in records:
        if canonical_species(record.get("species")) == "human":
            merged[str(record["record_id"])] = record
    sources = dict(manifest.get("sources", {}))
    sources["perturbseqr"] = {
        **source_info,
        "species": ["human"],
        "source_version": "live_graphql",
        "retrieved_at": refreshed_at,
        "status": status,
        "record_count": len(records),
        "records_sha256": records_checksum(records),
    }
    updated = write_evidence_cache(cache_root, list(merged.values()), {**manifest, "sources": sources})
    return updated, warnings
