"""Human-only KG disease and mimic evidence retrieval for the DSPIN interpreter.

This module intentionally does not score cytokines, DSPIN programs, DSPIN genes, or
cell types. It consumes those completed outputs and replaces only disease/mimic
hypothesis generation.
"""

from __future__ import annotations

import hashlib
import json
import math
import os
import re
import urllib.error
import urllib.request
import zipfile
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

try:
    from scripts import public_evidence
except ImportError:  # Supports direct execution via `python scripts/...`.
    import public_evidence


CACHE_SCHEMA_VERSION = 1
MANIFEST_FILE = "kg_manifest.json"
MONDO_TERMS_FILE = "mondo_terms.json"
ENRICHR_DOWNLOADS_URL = "https://s3.amazonaws.com/maayan-kg/enrichr-kg/assets/downloads.json"
MONDO_OBO_URL = "https://purl.obolibrary.org/obo/mondo.obo"
ADAPTIVE_PVALUE_MAX_PERMUTATIONS = 100_000
ADAPTIVE_PVALUE_BATCH_PERMUTATIONS = 2_000
ADAPTIVE_PVALUE_ROLES = {"direct_disease", "cancer_context"}

TERM_COLUMNS = [
    "species",
    "role",
    "library",
    "term_id",
    "term_label",
    "term_size",
    "available_gene_count",
    "direct_score",
    "direct_effect_z",
    "direct_empirical_p",
    "direct_bh_fdr",
    "direct_permutations",
    "direct_exceedance_count",
    "direct_pvalue_is_upper_bound",
    "dspin_depth_0_score",
    "dspin_depth_1_score",
    "dspin_depth_2_score",
    "dspin_depth_3_score",
    "dspin_score",
    "chea_depth_0_score",
    "chea_depth_1_score",
    "chea_depth_2_score",
    "chea_depth_3_score",
    "chea_score",
    "manifold_score",
    "signed_sensitivity_score",
    "ppi_score",
    "coexpression_score",
    "leading_genes",
    "mondo_id",
    "mondo_name",
    "disease_family",
    "hypothesis_class",
]

DISEASE_COLUMNS = [
    "rank",
    "disease_family",
    "disease_id",
    "disease_name",
    "hypothesis_class",
    "combined_p_value",
    "combined_fdr",
    "evidence_score",
    "evidence_tier",
    "direct_effect_z",
    "direct_permutations",
    "combined_p_value_is_upper_bound",
    "direct_resource_count",
    "supporting_sources",
    "matched_genes",
    "dspin_score",
    "chea_score",
    "manifold_score",
    "ppi_score",
    "coexpression_score",
    "network_primary_promoted",
]

MIMIC_COLUMNS = [
    "rank",
    "disease_family",
    "disease_id",
    "disease_name",
    "hypothesis_class",
    "combined_p_value",
    "combined_fdr",
    "evidence_score",
    "evidence_tier",
    "direct_effect_z",
    "direct_permutations",
    "combined_p_value_is_upper_bound",
    "direct_resource_count",
    "supporting_sources",
    "matched_genes",
    "dspin_score",
    "chea_score",
    "manifold_score",
    "ppi_score",
    "coexpression_score",
    "network_primary_promoted",
]

EDGE_COLUMNS = [
    "species",
    "evidence_kind",
    "role",
    "resource",
    "term_id",
    "term_label",
    "disease_id",
    "disease_name",
    "patient_gene",
    "patient_signed_score",
    "direct_score",
    "dspin_score",
    "chea_score",
    "manifold_score",
    "ppi_score",
    "coexpression_score",
]

CONCORDANCE_COLUMNS = [
    "row_type",
    "network_variant",
    "source_gene",
    "target_gene",
    "dspin_weight",
    "chea_weight",
    "source_signal",
    "target_signal",
    "sign_concordant",
    "overlapping_edges",
    "topology_correlation",
]

MODULE_COLUMNS = [
    "module_id",
    "gene_count",
    "genes",
    "mean_manifold_similarity",
    "mean_signed_sensitivity_similarity",
    "supporting_views",
]


# Enrichr-KG is a heterogeneous term-gene knowledge graph, not one GRN. The
# `label_human` resources contain mixed-species terms, so they are filtered per
# term. `source_human` resources have documented human source-library provenance.
ASSET_SPECS: dict[str, dict[str, Any]] = {
    "chea_node_draw": {
        "kind": "chea",
        "url": "https://s3.amazonaws.com/maayan-kg/chea-kg/node_draw_GRN.zip",
        "filename": "chea_node_draw.zip",
        "species_policy": "human_tf_grn",
    },
    "chea_target_set_swap": {
        "kind": "chea_qc",
        "url": "https://s3.amazonaws.com/maayan-kg/chea-kg/target_set_swap_GRN.zip",
        "filename": "chea_target_set_swap.zip",
        "species_policy": "human_tf_grn",
    },
    "chea_unfiltered": {
        "kind": "chea_qc",
        "url": "https://s3.amazonaws.com/maayan-kg/chea-kg/unfiltered_GRN.zip",
        "filename": "chea_unfiltered.zip",
        "species_policy": "human_tf_grn",
    },
    "disgenet": {
        "source": "DisGeNET",
        "relation": "DisGeNET_Association",
        "filename": "disgenet.csv",
        "role": "direct_disease",
        "species_policy": "source_human",
    },
    "jensen_diseases": {
        "source": "Jensen_DISEASES",
        "relation": "Jensen_Disease",
        "filename": "jensen_diseases.csv",
        "role": "direct_disease",
        "species_policy": "source_human",
    },
    "gwas_catalog": {
        "source": "GWAS_Catalog_2019",
        "relation": "GWAS_Catalog",
        "filename": "gwas_catalog.csv",
        "role": "direct_disease",
        "species_policy": "source_human",
    },
    "hpo": {
        "source": "Human_Phenotype_Ontology",
        "relation": "HPO",
        "filename": "hpo.csv",
        "role": "phenotype",
        "species_policy": "source_human",
    },
    "go_bp": {
        "source": "GO_Biological_Process_2021",
        "relation": "GO_BP",
        "filename": "go_bp.csv",
        "role": "pathway",
        "species_policy": "source_human",
    },
    "reactome": {
        "source": "Reactome_2022",
        "relation": "Reactome",
        "filename": "reactome.csv",
        "role": "pathway",
        "species_policy": "source_human",
    },
    "kegg_human": {
        "source": "KEGG_2021_Human",
        "relation": "KEGG_Pathway",
        "filename": "kegg_human.csv",
        "role": "pathway",
        "species_policy": "source_human",
    },
    "wikipathways_human": {
        "source": "WikiPathway_2021_Human",
        "relation": "WikiPathways",
        "filename": "wikipathways_human.csv",
        "role": "pathway",
        "species_policy": "source_human",
    },
    "chea_2022": {
        "source": "ChEA_2022",
        "relation": "ChEA_2022_TF",
        "filename": "chea_2022.csv",
        "role": "regulatory",
        "species_policy": "label_human",
    },
    "trrust": {
        "source": "TRRUST_Transcription_Factors_2019",
        "relation": "TRRUST_TF",
        "filename": "trrust.csv",
        "role": "regulatory",
        "species_policy": "label_human",
    },
    "archs4_tf": {
        "source": "ARCHS4_TFs_Coexp",
        "relation": "ARCHS4_TFs_Coexpression",
        "filename": "archs4_tf.csv",
        "role": "regulatory",
        "species_policy": "label_human",
    },
    "tabula_sapiens": {
        "source": "Tabula_Sapiens",
        "relation": "Tabula_Sapiens_Association",
        "filename": "tabula_sapiens.csv",
        "role": "cell_context",
        "species_policy": "source_human",
    },
    "hubmap": {
        "source": "HuBMAP_ASCTplusB_augmented_2022",
        "relation": "HuBMAP_ASCTplusB_augmented",
        "filename": "hubmap.csv",
        "role": "cell_context",
        "species_policy": "source_human",
    },
    "descartes": {
        "source": "Descartes_Cell_Types_and_Tissue_2021",
        "relation": "Descartes_Cell_Types_and_Tissue",
        "filename": "descartes.csv",
        "role": "cell_context",
        "species_policy": "source_human",
    },
    "human_gene_atlas": {
        "source": "Human_Gene_Atlas",
        "relation": "Human_Gene_Atlas",
        "filename": "human_gene_atlas.csv",
        "role": "cell_context",
        "species_policy": "source_human",
    },
    "lincs_crispr_up": {
        "source": "LINCS_L1000_CRISPR_KO_Consensus_Sigs",
        "relation": "LINCS_L1000_CRISPR_KO_Up-regulated",
        "filename": "lincs_crispr_up.csv",
        "role": "perturbation",
        "direction": 1.0,
        "species_policy": "source_human",
    },
    "lincs_crispr_down": {
        "source": "LINCS_L1000_CRISPR_KO_Consensus_Sigs",
        "relation": "LINCS_L1000_CRISPR_KO_Down-regulated",
        "filename": "lincs_crispr_down.csv",
        "role": "perturbation",
        "direction": -1.0,
        "species_policy": "source_human",
    },
    "string_ppi": {
        "source": "Gene",
        "relation": "STRING-db PPI",
        "filename": "string_ppi.csv",
        "kind": "ppi",
        "species_policy": "source_human",
    },
    "ccle": {
        "source": "CCLE_Proteomics_2020",
        "relation": "CCLE_Proteomics",
        "filename": "ccle.csv",
        "role": "cancer_context",
        "species_policy": "source_human",
    },
    "achilles_decrease": {
        "source": "Achilles_Cell_Line",
        "relation": "Achilles_(fitness_decrease)",
        "filename": "achilles_decrease.csv",
        "role": "cancer_context",
        "species_policy": "source_human",
    },
    "achilles_increase": {
        "source": "Achilles_Cell_Line",
        "relation": "Achilles_(fitness_increase)",
        "filename": "achilles_increase.csv",
        "role": "cancer_context",
        "species_policy": "source_human",
    },
}


def _now() -> str:
    return datetime.now().isoformat(timespec="seconds")


def _canonical_gene(value: Any) -> str:
    return str(value or "").strip().upper()


def _canonical_label(value: Any) -> str:
    text = str(value or "").casefold().strip()
    text = re.sub(r"\s*\([^)]*(?:mondo|mesh|omim|orpha|hp:)[^)]*\)", "", text)
    text = re.sub(r"[^a-z0-9]+", " ", text)
    return re.sub(r"\s+", " ", text).strip()


def _empty(columns: list[str]) -> pd.DataFrame:
    return pd.DataFrame(columns=columns)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _read_json(path: Path, default: dict[str, Any]) -> dict[str, Any]:
    if not path.exists():
        return default
    try:
        payload = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError):
        return default
    return payload if isinstance(payload, dict) else default


def _write_kg_manifest(cache_root: Path, manifest: dict[str, Any]) -> None:
    path = cache_root / MANIFEST_FILE
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(manifest, indent=2) + "\n")
    os.replace(temporary, path)


def _asset_path(cache_root: Path, manifest: dict[str, Any], asset_id: str) -> Path | None:
    entry = manifest.get("assets", {}).get(asset_id, {})
    raw_path = entry.get("path") if isinstance(entry, dict) else None
    if not raw_path:
        return None
    path = Path(raw_path)
    if not path.is_absolute():
        path = cache_root / path
    return path if path.exists() else None


def _fetch_json(url: str, timeout: int = 60) -> dict[str, Any]:
    request = urllib.request.Request(url, headers={"User-Agent": "DSPIN-human-kg/1.0"})
    with urllib.request.urlopen(request, timeout=timeout) as response:
        payload = json.loads(response.read().decode(errors="replace"))
    if not isinstance(payload, dict):
        raise ValueError(f"Expected an object from {url}")
    return payload


def _download(url: str, path: Path, timeout: int = 120) -> dict[str, Any]:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".part")
    request = urllib.request.Request(url, headers={"User-Agent": "DSPIN-human-kg/1.0"})
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response, temporary.open("wb") as handle:
            while True:
                chunk = response.read(1024 * 1024)
                if not chunk:
                    break
                handle.write(chunk)
            etag = response.headers.get("ETag", "").strip('"')
    except Exception:
        temporary.unlink(missing_ok=True)
        raise
    os.replace(temporary, path)
    return {"etag": etag, "bytes": path.stat().st_size, "sha256": _sha256(path)}


def _enrichr_asset_urls(downloads: dict[str, Any]) -> dict[tuple[str, str], dict[str, Any]]:
    lookup: dict[tuple[str, str], dict[str, Any]] = {}
    for entry in downloads.get("edges", []):
        if not isinstance(entry, dict):
            continue
        source = str(entry.get("source", ""))
        relation = str(entry.get("relation", ""))
        url = str(entry.get("url", ""))
        if source and relation and url:
            lookup[(source, relation)] = entry
    return lookup


def refresh_kg_cache(cache_root: Path, timeout: int = 120) -> tuple[dict[str, Any], list[str]]:
    """Refresh static public KG assets. No patient-derived value is transmitted."""
    cache_root.mkdir(parents=True, exist_ok=True)
    assets_dir = cache_root / "assets"
    existing = _read_json(cache_root / MANIFEST_FILE, {})
    manifest: dict[str, Any] = {
        "schema_version": CACHE_SCHEMA_VERSION,
        "engine": "human_kg_multiview",
        "species": "human",
        "patient_signature_transmitted": False,
        "refreshed_at": _now(),
        "assets": dict(existing.get("assets", {})),
        "network_primary_validation": existing.get(
            "network_primary_validation",
            {"status": "not_run", "approved": False},
        ),
    }
    warnings: list[str] = []
    try:
        enrichr_index = _enrichr_asset_urls(_fetch_json(ENRICHR_DOWNLOADS_URL, timeout=timeout))
    except (urllib.error.URLError, TimeoutError, OSError, ValueError) as exc:
        enrichr_index = {}
        warnings.append(f"Enrichr-KG download index refresh failed: {exc}")

    for asset_id, spec in ASSET_SPECS.items():
        url = str(spec.get("url", ""))
        if not url:
            entry = enrichr_index.get((str(spec["source"]), str(spec["relation"])))
            url = str(entry.get("url", "")) if entry else ""
        if not url:
            warnings.append(f"No Enrichr-KG download URL was found for {asset_id}; retained any cached copy.")
            continue
        destination = assets_dir / str(spec["filename"])
        if destination.exists():
            previous = manifest["assets"].get(asset_id, {})
            metadata = {
                "etag": previous.get("etag", ""),
                "bytes": destination.stat().st_size,
                "sha256": _sha256(destination),
            }
        else:
            try:
                metadata = _download(url, destination, timeout=timeout)
            except (urllib.error.URLError, TimeoutError, OSError) as exc:
                warnings.append(f"Static KG download failed for {asset_id}; retained any cached copy: {exc}")
                continue
        manifest["assets"][asset_id] = {
            "path": str(destination.relative_to(cache_root)),
            "source_url": url,
            "retrieved_at": _now(),
            "source_version": "current",
            "species": "human",
            "species_policy": spec["species_policy"],
            **metadata,
        }

    mondo_path = assets_dir / "mondo.obo"
    try:
        mondo_metadata = (
            {
                "etag": manifest.get("mondo", {}).get("etag", ""),
                "bytes": mondo_path.stat().st_size,
                "sha256": _sha256(mondo_path),
            }
            if mondo_path.exists()
            else _download(MONDO_OBO_URL, mondo_path, timeout=timeout)
        )
        terms = _parse_mondo_obo(mondo_path.read_text(errors="replace"))
        term_path = cache_root / MONDO_TERMS_FILE
        term_path.write_text(json.dumps({"terms": terms}, indent=2) + "\n")
        manifest["mondo"] = {
            "path": str(mondo_path.relative_to(cache_root)),
            "terms_path": MONDO_TERMS_FILE,
            "source_url": MONDO_OBO_URL,
            "retrieved_at": _now(),
            "source_version": "current",
            **mondo_metadata,
        }
    except (urllib.error.URLError, TimeoutError, OSError) as exc:
        warnings.append(f"MONDO refresh failed; retained any cached ontology: {exc}")

    _write_kg_manifest(cache_root, manifest)
    return manifest, warnings


def read_kg_cache(cache_root: Path) -> tuple[dict[str, Any], list[str]]:
    manifest_path = cache_root / MANIFEST_FILE
    if not manifest_path.exists():
        return {}, [f"Human KG cache is missing at {cache_root}; run once with --refresh-kg."]
    manifest = _read_json(manifest_path, {})
    if not manifest:
        return {}, [f"Human KG manifest is unreadable at {manifest_path}."]
    if manifest.get("species") != "human":
        return {}, [f"Human KG manifest at {manifest_path} is not human-only and was rejected."]
    assets = manifest.get("assets")
    if isinstance(assets, dict) and assets:
        return manifest, []

    recovered: dict[str, dict[str, Any]] = {}
    for asset_id, spec in ASSET_SPECS.items():
        path = cache_root / "assets" / str(spec["filename"])
        if not path.exists():
            continue
        recovered[asset_id] = {
            "path": str(path.relative_to(cache_root)),
            "source_url": str(spec.get("url", "cached Enrichr-KG asset")),
            "source_version": "recovered_local_file",
            "species": "human",
            "species_policy": spec["species_policy"],
            "bytes": path.stat().st_size,
            "sha256": _sha256(path),
        }
    if recovered:
        recovered_manifest = dict(manifest)
        recovered_manifest["assets"] = recovered
        return recovered_manifest, [
            "KG manifest had no asset entries; recovered verified local human assets from the cache directory."
        ]
    return manifest, []


def _human_term_allowed(label: str, policy: str) -> bool:
    lowered = str(label or "").casefold()
    if re.search(r"\b(mouse|murine|mus musculus|mixed|unknown)\b", lowered):
        return False
    if policy == "label_human":
        return bool(re.search(r"\b(human|homo sapiens)\b", lowered))
    return True


def _explicitly_nonhuman(value: Any) -> bool:
    return bool(re.search(r"\b(mouse|murine|mus musculus|mixed|unknown)\b", str(value or "").casefold()))


def _csv_columns(path: Path) -> set[str]:
    return set(pd.read_csv(path, nrows=0).columns.astype(str))


def load_term_sets(
    path: Path,
    library: str,
    role: str,
    gene_index: pd.Index,
    species_policy: str,
    direction: float = 0.0,
    chunk_size: int = 100_000,
) -> tuple[list[dict[str, Any]], dict[str, int]]:
    """Read a term-gene CSV without loading genes outside the DSPIN gene universe."""
    columns = _csv_columns(path)
    required = {"source", "source_label", "target_label"}
    missing = required - columns
    if missing:
        raise ValueError(f"{path} is missing required columns: {', '.join(sorted(missing))}")
    aliases = {_canonical_gene(gene): str(gene) for gene in gene_index.astype(str)}
    terms: dict[str, dict[str, Any]] = {}
    stats = {"accepted_edges": 0, "rejected_species_edges": 0, "universe_edges": 0}
    usecols = ["source", "source_label", "target_label"]
    for chunk in pd.read_csv(path, usecols=usecols, dtype=str, chunksize=chunk_size):
        for row in chunk.itertuples(index=False):
            term_id, label, target = (str(value) for value in row)
            if not _human_term_allowed(label, species_policy):
                stats["rejected_species_edges"] += 1
                continue
            stats["accepted_edges"] += 1
            term = terms.setdefault(
                term_id,
                {
                    "library": library,
                    "role": role,
                    "term_id": term_id,
                    "term_label": label,
                    "term_size": 0,
                    "genes": set(),
                    "direction": float(direction),
                },
            )
            term["term_size"] += 1
            gene = aliases.get(_canonical_gene(target))
            if gene:
                term["genes"].add(gene)
                stats["universe_edges"] += 1
    return list(terms.values()), stats


def _aligned_matrix(network: pd.DataFrame | None, genes: list[str]) -> np.ndarray:
    if network is None or network.empty:
        return np.zeros((len(genes), len(genes)), dtype=float)
    frame = network.copy()
    frame.index = frame.index.astype(str)
    frame.columns = frame.columns.astype(str)
    matrix = frame.reindex(index=genes, columns=genes).apply(pd.to_numeric, errors="coerce").fillna(0.0)
    values = matrix.to_numpy(dtype=float, copy=True)
    np.fill_diagonal(values, 0.0)
    return values


def _load_chea_matrix(path: Path, genes: list[str]) -> np.ndarray:
    aliases = {_canonical_gene(gene): index for index, gene in enumerate(genes)}
    matrix = np.zeros((len(genes), len(genes)), dtype=float)
    with zipfile.ZipFile(path) as archive:
        for member in archive.namelist():
            lowered = member.casefold()
            if not lowered.endswith(".csv") or "node" in lowered or "/._" in lowered:
                continue
            if "upregulat" in lowered:
                edge_sign = 1.0
            elif "downregulat" in lowered:
                edge_sign = -1.0
            else:
                continue
            with archive.open(member) as handle:
                header = pd.read_csv(handle, nrows=0)
            usecols = [name for name in ["source_label", "target_label", "z_score", "species"] if name in header.columns]
            if not {"source_label", "target_label"}.issubset(usecols):
                continue
            with archive.open(member) as handle:
                for chunk in pd.read_csv(handle, usecols=usecols, chunksize=100_000):
                    for row in chunk.itertuples(index=False):
                        if _explicitly_nonhuman(getattr(row, "species", "")):
                            continue
                        source = aliases.get(_canonical_gene(getattr(row, "source_label")))
                        target = aliases.get(_canonical_gene(getattr(row, "target_label")))
                        if source is None or target is None or source == target:
                            continue
                        z_score = getattr(row, "z_score", 1.0)
                        try:
                            weight = abs(float(z_score))
                        except (TypeError, ValueError):
                            weight = 1.0
                        # Rows are targets and columns are sources, matching DSPIN J @ state.
                        matrix[target, source] += edge_sign * weight
    return matrix


def _load_ppi_matrix(path: Path, genes: list[str]) -> np.ndarray:
    columns = _csv_columns(path)
    required = {"source_label", "target_label"}
    if not required.issubset(columns):
        raise ValueError(f"{path} is not a source/target gene edge table.")
    aliases = {_canonical_gene(gene): index for index, gene in enumerate(genes)}
    matrix = np.zeros((len(genes), len(genes)), dtype=float)
    usecols = [name for name in ["source_label", "target_label", "combined_score"] if name in columns]
    for chunk in pd.read_csv(path, usecols=usecols, chunksize=100_000):
        for row in chunk.itertuples(index=False):
            source = aliases.get(_canonical_gene(getattr(row, "source_label")))
            target = aliases.get(_canonical_gene(getattr(row, "target_label")))
            if source is None or target is None or source == target:
                continue
            raw_weight = getattr(row, "combined_score", 1.0)
            try:
                weight = float(raw_weight)
            except (TypeError, ValueError):
                weight = 1.0
            if weight > 1.0:
                weight /= 1000.0
            if not math.isfinite(weight) or weight <= 0:
                continue
            matrix[source, target] = max(matrix[source, target], weight)
            matrix[target, source] = max(matrix[target, source], weight)
    return matrix


def _archs4_tf_coexpression_matrix(terms: list[dict[str, Any]], genes: list[str]) -> np.ndarray:
    """Project ARCHS4 TF coexpression term-gene edges to a symmetric support graph.

    This is deliberately unsigned and non-causal: it is a TF-gene coexpression
    projection, not a regulatory edge assertion.
    """
    aliases = {_canonical_gene(gene): index for index, gene in enumerate(genes)}
    matrix = np.zeros((len(genes), len(genes)), dtype=float)
    for term in terms:
        label = str(term["term_label"])
        tf = aliases.get(_canonical_gene(label.split()[0] if label else ""))
        members = [aliases[_canonical_gene(gene)] for gene in term["genes"] if _canonical_gene(gene) in aliases]
        if tf is None or not members:
            continue
        weight = 1.0 / math.sqrt(len(members))
        for target in members:
            if tf == target:
                continue
            matrix[tf, target] = max(matrix[tf, target], weight)
            matrix[target, tf] = max(matrix[target, tf], weight)
    return matrix


def _row_normalize_absolute(matrix: np.ndarray, symmetric: bool = False) -> np.ndarray:
    values = np.asarray(matrix, dtype=float).copy()
    np.fill_diagonal(values, 0.0)
    values = np.abs(values)
    if symmetric:
        values = np.maximum(values, values.T)
    totals = values.sum(axis=1, keepdims=True)
    return np.divide(values, totals, out=np.zeros_like(values), where=totals > 0)


def propagate_signed(seed: np.ndarray, matrix: np.ndarray, max_depth: int = 3) -> list[np.ndarray]:
    values = np.asarray(matrix, dtype=float).copy()
    np.fill_diagonal(values, 0.0)
    totals = np.abs(values).sum(axis=1, keepdims=True)
    normalized = np.divide(values, totals, out=np.zeros_like(values), where=totals > 0)
    states = [np.asarray(seed, dtype=float).copy()]
    for _ in range(max_depth):
        states.append(normalized.dot(states[-1]))
    return states


def _unsigned_proximity(seed: np.ndarray, matrix: np.ndarray) -> np.ndarray:
    normalized = _row_normalize_absolute(matrix, symmetric=True)
    return 0.5 * np.abs(seed) + 0.5 * normalized.dot(np.abs(seed))


def _randomized_embedding(features: np.ndarray, dimensions: int = 16, seed: int = 17) -> np.ndarray:
    if features.size == 0 or not np.any(features):
        return np.zeros((features.shape[0], 0), dtype=float)
    rank = min(dimensions, features.shape[0] - 1, features.shape[1])
    if rank <= 0:
        return np.zeros((features.shape[0], 0), dtype=float)
    oversample = min(8, max(2, features.shape[1] - rank))
    random = np.random.default_rng(seed).standard_normal((features.shape[1], rank + oversample))
    q, _ = np.linalg.qr(features.dot(random), mode="reduced")
    compressed = q.T.dot(features)
    left, singular, _ = np.linalg.svd(compressed, full_matrices=False)
    embedding = q.dot(left[:, :rank] * singular[:rank])
    norms = np.linalg.norm(embedding, axis=1, keepdims=True)
    return np.divide(embedding, norms, out=np.zeros_like(embedding), where=norms > 0)


def _query_embedding_similarity(embedding: np.ndarray, seed: np.ndarray) -> np.ndarray:
    if embedding.size == 0 or not np.any(embedding) or not np.any(seed):
        return np.zeros(seed.size, dtype=float)
    weights = np.abs(seed)
    query = weights.dot(embedding) / weights.sum()
    norm = float(np.linalg.norm(query))
    if norm <= 0:
        return np.zeros(seed.size, dtype=float)
    return embedding.dot(query / norm)


def _multi_view_embeddings(
    dspin: np.ndarray,
    chea: np.ndarray,
    ppi: np.ndarray,
    coexpression: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    unsigned_views = [
        _row_normalize_absolute(dspin, symmetric=True),
        _row_normalize_absolute(chea, symmetric=True),
        _row_normalize_absolute(ppi, symmetric=True),
        _row_normalize_absolute(coexpression, symmetric=True),
    ]
    unsigned = _randomized_embedding(np.concatenate(unsigned_views, axis=1), seed=17)
    signed_features: list[np.ndarray] = []
    for matrix in [dspin, chea]:
        signed_features.extend(
            [
                np.maximum(matrix, 0.0),
                np.maximum(-matrix, 0.0),
                np.maximum(matrix.T, 0.0),
                np.maximum(-matrix.T, 0.0),
            ]
        )
    sensitivity = _randomized_embedding(np.concatenate(signed_features, axis=1), seed=31)
    return unsigned, sensitivity


def _term_score(values: np.ndarray, member_indices: np.ndarray, direction: float = 0.0) -> float:
    if member_indices.size == 0:
        return 0.0
    selected = values[member_indices]
    if direction:
        selected = selected * direction
    else:
        selected = np.abs(selected)
    return float(selected.sum() / math.sqrt(member_indices.size))


def _degree_bins(degrees: np.ndarray, n_bins: int = 5) -> np.ndarray:
    if degrees.size == 0 or np.allclose(degrees, degrees[0]):
        return np.zeros(degrees.size, dtype=int)
    quantiles = np.unique(np.quantile(degrees, np.linspace(0, 1, n_bins + 1)))
    if quantiles.size <= 2:
        return np.zeros(degrees.size, dtype=int)
    return np.digitize(degrees, quantiles[1:-1], right=True)


def _bh_adjust(values: list[float]) -> list[float]:
    if not values:
        return []
    array = np.asarray(values, dtype=float)
    order = np.argsort(array)
    ranked = array[order] * len(array) / np.arange(1, len(array) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adjusted = np.empty_like(ranked)
    adjusted[order] = np.clip(ranked, 0.0, 1.0)
    return adjusted.tolist()


def _degree_matched_null_scores(
    observed: float,
    metric: np.ndarray,
    member_indices: np.ndarray,
    bins: np.ndarray,
    permutations: int,
    rng: np.random.Generator,
) -> np.ndarray:
    if member_indices.size == 0 or permutations <= 0:
        return np.empty(0, dtype=float)
    pools = {bin_id: np.flatnonzero(bins == bin_id) for bin_id in np.unique(bins)}
    null = np.zeros(permutations, dtype=float)
    for bin_id in np.unique(bins[member_indices]):
        count = int(np.count_nonzero(bins[member_indices] == bin_id))
        pool = pools[int(bin_id)]
        if count <= 0 or pool.size == 0:
            continue
        if count > pool.size:
            sampled = rng.choice(pool, size=(permutations, count), replace=True)
        else:
            # `Generator.choice(..., replace=False)` treats a 2D size as one
            # global draw. Generate one without-replacement sample per null set.
            random_keys = rng.random((permutations, pool.size))
            positions = np.argpartition(random_keys, count - 1, axis=1)[:, :count]
            sampled = pool[positions]
        null += metric[sampled].sum(axis=1)
    return null / math.sqrt(member_indices.size)


def _empirical_pvalue(
    observed: float,
    metric: np.ndarray,
    member_indices: np.ndarray,
    bins: np.ndarray,
    permutations: int,
    rng: np.random.Generator,
    adaptive: bool = False,
    max_permutations: int = ADAPTIVE_PVALUE_MAX_PERMUTATIONS,
    batch_permutations: int = ADAPTIVE_PVALUE_BATCH_PERMUTATIONS,
) -> dict[str, Any]:
    if member_indices.size == 0 or permutations <= 0:
        return {
            "pvalue": float("nan"),
            "effect_z": float("nan"),
            "permutations": 0,
            "exceedances": 0,
            "pvalue_is_upper_bound": False,
        }

    total = 0
    total_sum = 0.0
    total_sum_squares = 0.0
    exceedances = 0

    def add_batch(draws: int) -> None:
        nonlocal total, total_sum, total_sum_squares, exceedances
        null = _degree_matched_null_scores(
            observed,
            metric,
            member_indices,
            bins,
            draws,
            rng,
        )
        total += int(null.size)
        total_sum += float(null.sum())
        total_sum_squares += float(np.dot(null, null))
        exceedances += int(np.count_nonzero(null >= observed))

    add_batch(permutations)
    ceiling = max(permutations, max_permutations)
    if adaptive:
        while exceedances == 0 and total < ceiling:
            add_batch(min(batch_permutations, ceiling - total))

    variance = (
        max(0.0, (total_sum_squares - (total_sum * total_sum / total)) / (total - 1))
        if total > 1
        else 0.0
    )
    null_std = math.sqrt(variance)
    effect_z = (observed - total_sum / total) / null_std if null_std > 0 else 0.0
    return {
        "pvalue": float((1 + exceedances) / (total + 1)),
        "effect_z": float(effect_z),
        "permutations": total,
        "exceedances": exceedances,
        "pvalue_is_upper_bound": bool(exceedances == 0),
    }


def _validate_term_rows(
    rows: list[dict[str, Any]],
    gene_values: np.ndarray,
    bins: np.ndarray,
    permutations: int,
) -> None:
    by_library: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        if row["available_gene_count"] > 0:
            by_library[row["library"]].append(row)
    for library, library_rows in by_library.items():
        selected = sorted(library_rows, key=lambda row: row["direct_score"], reverse=True)[:200]
        seed = int(hashlib.sha256(library.encode()).hexdigest()[:8], 16)
        rng = np.random.default_rng(seed)
        tested: list[dict[str, Any]] = []
        pvalues: list[float] = []
        for row in selected:
            direction = float(row["_direction"])
            metric = gene_values * direction if direction else np.abs(gene_values)
            result = _empirical_pvalue(
                float(row["direct_score"]),
                metric,
                np.asarray(row["_member_indices"], dtype=int),
                bins,
                permutations,
                rng,
                adaptive=row["role"] in ADAPTIVE_PVALUE_ROLES,
            )
            row["direct_empirical_p"] = result["pvalue"]
            row["direct_effect_z"] = result["effect_z"]
            row["direct_permutations"] = result["permutations"]
            row["direct_exceedance_count"] = result["exceedances"]
            row["direct_pvalue_is_upper_bound"] = result["pvalue_is_upper_bound"]
            tested.append(row)
            pvalues.append(result["pvalue"])
        for row, adjusted in zip(tested, _bh_adjust(pvalues)):
            row["direct_bh_fdr"] = adjusted


def _score_term_sets(
    term_sets: list[dict[str, Any]],
    gene_lookup: dict[str, int],
    seed: np.ndarray,
    dspin_states: list[np.ndarray],
    chea_states: list[np.ndarray],
    manifold_similarity: np.ndarray,
    sensitivity_similarity: np.ndarray,
    ppi_proximity: np.ndarray,
    coexpression_proximity: np.ndarray,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for term in term_sets:
        members = sorted(
            {
                gene_lookup[_canonical_gene(gene)]
                for gene in term["genes"]
                if _canonical_gene(gene) in gene_lookup
            }
        )
        indices = np.asarray(members, dtype=int)
        direction = float(term.get("direction", 0.0))
        dspin_depths = [_term_score(state, indices, direction) for state in dspin_states]
        chea_depths = [_term_score(state, indices, direction) for state in chea_states]
        leading = sorted(indices.tolist(), key=lambda index: abs(seed[index]), reverse=True)[:25]
        rows.append(
            {
                "species": "human",
                "role": term["role"],
                "library": term["library"],
                "term_id": term["term_id"],
                "term_label": term["term_label"],
                "term_size": int(term["term_size"]),
                "available_gene_count": int(indices.size),
                "direct_score": _term_score(seed, indices, direction),
                "direct_effect_z": float("nan"),
                "direct_empirical_p": float("nan"),
                "direct_bh_fdr": float("nan"),
                "direct_permutations": 0,
                "direct_exceedance_count": 0,
                "direct_pvalue_is_upper_bound": False,
                "dspin_depth_0_score": dspin_depths[0],
                "dspin_depth_1_score": dspin_depths[1],
                "dspin_depth_2_score": dspin_depths[2],
                "dspin_depth_3_score": dspin_depths[3],
                "dspin_score": float(np.mean(dspin_depths[1:])),
                "chea_depth_0_score": chea_depths[0],
                "chea_depth_1_score": chea_depths[1],
                "chea_depth_2_score": chea_depths[2],
                "chea_depth_3_score": chea_depths[3],
                "chea_score": float(np.mean(chea_depths[1:])),
                "manifold_score": float(manifold_similarity[indices].mean()) if indices.size else 0.0,
                "signed_sensitivity_score": float(sensitivity_similarity[indices].mean()) if indices.size else 0.0,
                "ppi_score": _term_score(ppi_proximity, indices, 0.0),
                "coexpression_score": _term_score(coexpression_proximity, indices, 0.0),
                "_member_indices": indices.tolist(),
                "_leading_indices": leading,
                "_direction": direction,
            }
        )
    return rows


def _parse_mondo_obo(text: str) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    current: dict[str, Any] = {}
    for line in text.splitlines() + ["[Term]"]:
        if line == "[Term]":
            if current.get("id") and current.get("name") and not current.get("obsolete"):
                records.append(current)
            current = {"synonyms": [], "parent_ids": [], "obsolete": False}
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
            current["parent_ids"].append(line.split("!", 1)[0].split("is_a:", 1)[1].strip())
        elif line == "is_obsolete: true":
            current["obsolete"] = True
    by_id = {record["id"]: record for record in records}
    for record in records:
        parents = [by_id[parent]["name"] for parent in record["parent_ids"] if parent in by_id]
        record["family"] = parents[0] if parents else record["name"]
        record.pop("obsolete", None)
    return records


def _load_mondo(cache_root: Path, manifest: dict[str, Any]) -> list[dict[str, Any]]:
    raw_path = manifest.get("mondo", {}).get("terms_path", MONDO_TERMS_FILE)
    path = Path(raw_path)
    if not path.is_absolute():
        path = cache_root / path
    payload = _read_json(path, {})
    terms = payload.get("terms", [])
    return [term for term in terms if isinstance(term, dict) and term.get("id") and term.get("name")]


def _mondo_index(terms: list[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    index: dict[str, dict[str, Any]] = {}
    for term in terms:
        if not public_evidence.is_specific_disease_term(str(term["id"]), str(term["name"])):
            continue
        for label in [term["name"], *term.get("synonyms", [])]:
            normalized = _canonical_label(label)
            if normalized and normalized not in index:
                index[normalized] = term
    return index


def _lineage_names(term: dict[str, Any], by_id: dict[str, dict[str, Any]]) -> list[str]:
    names: list[str] = []
    seen: set[str] = set()
    pending = list(term.get("parent_ids", []))
    while pending:
        term_id = pending.pop()
        if term_id in seen:
            continue
        seen.add(term_id)
        parent = by_id.get(term_id)
        if not parent:
            continue
        names.append(str(parent.get("name", "")))
        pending.extend(parent.get("parent_ids", []))
    return names


def classify_hypothesis(term: dict[str, Any]) -> str:
    """Ontology/name-based classification; no disease-gene list is used."""
    by_id = term.get("_by_id", {})
    text = " ".join([str(term.get("name", "")), str(term.get("family", "")), *_lineage_names(term, by_id)]).casefold()
    if re.search(r"\b(cancer|carcinoma|leukemia|lymphoma|tumou?r|neoplasm|malignan)", text):
        return "Cancer-associated mimic"
    if re.search(r"\b(sepsis|infection|infectious|viral|bacterial|fungal|covid|pneumonia)", text):
        return "Infection/sepsis-like mimic"
    if re.search(r"\b(autoimmune|inflammatory|arthritis|lupus|vasculitis|dermatitis|inflammatory bowel)", text):
        return "Acquired inflammatory or autoimmune mimic"
    if re.search(
        r"\b(monogenic|mendelian|single[- ]gene|inborn error|primary immunodeficiency|"
        r"inborn errors of immunity|congenital|genetic disease)\b",
        text,
    ):
        return "Monogenic/genetic disease"
    return "Other disease evidence"


def _acat(pvalues: list[float]) -> float:
    usable = np.clip(np.asarray(pvalues, dtype=float), 1e-15, 1.0 - 1e-15)
    if usable.size == 0:
        return 1.0
    statistic = float(np.tan((0.5 - usable) * math.pi).mean())
    return float(np.clip(0.5 - math.atan(statistic) / math.pi, 0.0, 1.0))


def _join_unique(values: list[str], limit: int = 25) -> str:
    unique = list(dict.fromkeys(value for value in values if value))
    return ";".join(unique[:limit])


def _build_disease_and_mimic_tables(
    supporting_terms: pd.DataFrame,
    mondo_terms: list[dict[str, Any]],
    network_primary_promoted: bool,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, list[str]]:
    if supporting_terms.empty:
        return _empty(DISEASE_COLUMNS), _empty(MIMIC_COLUMNS), supporting_terms, []
    warnings: list[str] = []
    by_id = {str(term["id"]): term for term in mondo_terms}
    index = _mondo_index(mondo_terms)
    scored = supporting_terms.copy()
    scored["mondo_id"] = ""
    scored["mondo_name"] = ""
    scored["disease_family"] = ""
    scored["hypothesis_class"] = ""
    direct = scored.loc[scored["role"] == "direct_disease"].copy()
    unmapped = 0
    for row_index, row in direct.iterrows():
        term = index.get(_canonical_label(row["term_label"]))
        if term is None:
            unmapped += 1
            continue
        annotated = dict(term)
        annotated["_by_id"] = by_id
        scored.loc[row_index, "mondo_id"] = str(term["id"])
        scored.loc[row_index, "mondo_name"] = str(term["name"])
        scored.loc[row_index, "disease_family"] = str(term.get("family", term["name"]))
        scored.loc[row_index, "hypothesis_class"] = classify_hypothesis(annotated)
    if unmapped:
        warnings.append(f"{unmapped} direct disease terms lacked an exact MONDO mapping and remain supporting evidence only.")

    rows: list[dict[str, Any]] = []
    valid = scored.loc[
        (scored["role"] == "direct_disease")
        & scored["mondo_id"].astype(str).str.len().gt(0)
        & scored["direct_empirical_p"].notna()
    ]
    for disease_id, values in valid.groupby("mondo_id", sort=False):
        per_resource = values.sort_values("direct_empirical_p").groupby("library", sort=False).first().reset_index()
        pvalue = _acat(per_resource["direct_empirical_p"].astype(float).tolist())
        matched = []
        for leading in per_resource["leading_genes"].astype(str):
            matched.extend(gene for gene in leading.split(";") if gene)
        first = per_resource.iloc[0]
        rows.append(
            {
                "disease_family": first["disease_family"],
                "disease_id": disease_id,
                "disease_name": first["mondo_name"],
                "hypothesis_class": first["hypothesis_class"],
                "combined_p_value": pvalue,
                "combined_fdr": float("nan"),
                "evidence_score": -math.log10(max(pvalue, 1e-300)),
                "evidence_tier": "corroborated_direct_evidence" if int(per_resource["library"].nunique()) >= 2 else "single_source_direct_evidence",
                "direct_effect_z": float(per_resource["direct_effect_z"].mean()),
                "direct_permutations": _join_unique(
                    [f"{row.library}:{int(row.direct_permutations)}" for row in per_resource.itertuples()]
                ),
                "combined_p_value_is_upper_bound": bool(
                    per_resource["direct_pvalue_is_upper_bound"].fillna(False).astype(bool).any()
                ),
                "direct_resource_count": int(per_resource["library"].nunique()),
                "supporting_sources": _join_unique(per_resource["library"].astype(str).tolist(), limit=10),
                "matched_genes": _join_unique(matched),
                "dspin_score": float(per_resource["dspin_score"].mean()),
                "chea_score": float(per_resource["chea_score"].mean()),
                "manifold_score": float(per_resource["manifold_score"].mean()),
                "ppi_score": float(per_resource["ppi_score"].mean()),
                "coexpression_score": float(per_resource["coexpression_score"].mean()),
                "network_primary_promoted": bool(network_primary_promoted),
            }
        )
    classified = pd.DataFrame(rows)
    if not classified.empty:
        classified["combined_fdr"] = _bh_adjust(classified["combined_p_value"].astype(float).tolist())
    genetic = classified.loc[~classified["hypothesis_class"].str.contains("mimic", case=False, na=False)].copy()
    named_mimics = classified.loc[classified["hypothesis_class"].str.contains("mimic", case=False, na=False)].copy()

    context = scored.loc[
        (scored["role"] == "cancer_context") & scored["direct_empirical_p"].notna()
    ].copy()
    context_rows = []
    for _, row in context.iterrows():
        pvalue = float(row["direct_empirical_p"])
        context_rows.append(
            {
                "disease_family": "Cancer cell-line context",
                "disease_id": "",
                "disease_name": str(row["term_label"]),
                "hypothesis_class": "Cancer-associated mimic",
                "combined_p_value": pvalue,
                "combined_fdr": float("nan"),
                "evidence_score": -math.log10(max(pvalue, 1e-300)),
                "evidence_tier": "single_source_cancer_context",
                "direct_effect_z": float(row["direct_effect_z"]),
                "direct_permutations": f"{row['library']}:{int(row['direct_permutations'])}",
                "combined_p_value_is_upper_bound": bool(row["direct_pvalue_is_upper_bound"]),
                "direct_resource_count": 1,
                "supporting_sources": str(row["library"]),
                "matched_genes": str(row["leading_genes"]),
                "dspin_score": float(row["dspin_score"]),
                "chea_score": float(row["chea_score"]),
                "manifold_score": float(row["manifold_score"]),
                "ppi_score": float(row["ppi_score"]),
                "coexpression_score": float(row["coexpression_score"]),
            }
        )
    for row in context_rows:
        row["network_primary_promoted"] = bool(network_primary_promoted)
    mimic = pd.concat([named_mimics, pd.DataFrame(context_rows)], ignore_index=True)
    if not mimic.empty:
        mimic["combined_fdr"] = _bh_adjust(mimic["combined_p_value"].astype(float).tolist())
        if network_primary_promoted:
            mimic["_direct_rank"] = mimic["combined_p_value"].rank(method="average", ascending=True)
            mimic["_directed_rank"] = (
                (mimic["dspin_score"] + mimic["chea_score"]) / 2.0
            ).rank(method="average", ascending=False)
            mimic["_manifold_rank"] = mimic["manifold_score"].rank(method="average", ascending=False)
            mimic["_primary_rank"] = (
                mimic["_direct_rank"] + mimic["_directed_rank"] + mimic["_manifold_rank"]
            ) / 3.0
            mimic = mimic.sort_values(
                ["_primary_rank", "combined_p_value", "direct_resource_count", "direct_effect_z"],
                ascending=[True, True, False, False],
            )
        else:
            mimic = mimic.sort_values(
                ["combined_p_value", "direct_resource_count", "direct_effect_z"],
                ascending=[True, False, False],
            )
        mimic = mimic.drop(columns=["_direct_rank", "_directed_rank", "_manifold_rank", "_primary_rank"], errors="ignore").reset_index(drop=True)
        mimic.insert(0, "rank", np.arange(1, len(mimic) + 1))
    if not genetic.empty:
        if network_primary_promoted:
            genetic["_direct_rank"] = genetic["combined_p_value"].rank(method="average", ascending=True)
            genetic["_directed_rank"] = (
                (genetic["dspin_score"] + genetic["chea_score"]) / 2.0
            ).rank(method="average", ascending=False)
            genetic["_manifold_rank"] = genetic["manifold_score"].rank(method="average", ascending=False)
            genetic["_primary_rank"] = (
                genetic["_direct_rank"] + genetic["_directed_rank"] + genetic["_manifold_rank"]
            ) / 3.0
            genetic = genetic.sort_values(
                ["_primary_rank", "combined_p_value", "direct_resource_count", "direct_effect_z"],
                ascending=[True, True, False, False],
            )
        else:
            genetic = genetic.sort_values(
                ["combined_p_value", "direct_resource_count", "direct_effect_z"],
                ascending=[True, False, False],
            )
        genetic = genetic.drop(columns=["_direct_rank", "_directed_rank", "_manifold_rank", "_primary_rank"], errors="ignore").reset_index(drop=True)
        genetic.insert(0, "rank", np.arange(1, len(genetic) + 1))
    return genetic.reindex(columns=DISEASE_COLUMNS), mimic.reindex(columns=MIMIC_COLUMNS), scored, warnings


def _build_evidence_edges(supporting_terms: pd.DataFrame, genes: list[str], seed: np.ndarray) -> pd.DataFrame:
    if supporting_terms.empty:
        return _empty(EDGE_COLUMNS)
    candidates = supporting_terms.copy()
    candidates["_priority"] = candidates["direct_empirical_p"].fillna(1.0)
    candidates = candidates.sort_values(["_priority", "direct_score"], ascending=[True, False]).head(100)
    rows = []
    for _, row in candidates.iterrows():
        for index in row.get("_leading_indices", [])[:25]:
            rows.append(
                {
                    "species": "human",
                    "evidence_kind": "kg_term_gene",
                    "role": row["role"],
                    "resource": row["library"],
                    "term_id": row["term_id"],
                    "term_label": row["term_label"],
                    "disease_id": row.get("mondo_id", ""),
                    "disease_name": row.get("mondo_name", ""),
                    "patient_gene": genes[int(index)],
                    "patient_signed_score": float(seed[int(index)]),
                    "direct_score": float(row["direct_score"]),
                    "dspin_score": float(row["dspin_score"]),
                    "chea_score": float(row["chea_score"]),
                    "manifold_score": float(row["manifold_score"]),
                    "ppi_score": float(row["ppi_score"]),
                    "coexpression_score": float(row["coexpression_score"]),
                }
            )
    return pd.DataFrame(rows, columns=EDGE_COLUMNS)


def _network_concordance(
    genes: list[str],
    dspin: np.ndarray,
    chea: np.ndarray,
    chea_variants: dict[str, np.ndarray],
    seed: np.ndarray,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    nonzero = np.argwhere(np.abs(chea) > 0)
    ranked = sorted(nonzero.tolist(), key=lambda pair: abs(chea[pair[0], pair[1]]) * (abs(seed[pair[0]]) + abs(seed[pair[1]])), reverse=True)
    for target, source in ranked[:250]:
        dspin_weight = float(dspin[target, source])
        chea_weight = float(chea[target, source])
        rows.append(
            {
                "row_type": "patient_edge",
                "network_variant": "node_draw",
                "source_gene": genes[source],
                "target_gene": genes[target],
                "dspin_weight": dspin_weight,
                "chea_weight": chea_weight,
                "source_signal": float(seed[source]),
                "target_signal": float(seed[target]),
                "sign_concordant": bool(dspin_weight * chea_weight > 0),
                "overlapping_edges": np.nan,
                "topology_correlation": np.nan,
            }
        )
    primary_mask = np.abs(chea) > 0
    for name, variant in chea_variants.items():
        variant_mask = np.abs(variant) > 0
        overlap = primary_mask & variant_mask
        if overlap.any():
            first = np.abs(chea[overlap])
            second = np.abs(variant[overlap])
            correlation = (
                float(np.corrcoef(first, second)[0, 1])
                if first.size > 1 and float(first.std()) > 0 and float(second.std()) > 0
                else float("nan")
            )
        else:
            correlation = float("nan")
        rows.append(
            {
                "row_type": "variant_qc",
                "network_variant": name,
                "source_gene": "",
                "target_gene": "",
                "dspin_weight": np.nan,
                "chea_weight": np.nan,
                "source_signal": np.nan,
                "target_signal": np.nan,
                "sign_concordant": np.nan,
                "overlapping_edges": int(overlap.sum()),
                "topology_correlation": correlation,
            }
        )
    return pd.DataFrame(rows, columns=CONCORDANCE_COLUMNS)


def _manifold_modules(
    genes: list[str],
    similarity: np.ndarray,
    sensitivity: np.ndarray,
    views: dict[str, np.ndarray],
) -> pd.DataFrame:
    if similarity.size == 0 or not np.any(similarity):
        return _empty(MODULE_COLUMNS)
    selected = np.argsort(-similarity)[: min(80, len(genes))]
    strength = np.zeros((len(genes), len(genes)), dtype=float)
    for matrix in views.values():
        strength += np.maximum(np.abs(matrix), np.abs(matrix.T))
    submatrix = strength[np.ix_(selected, selected)]
    positive = submatrix[submatrix > 0]
    threshold = float(np.quantile(positive, 0.75)) if positive.size else float("inf")
    adjacency = {int(index): set() for index in selected}
    for left, index in enumerate(selected):
        neighbours = selected[np.flatnonzero(submatrix[left] >= threshold)] if math.isfinite(threshold) else []
        for neighbour in neighbours:
            if index != neighbour:
                adjacency[int(index)].add(int(neighbour))
    components: list[list[int]] = []
    remaining = set(adjacency)
    while remaining:
        start = remaining.pop()
        component = {start}
        pending = [start]
        while pending:
            current = pending.pop()
            newly = adjacency[current] & remaining
            pending.extend(newly)
            component.update(newly)
            remaining -= newly
        components.append(sorted(component, key=lambda index: similarity[index], reverse=True))
    rows = []
    for position, component in enumerate(sorted(components, key=lambda group: max(similarity[index] for index in group), reverse=True), start=1):
        active = []
        for name, matrix in views.items():
            if len(component) > 1 and np.any(matrix[np.ix_(component, component)]):
                active.append(name)
        rows.append(
            {
                "module_id": f"M{position}",
                "gene_count": len(component),
                "genes": ";".join(genes[index] for index in component[:30]),
                "mean_manifold_similarity": float(np.mean(similarity[component])),
                "mean_signed_sensitivity_similarity": float(np.mean(sensitivity[component])) if sensitivity.size else 0.0,
                "supporting_views": ";".join(active),
            }
        )
    return pd.DataFrame(rows, columns=MODULE_COLUMNS)


def bootstrap_mrr_improvement(
    direct_ranks: list[int],
    hybrid_ranks: list[int],
    bootstrap_iterations: int = 2_000,
    seed: int = 17,
) -> dict[str, Any]:
    """Conservative promotion gate used by leave-one-resource-out validation."""
    if len(direct_ranks) != len(hybrid_ranks) or not direct_ranks:
        return {"status": "insufficient_cases", "approved": False}
    direct = 1.0 / np.asarray(direct_ranks, dtype=float)
    hybrid = 1.0 / np.asarray(hybrid_ranks, dtype=float)
    rng = np.random.default_rng(seed)
    indices = rng.integers(0, len(direct), size=(bootstrap_iterations, len(direct)))
    deltas = (hybrid[indices] - direct[indices]).mean(axis=1)
    lower, upper = np.quantile(deltas, [0.025, 0.975])
    return {
        "status": "complete",
        "direct_mrr": float(direct.mean()),
        "hybrid_mrr": float(hybrid.mean()),
        "mrr_improvement": float(hybrid.mean() - direct.mean()),
        "bootstrap_ci_95": [float(lower), float(upper)],
        "approved": bool(lower > 0),
    }


def _member_index_set(value: Any) -> set[int]:
    if isinstance(value, (list, tuple, np.ndarray)):
        return {int(index) for index in value}
    return set()


def _descending_rank(values: dict[str, float], target: str) -> int:
    return _descending_rank_map(values)[target]


def _ascending_rank(values: dict[str, float], target: str) -> int:
    return _ascending_rank_map(values)[target]


def _descending_rank_map(values: dict[str, float]) -> dict[str, int]:
    return {key: rank for rank, key in enumerate(sorted(values, key=lambda key: (-values[key], key)), start=1)}


def _ascending_rank_map(values: dict[str, float]) -> dict[str, int]:
    return {key: rank for rank, key in enumerate(sorted(values, key=lambda key: (values[key], key)), start=1)}


def validate_leave_one_resource_out(
    annotated_terms: pd.DataFrame,
    dspin: np.ndarray,
    chea: np.ndarray,
    manifold_embedding: np.ndarray,
    max_cases_per_resource: int = 250,
    minimum_cases: int = 20,
    bootstrap_iterations: int = 2_000,
) -> dict[str, Any]:
    """Evaluate direct and network retrieval without using a patient query.

    For each direct disease resource, a disease gene set is held out. The same
    MONDO concept must be represented by one or more other direct resources to
    become a retrieval case. Direct overlap, directed DSPIN/ChEA propagation,
    unsigned-manifold proximity, and equal-rank hybrid retrieval are compared.
    """
    required = {"role", "library", "mondo_id", "_member_indices"}
    if annotated_terms.empty or not required.issubset(annotated_terms.columns):
        return {"status": "insufficient_cases", "approved": False, "case_count": 0}
    terms = annotated_terms.loc[
        (annotated_terms["role"] == "direct_disease")
        & annotated_terms["mondo_id"].astype(str).str.len().gt(0)
    ].copy()
    terms["_members"] = terms["_member_indices"].apply(_member_index_set)
    terms = terms.loc[terms["_members"].map(bool)]
    libraries = sorted(terms["library"].astype(str).unique())
    if len(libraries) < 2:
        return {"status": "insufficient_cases", "approved": False, "case_count": 0}

    cases: list[tuple[str, set[int], dict[str, set[int]]]] = []
    case_counts: dict[str, int] = {}
    for held_out in libraries:
        held = terms.loc[terms["library"] == held_out]
        candidate_rows = terms.loc[terms["library"] != held_out]
        candidates: dict[str, set[int]] = {}
        for disease_id, group in candidate_rows.groupby("mondo_id", sort=False):
            members: set[int] = set()
            for value in group["_members"]:
                members.update(value)
            if members:
                candidates[str(disease_id)] = members
        eligible = []
        for disease_id, group in held.groupby("mondo_id", sort=True):
            target = str(disease_id)
            if target not in candidates or len(candidates) < 2:
                continue
            members: set[int] = set()
            for value in group["_members"]:
                members.update(value)
            if members:
                eligible.append((target, members))
        eligible = eligible[:max_cases_per_resource]
        case_counts[held_out] = len(eligible)
        for target, query_members in eligible:
            cases.append((target, query_members, candidates))

    if len(cases) < minimum_cases:
        return {
            "status": "insufficient_cases",
            "approved": False,
            "case_count": len(cases),
            "held_out_resource_case_counts": case_counts,
            "minimum_cases_required": minimum_cases,
        }

    queries = np.zeros((dspin.shape[0], len(cases)), dtype=float)
    for column, (_, members, _) in enumerate(cases):
        queries[list(members), column] = 1.0

    def batch_propagation_signal(matrix: np.ndarray) -> np.ndarray:
        values = np.asarray(matrix, dtype=float).copy()
        np.fill_diagonal(values, 0.0)
        totals = np.abs(values).sum(axis=1, keepdims=True)
        transition = np.divide(values, totals, out=np.zeros_like(values), where=totals > 0)
        state = queries
        accumulated = np.zeros_like(queries)
        for _ in range(3):
            state = transition.dot(state)
            accumulated += np.abs(state)
        return accumulated / 3.0

    dspin_signals = batch_propagation_signal(dspin)
    chea_signals = batch_propagation_signal(chea)
    manifold_signals = np.column_stack(
        [
            _query_embedding_similarity(manifold_embedding, queries[:, column])
            for column in range(queries.shape[1])
        ]
    )

    direct_ranks: list[int] = []
    directed_ranks: list[int] = []
    manifold_ranks: list[int] = []
    hybrid_ranks: list[int] = []
    for column, (target, query_members, candidates) in enumerate(cases):
        directed_signal = 0.5 * (dspin_signals[:, column] + chea_signals[:, column])
        manifold_signal = manifold_signals[:, column]
        direct_scores: dict[str, float] = {}
        directed_scores: dict[str, float] = {}
        manifold_scores: dict[str, float] = {}
        for candidate_id, candidate_members in candidates.items():
            indices = np.fromiter(candidate_members, dtype=int)
            direct_scores[candidate_id] = len(query_members & candidate_members) / math.sqrt(
                len(query_members) * len(candidate_members)
            )
            directed_scores[candidate_id] = _term_score(directed_signal, indices)
            manifold_scores[candidate_id] = float(manifold_signal[indices].mean())
        direct_rank_map = _descending_rank_map(direct_scores)
        directed_rank_map = _descending_rank_map(directed_scores)
        manifold_rank_map = _descending_rank_map(manifold_scores)
        direct_ranks.append(direct_rank_map[target])
        directed_ranks.append(directed_rank_map[target])
        manifold_ranks.append(manifold_rank_map[target])
        hybrid = {
            candidate: (direct_rank_map[candidate] + directed_rank_map[candidate] + manifold_rank_map[candidate]) / 3.0
            for candidate in candidates
        }
        hybrid_ranks.append(_ascending_rank_map(hybrid)[target])

    promotion = bootstrap_mrr_improvement(
        direct_ranks,
        hybrid_ranks,
        bootstrap_iterations=bootstrap_iterations,
    )
    return {
        **promotion,
        "case_count": len(direct_ranks),
        "held_out_resource_case_counts": case_counts,
        "retrieval_modes": {
            "direct_mrr": float(np.mean(1.0 / np.asarray(direct_ranks, dtype=float))),
            "directed_mrr": float(np.mean(1.0 / np.asarray(directed_ranks, dtype=float))),
            "unsigned_manifold_mrr": float(np.mean(1.0 / np.asarray(manifold_ranks, dtype=float))),
            "hybrid_mrr": float(np.mean(1.0 / np.asarray(hybrid_ranks, dtype=float))),
        },
        "hybrid": "equal-rank aggregation of direct, directed DSPIN/ChEA, and unsigned-manifold retrieval",
        "max_cases_per_resource": max_cases_per_resource,
        "minimum_cases_required": minimum_cases,
    }


def build_human_perturbseqr_mechanisms(
    evidence_cache: Path,
    signed_gene_scores: pd.Series,
    signed_program_scores: pd.Series,
    gene_network: pd.DataFrame | None,
    gene_to_program: pd.DataFrame | None,
    cell_type_table: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, list[str], dict[str, Any]]:
    """Preserve the existing human Perturb-Seqr algorithm without mouse work."""
    records, cache_manifest, warnings = public_evidence.read_evidence_cache(evidence_cache)
    patient_cell_types = set(cell_type_table.get("cell_type", pd.Series(dtype=str)).astype(str).str.casefold())
    edges: list[dict[str, Any]] = []
    rejected = 0
    for raw_record in records:
        record, reason = public_evidence.normalize_evidence_record(raw_record)
        if record is None:
            rejected += 1
            continue
        if record.get("species") != "human" or record.get("resource") != "perturbseqr":
            continue
        if record.get("kind") != "gene_perturbation":
            rejected += 1
            continue
        score = public_evidence.score_record(
            record,
            signed_gene_scores,
            signed_program_scores,
            gene_network,
            gene_to_program,
            patient_cell_types,
        )
        role = "mechanistic_phenocopy" if score["evidence_score"] >= score["counterevidence_score"] else "mechanistic_reversal"
        edges.append(public_evidence._edge_row(record, score, role))
    edge_table = pd.DataFrame(edges, columns=public_evidence.EDGE_COLUMNS)
    table = public_evidence.aggregate_mechanisms(edge_table, "human")
    raw_source = cache_manifest.get("sources", {}).get("perturbseqr", {})
    source_fields = [
        "source_url",
        "source_version",
        "retrieved_at",
        "status",
        "record_count",
        "records_sha256",
        "chemical_data",
        "query_mode",
    ]
    source_manifest = {field: raw_source[field] for field in source_fields if field in raw_source}
    source_manifest["species"] = ["human"]
    metadata = {
        "cache_root": str(evidence_cache),
        "cache_schema_version": cache_manifest.get("schema_version"),
        "human_perturbseqr_source": source_manifest,
        "human_perturbseqr_records_scored": int(len(edge_table)),
        "rejected_record_count": rejected,
    }
    return table, edge_table, warnings, metadata


def build_human_kg_hypotheses(
    kg_cache: Path,
    evidence_cache: Path,
    signed_gene_scores: pd.Series,
    signed_program_scores: pd.Series,
    gene_network: pd.DataFrame | None,
    gene_to_program: pd.DataFrame | None,
    cell_type_table: pd.DataFrame,
    permutations: int = 2_000,
) -> dict[str, Any]:
    if permutations <= 0:
        raise ValueError("--kg-permutations must be positive.")
    manifest, warnings = read_kg_cache(kg_cache)
    mechanisms, mechanism_edges, mechanism_warnings, mechanism_manifest = build_human_perturbseqr_mechanisms(
        evidence_cache,
        signed_gene_scores,
        signed_program_scores,
        gene_network,
        gene_to_program,
        cell_type_table,
    )
    warnings.extend(mechanism_warnings)
    genes = signed_gene_scores.index.astype(str).tolist()
    if not genes or not manifest:
        empty_terms = _empty(TERM_COLUMNS)
        empty_edges = _empty(EDGE_COLUMNS)
        combined_edges = pd.concat([empty_edges, mechanism_edges.assign(evidence_kind="human_perturbseqr")], ignore_index=True, sort=False)
        return {
            "human_diseases": _empty(DISEASE_COLUMNS),
            "human_mimics": _empty(MIMIC_COLUMNS),
            "human_mechanisms": mechanisms,
            "supporting_terms": empty_terms,
            "network_concordance": _empty(CONCORDANCE_COLUMNS),
            "manifold_modules": _empty(MODULE_COLUMNS),
            "edges": combined_edges,
            "manifest": {
                "engine": "human_kg_multiview",
                "species": "human",
                "kg_cache": str(kg_cache),
                "kg_cache_manifest": manifest,
                "human_perturbseqr": mechanism_manifest,
                "network_primary_validation": {"status": "not_run", "approved": False},
            },
            "warnings": list(dict.fromkeys(warnings)),
        }

    seed = signed_gene_scores.reindex(genes).fillna(0.0).astype(float).to_numpy()
    maximum = float(np.max(np.abs(seed))) if seed.size else 0.0
    seed = seed / maximum if maximum > 0 else seed
    gene_lookup = {_canonical_gene(gene): index for index, gene in enumerate(genes)}
    dspin = _aligned_matrix(gene_network, genes)

    def cached(asset_id: str) -> Path | None:
        return _asset_path(kg_cache, manifest, asset_id)

    chea_path = cached("chea_node_draw")
    chea = _load_chea_matrix(chea_path, genes) if chea_path else np.zeros_like(dspin)
    if chea_path is None:
        warnings.append("Primary ChEA-KG node-draw GRN is unavailable; ChEA support is zero.")
    chea_variants: dict[str, np.ndarray] = {}
    for asset_id in ["chea_target_set_swap", "chea_unfiltered"]:
        path = cached(asset_id)
        if path:
            chea_variants[asset_id.removeprefix("chea_")] = _load_chea_matrix(path, genes)

    ppi_path = cached("string_ppi")
    ppi = _load_ppi_matrix(ppi_path, genes) if ppi_path else np.zeros_like(dspin)
    if ppi_path is None:
        warnings.append("STRING PPI is unavailable; unsigned PPI support is zero.")

    term_sets: list[dict[str, Any]] = []
    cache_stats: dict[str, Any] = {}
    archs4_terms: list[dict[str, Any]] = []
    for asset_id, spec in ASSET_SPECS.items():
        if spec.get("kind") or asset_id == "string_ppi":
            continue
        path = cached(asset_id)
        if path is None:
            cache_stats[asset_id] = {"status": "missing"}
            continue
        try:
            terms, stats = load_term_sets(
                path,
                asset_id,
                str(spec["role"]),
                pd.Index(genes),
                str(spec["species_policy"]),
                direction=float(spec.get("direction", 0.0)),
            )
        except (OSError, ValueError, pd.errors.ParserError) as exc:
            warnings.append(f"Could not read {asset_id}: {exc}")
            cache_stats[asset_id] = {"status": "error", "message": str(exc)}
            continue
        cache_stats[asset_id] = {"status": "loaded", **stats, "term_count": len(terms)}
        term_sets.extend(terms)
        if asset_id == "archs4_tf":
            archs4_terms = terms
    coexpression = _archs4_tf_coexpression_matrix(archs4_terms, genes)
    if not archs4_terms:
        warnings.append("ARCHS4 human TF coexpression is unavailable; unsigned coexpression support is zero.")

    dspin_states = propagate_signed(seed, dspin)
    chea_states = propagate_signed(seed, chea)
    manifold_embedding, sensitivity_embedding = _multi_view_embeddings(dspin, chea, ppi, coexpression)
    manifold_similarity = _query_embedding_similarity(manifold_embedding, seed)
    sensitivity_similarity = _query_embedding_similarity(sensitivity_embedding, seed)
    ppi_proximity = _unsigned_proximity(seed, ppi)
    coexpression_proximity = _unsigned_proximity(seed, coexpression)
    rows = _score_term_sets(
        term_sets,
        gene_lookup,
        seed,
        dspin_states,
        chea_states,
        manifold_similarity,
        sensitivity_similarity,
        ppi_proximity,
        coexpression_proximity,
    )
    degrees = (
        np.abs(dspin).sum(axis=0)
        + np.abs(dspin).sum(axis=1)
        + np.abs(chea).sum(axis=0)
        + np.abs(chea).sum(axis=1)
        + ppi.sum(axis=0)
        + coexpression.sum(axis=0)
    )
    _validate_term_rows(rows, seed, _degree_bins(degrees), permutations)
    for row in rows:
        row["leading_genes"] = ";".join(genes[index] for index in row["_leading_indices"])
        row["mondo_id"] = ""
        row["mondo_name"] = ""
        row["disease_family"] = ""
        row["hypothesis_class"] = ""
    supporting = pd.DataFrame(rows)
    if supporting.empty:
        supporting = _empty(TERM_COLUMNS)
    mondo_terms = _load_mondo(kg_cache, manifest)
    if not mondo_terms:
        warnings.append("MONDO ontology cache is unavailable; named disease hypotheses cannot be normalized.")
    # The validation corpus is source-only. It must not depend on the current
    # patient's top terms or values before it can promote network-aware ranking.
    _, _, annotated_supporting, _ = _build_disease_and_mimic_tables(
        supporting,
        mondo_terms,
        False,
    )
    validation = manifest.get("network_primary_validation", {"status": "not_run", "approved": False})
    if not isinstance(validation, dict) or validation.get("status") == "not_run":
        validation = validate_leave_one_resource_out(
            annotated_supporting,
            dspin,
            chea,
            manifold_embedding,
        )
        manifest = dict(manifest)
        manifest["network_primary_validation"] = validation
        _write_kg_manifest(kg_cache, manifest)
    approved = bool(validation.get("approved"))
    diseases, mimics, supporting, disease_warnings = _build_disease_and_mimic_tables(
        annotated_supporting,
        mondo_terms,
        approved,
    )
    warnings.extend(disease_warnings)
    kg_edges = _build_evidence_edges(supporting, genes, seed)
    mechanism_edges = mechanism_edges.copy()
    if not mechanism_edges.empty:
        mechanism_edges["evidence_kind"] = "human_perturbseqr"
    edges = pd.concat([kg_edges, mechanism_edges], ignore_index=True, sort=False)
    concordance = _network_concordance(genes, dspin, chea, chea_variants, seed)
    modules = _manifold_modules(
        genes,
        manifold_similarity,
        sensitivity_similarity,
        {"DSPIN": dspin, "ChEA-KG": chea, "STRING PPI": ppi, "ARCHS4 TF coexpression": coexpression},
    )
    output_manifest = {
        "engine": "human_kg_multiview",
        "species": "human",
        "kg_cache": str(kg_cache),
        "kg_cache_manifest": manifest,
        "asset_load_stats": cache_stats,
        "scoring": {
            "direct_gene_score": "absolute signed DSPIN gene magnitude for unsigned terms",
            "directed_propagation_depths": [0, 1, 2, 3],
            "term_validation": {
                "top_terms_per_library": 200,
                "initial_permutations": permutations,
                "adaptive_max_permutations": ADAPTIVE_PVALUE_MAX_PERMUTATIONS,
                "adaptive_batch_permutations": ADAPTIVE_PVALUE_BATCH_PERMUTATIONS,
                "adaptive_roles": sorted(ADAPTIVE_PVALUE_ROLES),
                "correction": "BH",
            },
            "disease_primary": "equal-weight ACAT of direct disease-library empirical p-values",
            "network_primary_validation": validation,
            "unsigned_networks": ["STRING PPI", "symmetrized ARCHS4 TF coexpression projection"],
            "patient_signature_transmitted_for_kg": False,
        },
        "human_perturbseqr": mechanism_manifest,
    }
    return {
        "human_diseases": diseases,
        "human_mimics": mimics,
        "human_mechanisms": mechanisms,
        "supporting_terms": supporting.drop(columns=["_member_indices", "_leading_indices", "_direction"], errors="ignore").reindex(columns=TERM_COLUMNS),
        "network_concordance": concordance,
        "manifold_modules": modules,
        "edges": edges,
        "manifest": output_manifest,
        "warnings": list(dict.fromkeys(warnings)),
    }
