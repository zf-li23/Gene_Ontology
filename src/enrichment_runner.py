"""Enrichr-based enrichment using gseapy for detected communities.

This module is optional; if gseapy is not installed, functions will raise ImportError.
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, Iterable, List, Optional

import pandas as pd
import warnings
import requests
from pathlib import Path
from typing import Union
from scipy.stats import combine_pvalues

try:
    import gseapy as gp
except ImportError as exc:  # pragma: no cover - optional dependency
    gp = None
    _gseapy_import_error = exc
else:
    _gseapy_import_error = None

LOGGER = logging.getLogger(__name__)


def require_gseapy() -> None:
    if gp is None:
        raise ImportError("gseapy is required for Enrichr-based enrichment") from _gseapy_import_error


def enrich_communities(
    communities: Dict[str, List[str]],
    gene_sets: List[str],
    organism: str,
    cutoff: float = 0.05,
    top_terms: int = 10,
    out_file: Optional[Path] = None,
    out_dir: Optional[Path] = None,
) -> pd.DataFrame:
    """Run Enrichr on each community and return a concatenated DataFrame.

    Columns include community_id, Gene_set, Term, P-value, Adjusted P-value, Odds Ratio, Combined Score, Genes.
    """
    require_gseapy()
    all_rows: List[pd.DataFrame] = []
    # try to use tqdm for progress if available
    try:
        from tqdm import tqdm
        iterator = tqdm(communities.items(), desc="Enrichr communities")
    except Exception:
        iterator = communities.items()

    for cid, genes in iterator:
        gene_list = list(dict.fromkeys(genes))
        if not gene_list:
            continue
        try:
            # suppress known gseapy FutureWarning about DataFrame concatenation
            with warnings.catch_warnings():
                warnings.filterwarnings(
                    "ignore",
                    message=".*DataFrame concatenation with empty or all-NA entries is deprecated.*",
                    category=FutureWarning,
                )
                enr = gp.enrichr(
                    gene_list=gene_list,
                    gene_sets=gene_sets,
                    organism=organism,
                    outdir=None,
                    cutoff=cutoff,
                )
        except Exception as exc:  # pragma: no cover - runtime safety
            LOGGER.warning("Enrichr failed for %s: %s", cid, exc)
            continue
        if enr is None or not hasattr(enr, "results"):
            continue
        df = enr.results.copy()

        # Normalize column names: produce a standardized set of columns
        def _find_col(df, candidates):
            for c in candidates:
                if c in df.columns:
                    return c
            # try lower-case match
            lower_map = {col.lower(): col for col in df.columns}
            for c in candidates:
                key = c.lower()
                if key in lower_map:
                    return lower_map[key]
            return None

        term_col = _find_col(df, ["Term", "term"])
        gene_set_col = _find_col(df, ["Gene_set", "gene_set", "Gene Set"]) or _find_col(df, ["Gene_set"])
        p_col = _find_col(df, ["P-value", "Pvalue", "pvalue", "P-value "])
        adjp_col = _find_col(df, ["Adjusted P-value", "Adjusted Pvalue", "adj_p", "Adjusted P-value "])
        overlap_col = _find_col(df, ["Overlap", "Overlap "])
        genes_col = _find_col(df, ["Genes", "genes"])
        combined_col = _find_col(df, ["Combined Score", "CombinedScore", "Combined Score "])

        renamed = pd.DataFrame()
        renamed["community_id"] = [cid] * len(df)
        renamed["gene_set"] = df[gene_set_col] if gene_set_col in df.columns else (df[gene_set_col] if gene_set_col else None)
        renamed["term"] = df[term_col] if term_col in df.columns else df.iloc[:, 0]
        renamed["p_value"] = df[p_col] if p_col in df.columns else None
        renamed["adjusted_p_value"] = df[adjp_col] if adjp_col in df.columns else None
        renamed["overlap"] = df[overlap_col] if overlap_col in df.columns else None
        renamed["combined_score"] = df[combined_col] if combined_col in df.columns else None
        if genes_col and genes_col in df.columns:
            renamed["genes"] = df[genes_col]
        else:
            # try to reconstruct genes column from 'Overlap Genes' style columns
            if "Overlap" in df.columns and "Genes" in df.columns:
                renamed["genes"] = df.get(genes_col)
            else:
                renamed["genes"] = None

        # add community size for downstream filtering/inspection
        renamed["community_size"] = len(gene_list)

        if top_terms and "adjusted_p_value" in renamed.columns:
            try:
                renamed = renamed.sort_values("adjusted_p_value").groupby("gene_set").head(top_terms)
            except Exception:
                pass

        # write per-community file immediately if requested (robust against interruption)
        if out_dir is not None:
            try:
                out_dir.mkdir(parents=True, exist_ok=True)
                per_path = out_dir / f"{cid}_enrich.csv"
                renamed.to_csv(per_path, index=False)
            except Exception:
                LOGGER.warning("Failed to write per-community file for %s", cid)

        all_rows.append(renamed)
    if not all_rows:
        return pd.DataFrame()
    out_df = pd.concat(all_rows, axis=0, ignore_index=True)
    if out_file is not None:
        out_file.parent.mkdir(parents=True, exist_ok=True)
        out_df.to_csv(out_file, index=False)

    # write individual community CSVs if requested
    if out_dir is not None and not out_df.empty:
        out_dir.mkdir(parents=True, exist_ok=True)
        for cid, group in out_df.groupby("community_id"):
            path = out_dir / f"{cid}_enrich.csv"
            group.to_csv(path, index=False)
    return out_df


def download_enrichr_library(library_name: str, out_path: Union[str, Path], timeout: int = 60) -> Path:
    """Download an Enrichr gene-set library and save as GMT-like file.

    The Enrichr text endpoint returns tab-separated lines: term\tdescription\tcomma-separated-genes.
    We convert that to GMT (term\tdescription\tgene1\tgene2...).
    Returns the path to the saved GMT file.
    """
    require_gseapy()  # ensure requests available in environment; still safe if not
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    url = f"https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName={library_name}"
    resp = requests.get(url, timeout=timeout)
    resp.raise_for_status()
    lines = resp.text.strip().splitlines()
    with open(out_path, "w") as fh:
        for line in lines:
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            term = parts[0]
            desc = parts[1]
            genes = parts[2].split(",") if parts[2] else []
            gh = "\t".join([term, desc] + genes)
            fh.write(gh + "\n")
    return out_path


def summarize_enrichment(df: pd.DataFrame, method: str = "fisher") -> pd.DataFrame:
    """Aggregate enrichment results across communities into a summary table.

    - df: merged enrichment DataFrame produced by `enrich_communities` (must contain `term` and `p_value`/`adjusted_p_value`).
    - method: currently supports 'fisher' (Fisher's combined probability test).

    Returns DataFrame with columns: term, gene_set, count, combined_pvalue.
    """
    if df is None or df.empty:
        return pd.DataFrame()
    # prefer raw p_value if available, otherwise adjusted_p_value
    pcol = None
    for cand in ["p_value", "pvalue", "P-value", "pvalue"]:
        if cand in df.columns:
            pcol = cand
            break
    if pcol is None and "adjusted_p_value" in df.columns:
        pcol = "adjusted_p_value"
    if pcol is None:
        raise ValueError("No p-value column found in enrichment DataFrame")

    records = []
    grouped = df.groupby(["gene_set", "term"]) if "gene_set" in df.columns else df.groupby("term")
    for key, group in grouped:
        pvals = group[pcol].dropna().astype(float).tolist()
        if not pvals:
            continue
        if method == "fisher":
            stat, combined_p = combine_pvalues(pvals, method="fisher")
        else:
            raise ValueError(f"Unknown combine method: {method}")
        if isinstance(key, tuple):
            gene_set, term = key
        else:
            gene_set, term = (None, key)
        records.append({"gene_set": gene_set, "term": term, "count": len(pvals), "combined_pvalue": float(combined_p)})
    out = pd.DataFrame(records)
    if not out.empty:
        out = out.sort_values("combined_pvalue")
    return out
