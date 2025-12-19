"""Enrichr-based enrichment using gseapy for detected communities.

This module is optional; if gseapy is not installed, functions will raise ImportError.
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, Iterable, List, Optional

import pandas as pd

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
) -> pd.DataFrame:
    """Run Enrichr on each community and return a concatenated DataFrame.

    Columns include community_id, Gene_set, Term, P-value, Adjusted P-value, Odds Ratio, Combined Score, Genes.
    """
    require_gseapy()
    all_rows: List[pd.DataFrame] = []
    for cid, genes in communities.items():
        gene_list = list(dict.fromkeys(genes))
        if not gene_list:
            continue
        try:
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

        all_rows.append(renamed)
    if not all_rows:
        return pd.DataFrame()
    out_df = pd.concat(all_rows, axis=0, ignore_index=True)
    if out_file is not None:
        out_file.parent.mkdir(parents=True, exist_ok=True)
        out_df.to_csv(out_file, index=False)
    return out_df
