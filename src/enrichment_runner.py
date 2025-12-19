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
        df.insert(0, "community_id", cid)
        if top_terms:
            df = df.sort_values("Adjusted P-value").groupby("Gene_set").head(top_terms)
        all_rows.append(df)
    if not all_rows:
        return pd.DataFrame()
    out_df = pd.concat(all_rows, axis=0, ignore_index=True)
    if out_file is not None:
        out_file.parent.mkdir(parents=True, exist_ok=True)
        out_df.to_csv(out_file, index=False)
    return out_df
