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
import numpy as np
import contextlib
import os

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
    disable_progress: bool = False,
    silence_gseapy: bool = False,
) -> pd.DataFrame:
    """Run Enrichr on each community and return a concatenated DataFrame.

    Columns include community_id, Gene_set, Term, P-value, Adjusted P-value, Odds Ratio, Combined Score, Genes.
    """
    require_gseapy()
    all_rows: List[pd.DataFrame] = []
    # try to use tqdm for progress if available
    # Always try to show a tqdm progress bar when available (even for local GMTs)
    try:
        from tqdm import tqdm

        iterator = tqdm(communities.items(), desc="Enrichr communities")
    except Exception:
        iterator = communities.items()

    # Optionally silence noisy gseapy logger output (useful for local GMTs)
    gseapy_logger = logging.getLogger("gseapy")
    old_level = None
    if silence_gseapy:
        try:
            old_level = gseapy_logger.level
            gseapy_logger.setLevel(logging.CRITICAL)
        except Exception:
            old_level = None

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
                # gp.enrichr may emit logs/prints even when logging level is adjusted; when requested,
                # redirect stderr to suppress noisy messages coming from gseapy internals.
                if silence_gseapy:
                    with contextlib.redirect_stderr(open(os.devnull, "w")):
                        enr = gp.enrichr(
                            gene_list=gene_list,
                            gene_sets=gene_sets,
                            organism=organism,
                            outdir=None,
                            cutoff=cutoff,
                        )
                else:
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
        # `enr.results` can be a DataFrame, a list of dicts, or other types depending on gseapy usage
        res = enr.results
        if isinstance(res, pd.DataFrame):
            df = res.copy()
        else:
            try:
                df = pd.DataFrame(res)
            except Exception:
                LOGGER.warning("Unrecognized enrichr results format for %s; skipping", cid)
                continue

        # skip empty results
        if df is None or (hasattr(df, "empty") and df.empty):
            LOGGER.info("No enrichment results for %s; skipping", cid)
            continue

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
        if term_col in df.columns:
            renamed["term"] = df[term_col]
        else:
            # fallback: use first column if available, otherwise skip this community
            if df.shape[1] > 0:
                renamed["term"] = df.iloc[:, 0]
            else:
                LOGGER.warning("No term column available for %s; skipping", cid)
                continue
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
        # restore gseapy logger level
        if silence_gseapy and old_level is not None:
            try:
                gseapy_logger.setLevel(old_level)
            except Exception:
                pass
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
    # restore gseapy logger level
    if silence_gseapy and old_level is not None:
        try:
            gseapy_logger.setLevel(old_level)
        except Exception:
            pass
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


def compute_concordance_scores(df: pd.DataFrame, communities: Dict[str, List[str]], cutoff: float = 0.05) -> pd.DataFrame:
    """Compute concordance metrics between detected communities and gene-set libraries.

    Returns a DataFrame with per-gene_set metrics plus an overall row where gene_set==None.

    Metrics:
    - n_terms: number of (term,community) observations
    - n_significant_communities: number of communities with at least one term with adjusted p < cutoff
    - frac_significant_communities: fraction of communities with at least one significant term
    - mean_score: mean per-community score s_c where s_c = max(0, -log10(min_q))
    - normalized_score: mean_score divided by -log10(cutoff) (range ~0..1)
    - weighted_score: community-size-weighted mean_score
    """
    if df is None or df.empty:
        return pd.DataFrame()

    # identify p-value column (prefer adjusted_p_value)
    pcol = None
    for cand in ["adjusted_p_value", "Adjusted P-value", "p_value", "P-value"]:
        if cand in df.columns:
            pcol = cand
            break
    if pcol is None:
        raise ValueError("No p-value column found for concordance scoring")

    # ensure community list includes all communities
    all_communities = list(communities.keys())
    comm_sizes = {k: len(v) for k, v in communities.items()}

    records = []
    gene_sets = df["gene_set"].dropna().unique().tolist() if "gene_set" in df.columns else [None]
    # compute for each gene_set separately, and overall combined
    def _compute(subdf, label):
        # per-community min p
        grouped = subdf.groupby("community_id")[pcol].min()
        scores = {}
        for cid in all_communities:
            min_p = None
            if cid in grouped.index:
                try:
                    min_p = float(grouped.loc[cid])
                except Exception:
                    min_p = None
            if min_p is None or np.isnan(min_p):
                s = 0.0
            else:
                s = max(0.0, -np.log10(min_p))
            scores[cid] = s
        mean_score = float(np.mean(list(scores.values()))) if scores else 0.0
        denom = -np.log10(cutoff) if cutoff and cutoff > 0 else 1.0
        normalized = float(mean_score / denom) if denom else mean_score
        # fraction of communities with at least one significant hit
        n_sig = sum(1 for v in scores.values() if v > 0 and (10 ** (-v)) < cutoff)
        frac_sig = n_sig / max(1, len(all_communities))
        # weighted by community size
        weights = np.array([comm_sizes.get(cid, 1) for cid in all_communities], dtype=float)
        vals = np.array([scores.get(cid, 0.0) for cid in all_communities], dtype=float)
        weighted_score = float(np.average(vals, weights=weights)) if weights.sum() > 0 else mean_score
        return {
            "gene_set": label,
            "n_terms": int(len(subdf)),
            "n_significant_communities": int(n_sig),
            "frac_significant_communities": float(frac_sig),
            "mean_score": float(mean_score),
            "normalized_score": float(normalized),
            "weighted_score": float(weighted_score),
        }

    for gs in gene_sets:
        sub = df[df["gene_set"] == gs] if gs is not None else df
        records.append(_compute(sub, gs))

    # overall
    records.append(_compute(df, None))
    out = pd.DataFrame(records)
    return out
