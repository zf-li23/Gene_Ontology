"""GO enrichment utilities based on the hypergeometric test."""
from __future__ import annotations

import logging
from collections import defaultdict
from typing import Dict, Iterable, List, Sequence

import numpy as np
from scipy.stats import hypergeom
from statsmodels.stats.multitest import fdrcorrection

LOGGER = logging.getLogger(__name__)


def invert_annotations(go_annotations: Dict[str, Sequence[str]]) -> Dict[str, List[str]]:
    """Create a GO term -> genes mapping for quick overlap counting."""
    inverted: Dict[str, List[str]] = defaultdict(list)
    for gene, terms in go_annotations.items():
        for term in terms:
            inverted[term].append(gene)
    return inverted


def hypergeometric_enrichment(genes: Iterable[str], background: Sequence[str], go_annotations: Dict[str, Sequence[str]],
                              fdr_threshold: float = 0.05, min_genes: int = 2) -> List[Dict]:
    """Return enriched GO terms for ``genes`` using a hypergeometric test + FDR correction."""
    gene_list = list(dict.fromkeys(genes))  # preserve order, drop duplicates
    background_set = list(dict.fromkeys(background))
    N = len(background_set)
    if N == 0 or not gene_list:
        return []
    go_to_genes = invert_annotations(go_annotations)
    p_values = []
    records = []
    community_set = set(gene_list)
    for term, annotated_genes in go_to_genes.items():
        overlap = community_set.intersection(annotated_genes)
        if len(overlap) < min_genes:
            continue
        K = len(set(annotated_genes).intersection(background_set))
        n = len(community_set)
        k = len(overlap)
        if K == 0:
            continue
        # Survival function gives P(X >= k)
        p_val = hypergeom.sf(k - 1, N, K, n)
        p_values.append(p_val)
        records.append({
            "term": term,
            "overlap": k,
            "community_size": n,
            "background_hits": K,
            "p_value": p_val,
            "genes": sorted(overlap),
            "enrichment": (k / n) / (K / N),
        })
    if not p_values:
        return []
    rejected, q_values = fdrcorrection(p_values, alpha=fdr_threshold)
    enriched = []
    for rec, reject, q_val in zip(records, rejected, q_values):
        if reject:
            rec["q_value"] = float(q_val)
            enriched.append(rec)
    enriched.sort(key=lambda x: x["q_value"])
    return enriched


def background_from_strategy(strategy: str, graph_genes: Iterable[str], go_annotations: Dict[str, Sequence[str]]) -> List[str]:
    """Provide a flexible way to define the background gene universe."""
    strategy = strategy.lower()
    if strategy == "union":
        go_union = set(go_annotations.keys())
        go_union.update(graph_genes)
        return sorted(go_union)
    if strategy == "graph":
        return sorted(set(graph_genes))
    raise ValueError(f"Unknown background strategy: {strategy}")


def enrich_all_communities(communities: Dict[str, List[str]], background: Sequence[str], go_annotations: Dict[str, Sequence[str]],
                           fdr_threshold: float, min_genes: int) -> Dict[str, List[Dict]]:
    """Run enrichment for every community keyed by identifier."""
    results: Dict[str, List[Dict]] = {}
    for community_id, genes in communities.items():
        enriched = hypergeometric_enrichment(genes, background, go_annotations, fdr_threshold=fdr_threshold, min_genes=min_genes)
        results[community_id] = enriched
    return results
