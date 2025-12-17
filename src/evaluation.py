"""Evaluation utilities for community quality and biological coherence."""
from __future__ import annotations

import itertools
from statistics import mean
from typing import Dict, List, Sequence


def enrichment_coverage(enrichment_results: Dict[str, List[Dict]]) -> float:
    """Fraction of communities that yielded at least one enriched GO term."""
    if not enrichment_results:
        return 0.0
    enriched = sum(1 for res in enrichment_results.values() if res)
    return enriched / len(enrichment_results)


def jaccard(a: Sequence[str], b: Sequence[str]) -> float:
    set_a, set_b = set(a), set(b)
    union = set_a | set_b
    if not union:
        return 0.0
    return len(set_a & set_b) / len(union)


def community_consistency(communities: Dict[str, List[str]], go_annotations: Dict[str, Sequence[str]]) -> Dict[str, float]:
    """Compute an average pairwise Jaccard similarity of GO annotations per community."""
    scores: Dict[str, float] = {}
    for community_id, genes in communities.items():
        term_sets = [go_annotations.get(gene, []) for gene in genes]
        combos = list(itertools.combinations(term_sets, 2))
        if not combos:
            scores[community_id] = 0.0
            continue
        similarities = [jaccard(a, b) for a, b in combos]
        scores[community_id] = mean(similarities)
    return scores


def summarize_metrics(enrichment_results: Dict[str, List[Dict]], consistency_scores: Dict[str, float]) -> Dict[str, float]:
    """Aggregate a handful of scalar metrics for quick reporting."""
    coverage = enrichment_coverage(enrichment_results)
    avg_consistency = mean(consistency_scores.values()) if consistency_scores else 0.0
    return {
        "enrichment_coverage": coverage,
        "avg_go_consistency": avg_consistency,
    }
