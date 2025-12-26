"""Transform raw edge lists into analysis-ready NetworkX graphs."""
from __future__ import annotations

import logging
from typing import Tuple

import networkx as nx
import pandas as pd
from tqdm import tqdm

LOGGER = logging.getLogger(__name__)


def build_graph_from_edgelist(edges: pd.DataFrame, weighted: bool = True) -> nx.Graph:
    """Create an undirected graph from the provided edge table.

    This implementation avoids slow `iterrows()` and potential pandas
    dtype-inference memory spikes by coercing the weight column to numeric
    with `pd.to_numeric` and using NetworkX's optimized
    `from_pandas_edgelist` when possible. For very large tables a simple
    progress log is emitted.
    """
    # Ensure expected columns exist
    if not {"gene_a", "gene_b"}.issubset(edges.columns):
        raise ValueError("Edge table must contain 'gene_a' and 'gene_b' columns")

    df = edges.copy()

    LOGGER.info("Preparing to build graph from %s edges (weighted=%s)", len(df), weighted)

    # Coerce weight column to float if present; fill missing with 1.0
    if weighted and "weight" in df.columns:
        df["weight"] = pd.to_numeric(df["weight"], errors="coerce")
        df["weight"] = df["weight"].fillna(1.0).astype(float)
    else:
        # If not weighted or no weight col, create a default
        df["weight"] = 1.0

    # Use networkx helper which is vectorized and efficient
    try:
        graph = nx.from_pandas_edgelist(df, source="gene_a", target="gene_b", edge_attr=["weight"], create_using=nx.Graph())
    except Exception:
        # Fallback: build incrementally with a progress indicator
        graph = nx.Graph()
        total = len(df)
        for _, row in tqdm(df.iterrows(), total=total, desc="Adding edges", unit="rows"):
            graph.add_edge(row["gene_a"], row["gene_b"], weight=float(row.get("weight", 1.0)))

    LOGGER.info("Graph assembled: %s nodes, %s edges", graph.number_of_nodes(), graph.number_of_edges())
    return graph


def extract_giant_component(graph: nx.Graph) -> nx.Graph:
    """Return the largest connected component; falls back to original graph when trivial."""
    if graph.number_of_nodes() == 0:
        return graph
    components = list(nx.connected_components(graph))
    if not components:
        return graph
    giant = max(components, key=len)
    LOGGER.info("Retaining giant component with %s nodes", len(giant))
    return graph.subgraph(giant).copy()


def normalize_weights(graph: nx.Graph) -> nx.Graph:
    """Scale edge weights into [0, 1] to stabilize downstream algorithms."""
    weights = [data.get("weight", 1.0) for _, _, data in graph.edges(data=True)]
    if not weights:
        return graph
    w_min, w_max = min(weights), max(weights)
    if w_max == w_min:
        LOGGER.info("Weights are constant; skipping normalization")
        return graph
    for u, v, data in graph.edges(data=True):
        data["weight"] = (data.get("weight", 1.0) - w_min) / (w_max - w_min)
    return graph


def prepare_graph(edges: pd.DataFrame, keep_giant_component: bool = True, weighted: bool = True) -> nx.Graph:
    """High-level convenience wrapper handling graph creation and cleanup."""
    graph = build_graph_from_edgelist(edges, weighted=weighted)
    if keep_giant_component:
        graph = extract_giant_component(graph)
    graph = normalize_weights(graph)
    return graph
