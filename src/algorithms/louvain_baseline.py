"""Baseline Louvain community detection wrapper."""
from __future__ import annotations

import logging
from typing import Dict, Optional

import networkx as nx

try:
    from community import community_louvain
except ImportError as exc:  # pragma: no cover - import guard
    raise ImportError("Please install python-louvain: pip install python-louvain") from exc

LOGGER = logging.getLogger(__name__)


def run_louvain(graph: nx.Graph, resolution: float = 1.0, random_state: Optional[int] = None) -> Dict[str, int]:
    """Execute the Louvain heuristic and return a node->community map."""
    if graph.number_of_nodes() == 0:
        return {}
    # Print minimal progress so users can see that Louvain started
    print(f"Starting Louvain community detection on {graph.number_of_nodes()} nodes...", flush=True)
    LOGGER.info("Running Louvain on %s nodes", graph.number_of_nodes())
    partition = community_louvain.best_partition(graph, weight="weight", resolution=resolution, random_state=random_state)
    print(f"Louvain finished; detected {len(set(partition.values()))} communities", flush=True)
    LOGGER.info("Detected %s communities", len(set(partition.values())))
    return partition
