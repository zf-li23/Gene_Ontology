"""Plotting helpers for communities, hierarchies, and enrichment summaries."""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, List

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

from .algorithms.hierarchical_louvain import HierarchyNode

LOGGER = logging.getLogger(__name__)


COLORS = [
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b",
    "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
]


def plot_network(graph: nx.Graph, partition: Dict[str, int], output_path: Path, layout: str = "spring") -> None:
    """Draw the interaction network with community colors."""
    if graph.number_of_nodes() == 0:
        LOGGER.warning("Graph empty; skipping network plot")
        return
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if layout == "spring":
        pos = nx.spring_layout(graph, seed=42, weight="weight")
    else:
        pos = nx.kamada_kawai_layout(graph, weight="weight")
    communities = partition or {node: 0 for node in graph.nodes()}
    colors = [COLORS[communities.get(node, 0) % len(COLORS)] for node in graph.nodes()]
    plt.figure(figsize=(8, 6))
    nx.draw_networkx(
        graph,
        pos=pos,
        node_color=colors,
        with_labels=False,
        node_size=120,
        width=0.8,
        edge_color="#C0C0C0",
    )
    plt.title("Gene interaction network with Louvain communities")
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()


def plot_hierarchy(root: HierarchyNode, output_path: Path) -> None:
    """Visualize the hierarchy as a tree plot."""
    if not root.members:
        LOGGER.warning("Hierarchy empty; skipping plot")
        return
    output_path.parent.mkdir(parents=True, exist_ok=True)
    tree = nx.DiGraph()

    def _add(node: HierarchyNode) -> None:
        label = f"{node.name}\n(n={len(node.members)})"
        tree.add_node(label, depth=node.depth)
        for child in node.children:
            child_label = f"{child.name}\n(n={len(child.members)})"
            tree.add_edge(label, child_label)
            _add(child)

    _add(root)
    try:
        pos = nx.nx_pydot.graphviz_layout(tree, prog="dot") if tree.number_of_nodes() < 200 else nx.spring_layout(tree)
    except (ImportError, nx.NetworkXException):  # pragma: no cover - graceful fallback
        LOGGER.warning("graphviz_layout unavailable; falling back to spring layout")
        pos = nx.spring_layout(tree)
    plt.figure(figsize=(8, 6))
    nx.draw(tree, pos, with_labels=True, arrows=False, node_size=3000, font_size=8)
    plt.title("Hierarchical community tree")
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()


def plot_enrichment_bar(enrichment_results: Dict[str, List[Dict]], output_path: Path, top_n: int = 5) -> None:
    """Plot the top-N enriched terms across all communities."""
    flattened = []
    for community_id, terms in enrichment_results.items():
        for term in terms:
            flattened.append({
                "community": community_id,
                "term": term["term"],
                "q_value": term["q_value"],
            })
    if not flattened:
        LOGGER.warning("No enriched terms available; skipping bar plot")
        return
    flattened.sort(key=lambda x: x["q_value"])
    top = flattened[:top_n]
    labels = [f"{item['community']}\n{item['term']}" for item in top]
    scores = [-np.log10(item["q_value"]) for item in top]
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.figure(figsize=(8, 4))
    plt.barh(labels, scores, color="#2ca02c")
    plt.xlabel("-log10(q-value)")
    plt.title("Top enriched GO terms")
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()
