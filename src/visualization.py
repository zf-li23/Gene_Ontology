"""Plotting helpers for communities, hierarchies, and enrichment summaries."""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, List

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import networkx as nx
import numpy as np
import warnings

from .algorithms.hierarchical_louvain import HierarchyNode

LOGGER = logging.getLogger(__name__)


COLORS = [
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b",
    "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
]


def _scale(values, vmin=None, vmax=None):
    arr = np.array(list(values), dtype=float)
    if arr.size == 0:
        return arr
    if vmin is None:
        vmin = np.percentile(arr, 5)
    if vmax is None:
        vmax = np.percentile(arr, 95)
    if vmax == vmin:
        vmax = vmin + 1e-9
    arr = np.clip(arr, vmin, vmax)
    return (arr - vmin) / (vmax - vmin)


def plot_network(graph: nx.Graph, partition: Dict[str, int], output_path: Path, layout: str = "spring") -> None:
    """Draw the interaction network with improved aesthetics (scaled sizes/colors)."""
    if graph.number_of_nodes() == 0:
        LOGGER.warning("Graph empty; skipping network plot")
        return
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if layout == "spring":
        pos = nx.spring_layout(graph, seed=42, weight="weight")
    else:
        pos = nx.kamada_kawai_layout(graph, weight="weight")
    communities = partition or {node: 0 for node in graph.nodes()}
    comm_ids = np.array([communities.get(node, 0) for node in graph.nodes()])
    cmap = cm.get_cmap("tab20")
    node_colors = [cmap(c % cmap.N) for c in comm_ids]

    degrees = dict(graph.degree(weight="weight"))
    node_sizes = 200 * (_scale(degrees.values()) + 0.2)

    weights = [data.get("weight", 1.0) for _, _, data in graph.edges(data=True)]
    edge_widths = 1.0 + 2.5 * _scale(weights)
    edge_alphas = 0.2 + 0.5 * _scale(weights)

    plt.figure(figsize=(9, 7))
    nx.draw_networkx_edges(
        graph,
        pos=pos,
        width=edge_widths,
        alpha=edge_alphas,
        edge_color="#7f7f7f",
    )
    nx.draw_networkx_nodes(
        graph,
        pos=pos,
        node_color=node_colors,
        node_size=node_sizes,
        linewidths=0.3,
        edgecolors="#222222",
        alpha=0.9,
    )
    plt.title("Interaction network with community colors", fontsize=12)
    plt.axis("off")
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()


def plot_hierarchy(root: HierarchyNode, output_path: Path) -> None:
    """Visualize the hierarchy as a tree plot with better spacing."""
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
        LOGGER.debug("graphviz_layout unavailable; falling back to spring layout")
        pos = nx.spring_layout(tree)
    depths = nx.get_node_attributes(tree, "depth")
    depth_vals = [depths.get(n, 0) for n in tree.nodes()]
    node_sizes = 800 + 200 * np.array(depth_vals)
    plt.figure(figsize=(9, 7))
    nx.draw(
        tree,
        pos,
        with_labels=True,
        arrows=False,
        node_size=node_sizes,
        font_size=8,
        node_color=[COLORS[d % len(COLORS)] for d in depth_vals],
        edge_color="#888888",
    )
    plt.title("Hierarchical community tree")
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message=".*tight_layout.*", category=UserWarning)
        try:
            plt.tight_layout()
        except Exception:
            pass
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
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message=".*tight_layout.*", category=UserWarning)
        try:
            plt.tight_layout()
        except Exception:
            pass
    plt.savefig(output_path, dpi=300)
    plt.close()
