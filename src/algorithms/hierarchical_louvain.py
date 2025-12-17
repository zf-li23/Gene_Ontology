"""Recursive Louvain-based community detection producing a hierarchy."""
from __future__ import annotations

import itertools
import logging
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, List, Optional

import networkx as nx

try:
    from community import community_louvain
except ImportError as exc:  # pragma: no cover
    raise ImportError("Please install python-louvain: pip install python-louvain") from exc

LOGGER = logging.getLogger(__name__)


@dataclass
class HierarchyNode:
    """Lightweight tree node describing a community and its children."""
    name: str
    members: List[str]
    depth: int
    children: List["HierarchyNode"] = field(default_factory=list)

    def to_dict(self) -> Dict:
        return {
            "name": self.name,
            "depth": self.depth,
            "members": self.members,
            "children": [child.to_dict() for child in self.children],
        }


class HierarchicalLouvain:
    """Encapsulates recursion logic and hyper-parameters."""

    def __init__(self, min_size: int = 5, max_depth: int = 3, modularity_threshold: float = 1e-4,
                 resolution: float = 1.0, random_state: Optional[int] = None) -> None:
        self.min_size = min_size
        self.max_depth = max_depth
        self.modularity_threshold = modularity_threshold
        self.resolution = resolution
        self.random_state = random_state
        self._counter = itertools.count()

    def _should_split(self, graph: nx.Graph, partition: Dict[str, int]) -> bool:
        if graph.number_of_nodes() < self.min_size:
            return False
        if len(set(partition.values())) <= 1:
            return False
        modularity = community_louvain.modularity(partition, graph, weight="weight")
        LOGGER.debug("Depth split modularity=%.4f", modularity)
        return modularity >= self.modularity_threshold

    def _build_node(self, graph: nx.Graph, depth: int) -> HierarchyNode:
        node_name = f"C{depth}_{next(self._counter)}"
        node = HierarchyNode(name=node_name, members=list(graph.nodes()), depth=depth)
        if depth >= self.max_depth:
            return node
        if graph.number_of_nodes() < self.min_size:
            return node
        partition = community_louvain.best_partition(
            graph,
            weight="weight",
            resolution=self.resolution,
            random_state=self.random_state,
        )
        if not self._should_split(graph, partition):
            return node
        grouped: Dict[int, List[str]] = defaultdict(list)
        for vertex, community_id in partition.items():
            grouped[community_id].append(vertex)
        for comm_id, members in grouped.items():
            if len(members) < self.min_size:
                continue
            subgraph = graph.subgraph(members).copy()
            LOGGER.debug("Recursing into community %s at depth %s with %s genes", comm_id, depth + 1, len(members))
            child = self._build_node(subgraph, depth + 1)
            node.children.append(child)
        return node

    def fit(self, graph: nx.Graph) -> HierarchyNode:
        LOGGER.info("Building hierarchical communities (max_depth=%s)", self.max_depth)
        if graph.number_of_nodes() == 0:
            return HierarchyNode(name="empty", members=[], depth=0)
        return self._build_node(graph, depth=0)


def flatten_hierarchy(root: HierarchyNode) -> Dict[int, List[List[str]]]:
    """Convert the hierarchy into level -> list of communities for downstream analysis."""
    levels: Dict[int, List[List[str]]] = defaultdict(list)

    def _walk(node: HierarchyNode) -> None:
        levels[node.depth].append(node.members)
        for child in node.children:
            _walk(child)

    _walk(root)
    return dict(levels)
