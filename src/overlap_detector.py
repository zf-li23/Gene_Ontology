"""Unsupervised, multi-level, overlapping community detection tailored for large PPIs.

Functions are ordered per user request and include lightweight assertions using
Karate club to guarantee sanity without training labels.
"""
from __future__ import annotations

import json
import logging
import math
from collections import defaultdict, deque
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Sequence, Set, Tuple

import networkx as nx
import numpy as np
from numba import njit
from community import community_louvain

LOGGER = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# 1) isolate_partition
# ---------------------------------------------------------------------------

def _neighbors_array(graph: nx.Graph) -> Tuple[np.ndarray, np.ndarray, Dict[int, int], List[int]]:
    """CSR-like adjacency to accelerate Numba-backed overlap checks."""
    nodes = list(graph.nodes())
    idx = {n: i for i, n in enumerate(nodes)}
    indptr = [0]
    indices: List[int] = []
    for n in nodes:
        nbrs = graph.neighbors(n)
        indices.extend(idx[v] for v in nbrs)
        indptr.append(len(indices))
    return np.array(indptr, dtype=np.int64), np.array(indices, dtype=np.int64), idx, nodes


@njit(cache=True)
def _mark_isolates(indptr: np.ndarray, indices: np.ndarray) -> np.ndarray:
    n = len(indptr) - 1
    used = np.zeros(n, dtype=np.uint8)
    sets = []  # list of lists of node indices
    for u in range(n):
        if used[u]:
            continue
        start = indptr[u]
        end = indptr[u + 1]
        members = [u]
        for k in range(start, end):
            v = indices[k]
            if used[v]:
                continue
            members.append(v)
        for m in members:
            used[m] = 1
        sets.append(members)
    max_len = 0
    for s in sets:
        if len(s) > max_len:
            max_len = len(s)
    out = -np.ones((len(sets), max_len), dtype=np.int64)
    for i, s in enumerate(sets):
        for j, node_idx in enumerate(s):
            out[i, j] = node_idx
    return out


def isolate_partition(graph: nx.Graph) -> List[Set]:
    """Partition G into disjoint isolate-sets (node + neighbors, no overlap).

    Complexity: O(m) time via linear neighbor scan; memory O(n + m).
    """
    indptr, indices, idx_map, nodes = _neighbors_array(graph)
    raw = _mark_isolates(indptr, indices)
    result: List[Set] = []
    for row in raw:
        members = {nodes[i] for i in row if i >= 0}
        if members:
            result.append(members)
    return result


# quick sanity checks (Karate club)
_GK = nx.karate_club_graph()
assert len(isolate_partition(_GK)) >= 1
assert sum(len(s) for s in isolate_partition(_GK)) == _GK.number_of_nodes()


# ---------------------------------------------------------------------------
# 2) parallel_louvain_layer
# ---------------------------------------------------------------------------

def _run_louvain_subgraph(args: Tuple[nx.Graph, int, float, int]) -> Dict:
    subgraph, offset, resolution, seed = args
    partition = community_louvain.best_partition(subgraph, weight="weight", resolution=resolution, random_state=seed)
    return {n: cid + offset for n, cid in partition.items()}


def parallel_louvain_layer(graph: nx.Graph, isolate_sets: List[Set], threads: int = 8,
                           resolution: float = 1.0, random_state: int = 42) -> Dict:
    """Run Louvain independently per isolate-set, concatenating labels.

    Complexity: O(m) across subsets; memory O(n).
    """
    # Fast path for single-thread or tiny inputs to avoid process spawn on import-time asserts
    if threads <= 1 or len(isolate_sets) <= 1:
        merged: Dict = {}
        offset = 0
        for subset in isolate_sets:
            sub = graph.subgraph(subset).copy()
            part = community_louvain.best_partition(
                sub, weight="weight", resolution=resolution, random_state=random_state
            )
            merged.update({n: cid + offset for n, cid in part.items()})
            offset += max(1, sub.number_of_nodes())
        return merged

    from multiprocessing import Pool

    tasks = []
    offset = 0
    for subset in isolate_sets:
        sub = graph.subgraph(subset).copy()
        tasks.append((sub, offset, resolution, random_state))
        offset += max(1, sub.number_of_nodes())
    with Pool(processes=threads) as pool:
        results = pool.map(_run_louvain_subgraph, tasks)
    merged: Dict = {}
    for part in results:
        merged.update(part)
    return merged


assert isinstance(parallel_louvain_layer(_GK, isolate_partition(_GK), threads=1), dict)
assert len(set(parallel_louvain_layer(_GK, isolate_partition(_GK), threads=1).values())) >= 2


# ---------------------------------------------------------------------------
# 3) edge_centrality_seed
# ---------------------------------------------------------------------------

def edge_centrality_seed(graph: nx.Graph) -> List[Tuple]:
    """Rank edges by mean local clustering of endpoints (Wang 2021 inspired).

    Complexity: O(m * <deg>) via clustering computation; memory O(n).
    """
    cc = nx.clustering(graph)
    scores = []
    for u, v in graph.edges():
        scores.append(((cc.get(u, 0.0) + cc.get(v, 0.0)) / 2.0, (u, v)))
    scores.sort(key=lambda x: x[0], reverse=True)
    return [e for _, e in scores]


assert len(edge_centrality_seed(_GK)) == _GK.number_of_edges()
assert edge_centrality_seed(_GK)[0] in _GK.edges()


# ---------------------------------------------------------------------------
# 4) expand_overlapping_communities
# ---------------------------------------------------------------------------

def _attachment_score(graph: nx.Graph, node, community: Set) -> float:
    neighbors = set(graph.neighbors(node))
    if not neighbors:
        return 0.0
    hit = len(neighbors & community)
    return hit / len(neighbors)


def expand_overlapping_communities(graph: nx.Graph, central_edges: List[Tuple], prune: float = 0.3) -> Dict:
    """Grow edge-seeded communities allowing overlaps; prune weak attachments.

    Complexity: O(m) expected; memory O(n + m).
    """
    communities: List[Set] = []
    for (u, v) in central_edges:
        comm: Set = {u, v}
        frontier = deque([u, v])
        while frontier:
            x = frontier.popleft()
            for nbr in graph.neighbors(x):
                if nbr in comm:
                    continue
                if _attachment_score(graph, nbr, comm) >= prune:
                    comm.add(nbr)
                    frontier.append(nbr)
        communities.append(comm)
    memberships: Dict = defaultdict(list)
    for idx, comm in enumerate(communities):
        for node in comm:
            memberships[node].append(idx)
    # prune weak memberships by attachment
    final: Dict = defaultdict(list)
    for idx, comm in enumerate(communities):
        filtered = []
        for node in comm:
            score = _attachment_score(graph, node, comm)
            if score >= prune:
                filtered.append(node)
        if filtered:
            for node in filtered:
                final[node].append(idx)
    return dict(final)


_comm_labels = expand_overlapping_communities(_GK, edge_centrality_seed(_GK)[:5], prune=0.2)
assert isinstance(_comm_labels, dict)
assert any(len(v) > 1 for v in _comm_labels.values())  # overlapping expected on karate


# ---------------------------------------------------------------------------
# 5) hierarchy_merge
# ---------------------------------------------------------------------------

@dataclass
class TreeNode:
    module_id: str
    depth: int
    members: List
    children: List["TreeNode"] = field(default_factory=list)

    def to_dict(self) -> Dict:
        return {
            "id": self.module_id,
            "depth": self.depth,
            "members": self.members,
            "children": [c.to_dict() for c in self.children],
        }


def hierarchy_merge(graph: nx.Graph, bottom_labels: Dict) -> TreeNode:
    """Build a multi-branch hierarchy by agglomerating overlapping modules.

    Complexity: O(k^2 * s) where k=#modules, s avg size; memory O(k * s).
    """
    # invert bottom_labels to module -> members
    module_to_nodes: Dict[int, Set] = defaultdict(set)
    for node, mods in bottom_labels.items():
        for m in mods:
            module_to_nodes[m].add(node)
    level_modules = list(module_to_nodes.values())
    level_nodes: List[TreeNode] = [TreeNode(f"L0M{i}", 0, sorted(mod)) for i, mod in enumerate(level_modules)]
    depth = 0
    while len(level_nodes) > 1:
        depth += 1
        # build overlap graph among current modules
        ov_graph = nx.Graph()
        ov_graph.add_nodes_from(range(len(level_nodes)))
        for i in range(len(level_nodes)):
            for j in range(i + 1, len(level_nodes)):
                a, b = set(level_nodes[i].members), set(level_nodes[j].members)
                if not a or not b:
                    continue
                jacc = len(a & b) / len(a | b)
                if jacc >= 0.25:
                    ov_graph.add_edge(i, j)
        components = list(nx.connected_components(ov_graph)) if ov_graph.number_of_edges() > 0 else [{i} for i in range(len(level_nodes))]
        next_level: List[TreeNode] = []
        for cid, comp in enumerate(components):
            merged_members: Set = set()
            children: List[TreeNode] = []
            for idx in comp:
                merged_members.update(level_nodes[idx].members)
                children.append(level_nodes[idx])
            parent = TreeNode(f"L{depth}M{cid}", depth, sorted(merged_members), children)
            next_level.append(parent)
        level_nodes = next_level
    return level_nodes[0]


_root = hierarchy_merge(_GK, _comm_labels)
assert _root.depth >= 1
assert len(_root.children) >= 1


# ---------------------------------------------------------------------------
# 6) topological_validate
# ---------------------------------------------------------------------------

def topological_validate(graph: nx.Graph, tree: TreeNode) -> Dict[str, float]:
    """Return Ravasz-style topology fingerprints: C(k) slope, clustering, modularity.

    Complexity: O(m) for clustering; memory O(n).
    """
    degrees = np.array([d for _, d in graph.degree()], dtype=float)
    clustering = np.array(list(nx.clustering(graph).values()), dtype=float)
    valid = (degrees > 0) & (clustering > 0)
    if valid.sum() >= 2:
        slope, _ = np.polyfit(np.log(degrees[valid]), np.log(clustering[valid]), 1)
    else:
        slope = float("nan")
    avg_c = float(clustering.mean()) if len(clustering) else 0.0
    # use first-level children as a hard partition for modularity estimation
    if tree.children:
        part_sets = [set(child.members) for child in tree.children]
    else:
        part_sets = [set(tree.members)]
    # Ensure a valid partition (no overlaps) for modularity by greedily assigning each node once
    disjoint_sets: List[Set] = []
    assigned: Set = set()
    for comm in part_sets:
        pruned = {n for n in comm if n not in assigned}
        if pruned:
            disjoint_sets.append(pruned)
            assigned.update(pruned)
    if not disjoint_sets:
        disjoint_sets = [set(graph.nodes())]
    modularity = nx.algorithms.community.quality.modularity(graph, disjoint_sets)
    return {
        "C(k)_slope": float(slope),
        "avg_clustering": avg_c,
        "modularity": float(modularity),
    }


_val = topological_validate(_GK, _root)
assert "modularity" in _val and isinstance(_val["modularity"], float)
assert math.isnan(_val["C(k)_slope"]) or _val["C(k)_slope"] < 0.5


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def run_pipeline(edge_file: str, threads: int = 8, resolution: float = 1.0, prune: float = 0.3) -> Dict:
    """Convenience runner to integrate all steps and emit JSON-ready dict."""
    df = None
    try:
        df = np.genfromtxt(edge_file, dtype=str)
    except Exception as exc:  # pragma: no cover
        raise RuntimeError(f"Failed to read edge file {edge_file}: {exc}")
    if df.ndim == 1:
        raise ValueError("Edge file must have at least two columns")
    u = df[:, 0].tolist()
    v = df[:, 1].tolist()
    w = df[:, 2].astype(float).tolist() if df.shape[1] > 2 else [1.0] * len(u)
    G = nx.Graph()
    for a, b, wt in zip(u, v, w):
        if a == b:
            continue
        if G.has_edge(a, b):
            if wt > G[a][b].get("weight", 1.0):
                G[a][b]["weight"] = wt
        else:
            G.add_edge(a, b, weight=wt)
    iso_sets = isolate_partition(G)
    labels = parallel_louvain_layer(G, iso_sets, threads=threads, resolution=resolution)
    seeds = edge_centrality_seed(G)
    bottom = expand_overlapping_communities(G, seeds, prune=prune)
    tree = hierarchy_merge(G, bottom)
    topo = topological_validate(G, tree)
    modules = []
    for node, mods in bottom.items():
        for m in mods:
            modules.append({"id": f"L0M{m}", "depth": 0, "members": node, "overlap": True})
    return {
        "levels": tree.depth + 1,
        "modules": modules,
        "tree": tree.to_dict(),
        "topology": topo,
    }


if __name__ == "__main__":
    import argparse
    logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")
    parser = argparse.ArgumentParser(description="Overlapping hierarchical community detection for PPIs")
    parser.add_argument("edge_file", help="TSV with at least two columns: gene_a gene_b [weight]")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--prune", type=float, default=0.3)
    parser.add_argument("--resolution", type=float, default=1.0)
    args = parser.parse_args()
    output = run_pipeline(args.edge_file, threads=args.threads, resolution=args.resolution, prune=args.prune)
    print(json.dumps(output, indent=2))
