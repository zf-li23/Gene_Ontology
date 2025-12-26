"""Divisive Louvain algorithm variants."""
from __future__ import annotations

import math
import logging
from typing import Dict, List, Tuple
import networkx as nx
from tqdm import tqdm

try:
    from community import community_louvain
except ImportError:
    community_louvain = None

LOGGER = logging.getLogger(__name__)

def run_algorithm_div(graph: nx.Graph) -> Dict[str, Dict[str, int]]:
    """
    Run 4 variants of divisive Louvain (0.5, 1, 2, Average).
    Returns a dictionary of {variant_name: partition}.
    """
    if community_louvain is None:
        raise ImportError("python-louvain not installed")

    # Pre-calculate degrees for all nodes
    # Map node -> degree
    degrees = dict(graph.degree(weight=None)) # Use unweighted degree for m?
    # The original code used 'm' which was incremented per edge.
    # sp[xi].m += 1. So it is unweighted degree.
    
    variants = {
        "div_0_5": lambda d_u, d_v: 1.0 / math.sqrt(d_u) if d_u > 0 else 0,
        "div_1":   lambda d_u, d_v: 1.0 / d_u if d_u > 0 else 0,
        "div_2":   lambda d_u, d_v: 1.0 / (d_u * d_u) if d_u > 0 else 0,
        "div_avg": lambda d_u, d_v: 1.0 / (math.sqrt(d_u) * math.sqrt(d_v)) if d_u > 0 and d_v > 0 else 0
    }

    results = {}
    
    print(f"Running Divisive Louvain variants on {graph.number_of_nodes()} nodes...", flush=True)
    
    for name, weight_func in tqdm(variants.items(), desc="Divisive Variants"):
        # Build weighted graph
        # Note: Original code produced asymmetric weights for 0.5, 1, 2.
        # Louvain requires undirected. We will symmetrize by averaging (w_uv + w_vu) / 2.
        # For div_avg, it is already symmetric.
        # For div_1 (1/d_u), w_uv = 1/d_u, w_vu = 1/d_v.
        # Symmetrized = 0.5 * (1/d_u + 1/d_v).
        
        g_var = nx.Graph()
        g_var.add_nodes_from(graph.nodes())
        
        edges_to_add = []
        for u, v in graph.edges():
            d_u = degrees[u]
            d_v = degrees[v]
            
            w_uv = weight_func(d_u, d_v)
            w_vu = weight_func(d_v, d_u)
            
            final_weight = (w_uv + w_vu) / 2.0
            edges_to_add.append((u, v, final_weight))
            
        g_var.add_weighted_edges_from(edges_to_add)
        
        # Run Louvain
        try:
            partition = community_louvain.best_partition(g_var, weight='weight', random_state=42)
            # Filter small communities (size < 3)
            # partition is node -> comm_id
            # We need to check sizes
            from collections import Counter
            counts = Counter(partition.values())
            valid_comms = {c for c, count in counts.items() if count >= 3}
            
            # Filter partition
            filtered_partition = {n: c for n, c in partition.items() if c in valid_comms}
            
            results[name] = filtered_partition
        except Exception as e:
            LOGGER.warning(f"Failed to run Louvain for {name}: {e}")
            results[name] = {}

    return results
