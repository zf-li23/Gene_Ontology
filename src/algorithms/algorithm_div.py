"""
Divisive Louvain Algorithm Implementation.
Re-implemented to integrate with the Intelligent Bio Pipeline.
Matches C++ logic with 4 weighting variants:
- div_0_5: weight = 1/sqrt(degree)
- div_1: weight = 1/degree
- div_2: weight = 1/degree^2
- div_avg: weight = 1/(sqrt(d_u)*sqrt(d_v))
"""
import networkx as nx
import community.community_louvain as community_louvain
import math
import logging
from tqdm import tqdm

LOGGER = logging.getLogger(__name__)

def run_algorithm_div(graph: nx.Graph):
    """
    Runs Divisive Louvain variants on the given graph.
    
    Args:
        graph: NetworkX graph object
        
    Returns:
        dict: A dictionary where keys are variant names ('div_0_5', 'div_1', 'div_2', 'div_avg')
              and values are partition dictionaries (node -> community_id).
    """
    results = {}
    
    # Pre-calculate degrees for efficiency
    degrees = dict(graph.degree())
    
    # Define single-node weighting functions for the first 3 variants
    # These calculate the directed weight contribution from a node based on its degree.
    # The final edge weight is the average of the directed weights: (w(u->v) + w(v->u)) / 2
    variants_single = {
        "div_0_5": lambda d: 1.0 / math.sqrt(d) if d > 0 else 0,
        "div_1":   lambda d: 1.0 / d if d > 0 else 0,
        "div_2":   lambda d: 1.0 / (d * d) if d > 0 else 0,
    }

    print(f"Running Divisive Louvain variants on {graph.number_of_nodes()} nodes...", flush=True)
    
    # Process the first 3 variants
    for name, weight_func in tqdm(variants_single.items(), desc="Divisive Variants"):
        H = nx.Graph()
        H.add_nodes_from(graph.nodes())
        
        edges_to_add = []
        for u, v in graph.edges():
            d_u = degrees[u]
            d_v = degrees[v]
            
            # Calculate directed weights
            w_uv = weight_func(d_u)
            w_vu = weight_func(d_v)
            
            # Symmetrize
            weight = (w_uv + w_vu) / 2.0
            
            edges_to_add.append((u, v, weight))
            
        H.add_weighted_edges_from(edges_to_add)
        
        try:
            partition = community_louvain.best_partition(H, weight='weight', random_state=42)
            # Filter small communities (< 3 nodes)
            from collections import Counter
            counts = Counter(partition.values())
            valid_comms = {c for c, count in counts.items() if count >= 3}
            filtered_partition = {n: c for n, c in partition.items() if c in valid_comms}
            results[name] = filtered_partition
        except Exception as e:
            LOGGER.error(f"Failed to run {name}: {e}")
            results[name] = {}

    # Process div_avg variant
    # Weight is 1 / (sqrt(d_u) * sqrt(d_v))
    name = "div_avg"
    H_avg = nx.Graph()
    H_avg.add_nodes_from(graph.nodes())
    edges_avg = []
    for u, v in graph.edges():
        d_u = degrees[u]
        d_v = degrees[v]
        if d_u > 0 and d_v > 0:
            weight = 1.0 / (math.sqrt(d_u) * math.sqrt(d_v))
        else:
            weight = 0
        edges_avg.append((u, v, weight))
    
    H_avg.add_weighted_edges_from(edges_avg)
    
    try:
        partition = community_louvain.best_partition(H_avg, weight='weight', random_state=42)
        # Filter small communities (< 3 nodes)
        from collections import Counter
        counts = Counter(partition.values())
        valid_comms = {c for c, count in counts.items() if count >= 3}
        filtered_partition = {n: c for n, c in partition.items() if c in valid_comms}
        results[name] = filtered_partition
    except Exception as e:
        LOGGER.error(f"Failed to run {name}: {e}")
        results[name] = {}
        
    return results
