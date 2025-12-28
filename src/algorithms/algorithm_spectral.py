"""Spectral Clustering algorithm."""
from __future__ import annotations

import logging
import warnings
from typing import Dict, List
from collections import Counter
import networkx as nx
import numpy as np
from sklearn.cluster import SpectralClustering
import scipy.sparse.linalg

LOGGER = logging.getLogger(__name__)

def run_algorithm_spectral(graph: nx.Graph, n_clusters: int = None) -> Dict[str, List[str]]:
    """
    Run Spectral Clustering.
    If n_clusters is None, estimates it using the eigengap heuristic.
    Returns a dictionary of community_id -> list of members.
    """
    nodes = list(graph.nodes())
    n_nodes = len(nodes)
    if n_nodes == 0:
        return {}
    
    # Create adjacency matrix
    adj_matrix = nx.to_numpy_array(graph, nodelist=nodes)
    
    if n_clusters is None:
        # Eigengap heuristic to determine number of clusters
        try:
            # Use normalized Laplacian
            L = nx.normalized_laplacian_matrix(graph, nodelist=nodes)
            
            # We need to find the first k eigenvalues
            # For small graphs, we can compute more, for large graphs, limit it.
            # Heuristic: Check up to min(N, 20) eigenvalues
            k_check = min(n_nodes - 1, 20)
            
            if k_check > 1:
                # 'SM' = Smallest Magnitude. 
                # Note: eigsh can be unstable for very small matrices or specific conditions
                eigenvalues = scipy.sparse.linalg.eigsh(L, k=k_check, which='SM', return_eigenvectors=False)
                eigenvalues = np.sort(eigenvalues)
                
                # Calculate gaps: lambda_{i+1} - lambda_i
                gaps = np.diff(eigenvalues)
                
                # Find the index of the largest gap
                # We usually prefer a smaller k, so if there are multiple similar gaps, pick the earlier one?
                # argmax returns the first occurrence of the max value
                best_k = np.argmax(gaps) + 1
                
                # Ensure at least 2 clusters if possible, and not too many
                n_clusters = max(2, int(best_k))
                print(f"Spectral Clustering: Estimated k={n_clusters} (max gap at index {best_k})", flush=True)
            else:
                n_clusters = 2
        except Exception as e:
            LOGGER.warning(f"Eigengap estimation failed: {e}. Defaulting to 8.")
            n_clusters = min(n_nodes, 8)

    print(f"Running Spectral Clustering with k={n_clusters} on {n_nodes} nodes...", flush=True)
    
    try:
        # affinity='precomputed' means we pass the adjacency matrix (or affinity matrix)
        # SpectralClustering expects an affinity matrix where entries are similarities.
        # Adjacency matrix works as affinity.
        # Suppress warnings about disconnected graph
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=UserWarning, module="sklearn")
            sc = SpectralClustering(n_clusters=n_clusters, affinity='precomputed', random_state=42, n_init=10)
            labels = sc.fit_predict(adj_matrix)
        
        # Convert to partition dict
        partition = {}
        for i, label in enumerate(labels):
            partition[nodes[i]] = label
            
        # Filter small communities (< 3)
        counts = Counter(partition.values())
        valid_comms = {c for c, count in counts.items() if count >= 3}
        
        # Format result: comm_id -> list of members
        final_communities = {}
        for node, comm_id in partition.items():
            if comm_id in valid_comms:
                cid_str = f"Spectral_{comm_id}"
                if cid_str not in final_communities:
                    final_communities[cid_str] = []
                final_communities[cid_str].append(node)
                
        return final_communities
        
    except Exception as e:
        LOGGER.error(f"Spectral Clustering failed during fit: {e}")
        return {}
