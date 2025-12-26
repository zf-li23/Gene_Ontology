"""Algorithm 2: DFS-based community detection (Loose/Tight variant)."""
from __future__ import annotations

import logging
from typing import Dict, List, Optional, Any
import networkx as nx
from tqdm import tqdm

LOGGER = logging.getLogger(__name__)

class Node:
    def __init__(self, name: str):
        self.name = name
        self.neighbors: List[Node] = []
        self.m = 0  # degree (or similar metric)
        self.ex = 0  # status: 0 = active, -1 = removed/processed

    def __repr__(self):
        return f"Node({self.name}, m={self.m})"

class Algorithm2:
    def __init__(self, graph: nx.Graph, k1: float = 1, k2: float = 2.3):
        self.graph = graph
        self.k1 = k1
        self.k2 = k2
        self.nodes: List[Node] = []
        self.node_map: Dict[str, Node] = {}
        self.c: Dict[Node, int] = {}
        self.add: List[List[Node]] = []
        self.ack: List[int] = []
        self.ec: Dict[Node, int] = {}
        
        self._build_internal_graph()

    def _build_internal_graph(self):
        # Convert nx.Graph to internal Node structure
        for n_name in self.graph.nodes():
            node = Node(str(n_name))
            self.nodes.append(node)
            self.node_map[str(n_name)] = node
        
        for u, v in self.graph.edges():
            nu = self.node_map[str(u)]
            nv = self.node_map[str(v)]
            nu.neighbors.append(nv)
            nv.neighbors.append(nu)
            nu.m += 1
            nv.m += 1

    def dfs(self, u: Node, dep: int, mp: Dict[Node, int], group: List[Node]):
        self.c[u] = dep
        res = []
        temp = []
        
        # Collect valid neighbors
        for v in u.neighbors:
            if v.ex == -1:
                continue
            if self.c.get(v, 0) == 0 or self.c.get(v, 0) > dep + 1:
                mp[v] = mp.get(v, 0) + 1
        
        # Filter neighbors based on depth condition
        for k, v in mp.items():
            if self.c.get(k, 0) == 0 or self.c.get(k, 0) > dep + 1:
                res.append((k, v))
        
        # Sort by count desc, then name (or object id)
        res.sort(key=lambda x: (-x[1], x[0].name))
        
        ck = 0
        for r in res:
            node_r, count_r = r
            threshold = dep - max(0, int((dep - self.k1) / self.k2))
            if count_r >= threshold:
                ck += 1
                temp.append(node_r)
                del mp[node_r]
                group.append(node_r)
                # Pass a copy of mp to avoid pollution from children if that's the issue,
                # or ensure we revert changes. The original code passed 'mp' directly,
                # which in Python modifies the same dict. This causes counts to accumulate
                # across branches, likely causing infinite loops/growth.
                # We will pass a copy to enforce proper scoping.
                self.dfs(node_r, dep + 1, mp.copy(), group)
                group.pop()
                mp[node_r] = count_r

            else:
                break
        
        for t in temp:
            self.c[t] = 0
            
        if ck == 0:
            # Found a group
            # Create a copy of the current group
            new_group = sorted(group[:], key=lambda x: x.name)
            self.add.append(new_group)

    def merge(self, a_idx: int, b_idx: int):
        # Merge group b into group a
        group_a = self.add[a_idx]
        group_b = self.add[b_idx]
        
        combined = group_a + group_b
        # Unique and sort
        combined = sorted(list(set(combined)), key=lambda x: x.name)
        self.add[a_idx] = combined
        self.add[b_idx] = [] # Clear b

    def compare(self, a_idx: int, b_idx: int):
        group_a = self.add[a_idx]
        group_b = self.add[b_idx]
        
        if not group_a or not group_b:
            return

        # Intersection
        set_a = set(group_a)
        set_b = set(group_b)
        intersection = set_a.intersection(set_b)
        simi = len(intersection)
        
        target_len = min(len(group_a), len(group_b))
        if target_len == 0:
            return
            
        if simi / target_len >= 0.7:
            self.ack[b_idx] = 1 # Mark b as merged
            self.merge(a_idx, b_idx)

    def run(self) -> Dict[str, List[str]]:
        # Sort nodes by degree desc, name asc
        self.nodes.sort(key=lambda x: (-x.m, x.name))
        
        n = len(self.nodes)
        self.ack = [0] * 1000005 # Arbitrary large size from original code, we can use dynamic
        
        final_communities = {} # community_id -> list of node names
        comm_counter = 0
        
        # Progress bar
        pbar = tqdm(total=n, desc="Algorithm 2 Processing")
        
        for i in range(n):
            node = self.nodes[i]
            if node.ex == -1:
                pbar.update(1)
                continue
                
            # Reset ack for current iteration groups
            self.add = []
            
            # Mark neighbors
            for neighbor in node.neighbors:
                if neighbor.ex == -1:
                    continue
                self.ec[neighbor] = 1
            
            # Reset c
            self.c = {}
            
            # DFS
            self.dfs(node, 1, {}, [])
            
            node.ex = -1 # Mark processed
            
            # Compare and merge groups found in this iteration
            self.ack = [0] * len(self.add)
            
            for j in range(len(self.add)):
                if self.ack[j] == 1:
                    continue
                for k in range(j + 1, len(self.add)):
                    if self.ack[k] == 1:
                        continue
                    
                    g1 = self.add[j]
                    g2 = self.add[k]
                    if not g1 or not g2: continue
                    
                    s1 = set(g1)
                    s2 = set(g2)
                    inter = len(s1.intersection(s2))
                    smaller = min(len(g1), len(g2))
                    if smaller > 0 and inter / smaller >= 0.7:
                        # Merge
                        combined = list(s1.union(s2))
                        combined.sort(key=lambda x: x.name)
                        self.add[j] = combined
                        self.add[k] = []
                        self.ack[k] = 1
            
            # Collect "Real groups"
            for j in range(len(self.add)):
                if self.ack[j] == 1:
                    continue
                group = self.add[j]
                if not group: continue
                
                # Filter small communities
                if len(group) < 3:
                    continue

                # Save this group
                comm_id = f"A2_{comm_counter}"
                final_communities[comm_id] = [n.name for n in group]
                comm_counter += 1
                
                # Add virtual node for overlaps
                new_node = Node(node.name) # Same name as current node `sp[i]`
                self.nodes.append(new_node)
                
                for member in group:
                    if self.ec.get(member, 0) == 1:
                        new_node.neighbors.append(member)
                        member.neighbors.append(new_node)
                        new_node.m += 1
                        member.m += 1
            
            # Unmark ec
            for neighbor in node.neighbors:
                self.ec[neighbor] = 0
            
            pbar.update(1)
            
        pbar.close()
        return final_communities


def run_algorithm2(graph: nx.Graph, k1: float = 1, k2: float = 2.3) -> Dict[str, List[str]]:
    algo = Algorithm2(graph, k1, k2)
    return algo.run()
