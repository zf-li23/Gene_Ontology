"""Algorithm 3: DFS-based community detection (Tight variant)."""
from __future__ import annotations

import logging
from typing import Dict, List, Optional, Any
import networkx as nx
from tqdm import tqdm

# Reuse Algorithm2 class but with different default parameters and slightly different logic if needed
# Algorithm 3 has k1=1, k2=3.
# Algorithm 2 has k1=1, k2=2.3.
# The logic in DFS seems identical in structure but let's check `dfs` in algo3.
# Algo3 DFS:
# if (c[r[0]] > dep + 1 or c[r[0]] == 0) and r[1] >= dep - max(0, int((dep - k1) / k2)):
# Algo2 DFS:
# if r[1] >= dep - max(0, int((dep - k1) / k2)):
# Algo3 adds `(c[r[0]] > dep + 1 or c[r[0]] == 0)` check inside the loop condition.
# Also Algo3 has extra logic after loop:
# if dep == 2 and rememadd < len(add): ... (reordering add?)

# So I should implement Algorithm3 as a subclass or separate class.

LOGGER = logging.getLogger(__name__)

class Node:
    def __init__(self, name: str):
        self.name = name
        self.neighbors: List[Node] = []
        self.m = 0
        self.ex = 0

    def __repr__(self):
        return f"Node({self.name}, m={self.m})"

class Algorithm3:
    def __init__(self, graph: nx.Graph, k1: float = 1, k2: float = 3):
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
        rememadd = len(self.add)
        self.c[u] = dep
        res = []
        temp = []
        
        for v in u.neighbors:
            if v.ex == -1:
                continue
            if self.c.get(v, 0) == 0 or self.c.get(v, 0) > dep + 1:
                mp[v] = mp.get(v, 0) + 1
        
        for k, v in mp.items():
            if self.c.get(k, 0) == 0 or self.c.get(k, 0) > dep + 1:
                res.append((k, v))
        
        res.sort(key=lambda x: (-x[1], x[0].name))
        
        ck = 0
        for r in res:
            node_r, count_r = r
            # Algo3 specific condition
            cond1 = (self.c.get(node_r, 0) > dep + 1 or self.c.get(node_r, 0) == 0)
            threshold = dep - max(0, int((dep - self.k1) / self.k2))
            cond2 = (count_r >= threshold)
            
            if cond1 and cond2:
                ck += 1
                temp.append(node_r)
                del mp[node_r]
                group.append(node_r)
                self.dfs(node_r, dep + 1, mp.copy(), group)
                group.pop()
                mp[node_r] = count_r

            # Note: Algo2 had 'else: break'. Algo3 does NOT have 'else: break' in the provided snippet?
            # Let's check provided code.
            # "for r in res: if ...: ... " (no else break visible in snippet provided in context?)
            # Wait, let me check the read_file output for algorithm3.py.
            # Line 84: `for r in res:`
            # Line 85: `if ...:`
            # Line 94: `mp[r[0]] = r[1]`
            # No `else: break`.
            # So Algo3 continues checking other neighbors even if one fails?
            # Yes.
        
        for t in temp:
            self.c[t] = 0
            
        # Algo3 specific logic
        if dep == 2 and rememadd < len(self.add):
            maxaddsize = -1
            maxaddsizeid = -1
            for i in range(rememadd, len(self.add)):
                if len(self.add[i]) > maxaddsize:
                    maxaddsize = len(self.add[i])
                    maxaddsizeid = i
            
            if maxaddsizeid != -1:
                # Swap
                self.add[rememadd], self.add[maxaddsizeid] = self.add[maxaddsizeid], self.add[rememadd]
                # Delete rest?
                # "del add[rememadd + 1:]"
                self.add = self.add[:rememadd + 1]
                
                # Update c for the kept group
                for member in self.add[rememadd]:
                    if self.ec.get(member, 0) == 1:
                        self.c[member] = 2

        if ck == 0:
            new_group = sorted(group[:], key=lambda x: x.name)
            self.add.append(new_group)

    def merge(self, a_idx: int, b_idx: int):
        group_a = self.add[a_idx]
        group_b = self.add[b_idx]
        combined = group_a + group_b
        combined = sorted(list(set(combined)), key=lambda x: x.name)
        self.add[a_idx] = combined
        self.add[b_idx] = []

    def compare(self, a_idx: int, b_idx: int):
        group_a = self.add[a_idx]
        group_b = self.add[b_idx]
        if not group_a or not group_b: return
        
        s1 = set(group_a)
        s2 = set(group_b)
        inter = len(s1.intersection(s2))
        target_len = min(len(group_a), len(group_b))
        
        if target_len > 0 and inter / target_len >= 0.7:
            self.ack[b_idx] = 1
            self.merge(a_idx, b_idx)

    def run(self) -> Dict[str, List[str]]:
        self.nodes.sort(key=lambda x: (-x.m, x.name))
        n = len(self.nodes)
        self.ack = [0] * 1000005
        
        final_communities = {}
        comm_counter = 0
        
        pbar = tqdm(total=n, desc="Algorithm 3 Processing")
        
        for i in range(n):
            node = self.nodes[i]
            if node.ex == -1:
                pbar.update(1)
                continue
                
            self.add = []
            for neighbor in node.neighbors:
                if neighbor.ex == -1: continue
                self.ec[neighbor] = 1
            
            self.c = {}
            self.dfs(node, 1, {}, [])
            
            node.ex = -1
            
            self.ack = [0] * len(self.add)
            for j in range(len(self.add)):
                if self.ack[j] == 1: continue
                for k in range(j + 1, len(self.add)):
                    if self.ack[k] == 1: continue
                    
                    g1 = self.add[j]
                    g2 = self.add[k]
                    if not g1 or not g2: continue
                    
                    s1 = set(g1)
                    s2 = set(g2)
                    inter = len(s1.intersection(s2))
                    smaller = min(len(g1), len(g2))
                    if smaller > 0 and inter / smaller >= 0.7:
                        combined = list(s1.union(s2))
                        combined.sort(key=lambda x: x.name)
                        self.add[j] = combined
                        self.add[k] = []
                        self.ack[k] = 1
            
            for j in range(len(self.add)):
                if self.ack[j] == 1: continue
                group = self.add[j]
                if not group: continue
                
                # Filter small communities
                if len(group) < 3:
                    continue

                comm_id = f"A3_{comm_counter}"
                final_communities[comm_id] = [n.name for n in group]
                comm_counter += 1
                
                new_node = Node(node.name)
                self.nodes.append(new_node)

                for member in group:
                    if self.ec.get(member, 0) == 1:
                        new_node.neighbors.append(member)
                        member.neighbors.append(new_node)
                        new_node.m += 1
                        member.m += 1
            
            for neighbor in node.neighbors:
                self.ec[neighbor] = 0
            
            pbar.update(1)
        
        pbar.close()
        
        # Post-processing: Global Merge of similar communities
        comm_list = []
        for cid, members in final_communities.items():
            comm_list.append(set(members))
            
        merged = True
        while merged:
            merged = False
            new_comm_list = []
            skip_indices = set()
            
            for i in range(len(comm_list)):
                if i in skip_indices: continue
                
                current_set = comm_list[i]
                
                for j in range(i + 1, len(comm_list)):
                    if j in skip_indices: continue
                    
                    other_set = comm_list[j]
                    intersection = len(current_set.intersection(other_set))
                    smaller = min(len(current_set), len(other_set))
                    
                    if smaller > 0 and intersection / smaller >= 0.7:
                        current_set = current_set.union(other_set)
                        skip_indices.add(j)
                        merged = True
                
                new_comm_list.append(current_set)
            
            if merged:
                comm_list = new_comm_list
                
        final_communities = {}
        all_assigned_nodes = set()
        for idx, members in enumerate(comm_list):
            if len(members) < 3: continue
            final_communities[f"A3_{idx}"] = sorted(list(members))
            all_assigned_nodes.update(members)
            
        # Post-processing: Assign unassigned nodes
        graph_nodes = set(str(n) for n in self.graph.nodes())
        unassigned = graph_nodes - all_assigned_nodes
        
        if unassigned:
            adj = {str(n): set(str(nbr) for nbr in self.graph.neighbors(n)) for n in self.graph.nodes()}
            
            for node in unassigned:
                neighbors = adj.get(node, set())
                best_comm = None
                max_overlap = 0
                
                for cid, members in final_communities.items():
                    member_set = set(members)
                    overlap = len(neighbors.intersection(member_set))
                    if overlap > max_overlap:
                        max_overlap = overlap
                        best_comm = cid
                
                if best_comm:
                    final_communities[best_comm].append(node)
                    final_communities[best_comm].sort()

        return final_communities

def run_algorithm3(graph: nx.Graph, k1: float = 1, k2: float = 3) -> Dict[str, List[str]]:
    algo = Algorithm3(graph, k1, k2)
    return algo.run()
