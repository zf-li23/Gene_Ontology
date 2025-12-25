import sys
import random
from collections import defaultdict

EPS = 1e-9

# 计算模块度增益
def modularity_gain(node, new_community, graph, communities, total_weight, node_degrees, community_degrees):
    ki_in = 0.0
    ki = node_degrees[node]
    
    # 计算节点与新社区内部的连接权重
    for neighbor, weight in graph[node]:
        if communities[neighbor] == new_community:
            ki_in += weight
    
    sum_tot = community_degrees[new_community]
    return (ki_in - sum_tot * ki / total_weight) / total_weight

# Louvain 算法
def louvain_algorithm(node_names, graph):
    n = len(node_names)
    
    # 计算节点度数和总权重
    node_degrees = [0.0] * n
    total_weight = 0.0
    
    for i in range(n):
        for _, weight in graph[i]:
            node_degrees[i] += weight
        total_weight += node_degrees[i]
    total_weight /= 2.0  # 每条边被计算了两次
    
    # 初始化每个节点为一个社区
    communities = list(range(n))
    
    improvement = True
    iteration = 0
    MAX_ITERATIONS = 100  # 防止无限循环
    
    while improvement and iteration < MAX_ITERATIONS:
        improvement = False
        iteration += 1
        
        # 第一阶段：节点移动
        local_improvement = True
        local_iteration = 0
        MAX_LOCAL_ITERATIONS = 10
        
        while local_improvement and local_iteration < MAX_LOCAL_ITERATIONS:
            local_improvement = False
            local_iteration += 1
            
            # 计算每个社区的度数总和
            community_degrees = [0.0] * n
            for i in range(n):
                community_degrees[communities[i]] += node_degrees[i]
            
            # 随机访问顺序
            node_order = list(range(n))
            random.shuffle(node_order)
            
            # 尝试移动每个节点
            for idx in range(n):
                node = node_order[idx]
                current_community = communities[node]
                
                # 找到节点的所有邻居社区
                neighbor_communities = set()
                for neighbor, _ in graph[node]:
                    neighbor_communities.add(communities[neighbor])
                
                # 尝试移动到每个邻居社区
                best_gain = 0.0
                best_community = current_community
                
                for new_community in neighbor_communities:
                    if new_community == current_community:
                        continue
                    
                    gain = modularity_gain(node, new_community, graph,
                                         communities, total_weight,
                                         node_degrees, community_degrees)
                    
                    if gain > best_gain + EPS:
                        best_gain = gain
                        best_community = new_community
                
                # 如果找到更好的社区，移动节点
                if best_community != current_community and best_gain > EPS:
                    communities[node] = best_community
                    local_improvement = True
                    improvement = True
        
        # 如果没有改进，结束算法
        if not improvement:
            break
        
        # 第二阶段：社区合并
        # 重新标记社区，使编号连续
        community_mapping = {}
        new_community_count = 0
        for i in range(n):
            old_community = communities[i]
            if old_community not in community_mapping:
                community_mapping[old_community] = new_community_count
                new_community_count += 1
            communities[i] = community_mapping[old_community]
        
        # 如果社区数量没有减少，结束算法
        if new_community_count >= n:
            break
        
        # 创建新的图（社区级别的图）
        new_graph = [[] for _ in range(new_community_count)]
        new_names = [""] * new_community_count
        community_nodes = [set() for _ in range(new_community_count)]
        
        # 收集每个社区的节点
        for i in range(n):
            community = communities[i]
            community_nodes[community].add(i)
            if new_names[community] == "":
                new_names[community] = node_names[i]
            else:
                new_names[community] += "," + node_names[i]
        
        # 重新计算社区之间的边权
        for i in range(new_community_count):
            edge_map = defaultdict(float)
            
            for node in community_nodes[i]:
                for neighbor, weight in graph[node]:
                    neighbor_community = communities[neighbor]
                    
                    if neighbor_community == i:
                        # 社区内部边，权重加倍（因为会被计算两次）
                        edge_map[i] += weight
                    else:
                        # 社区间边
                        edge_map[neighbor_community] += weight
            
            # 添加到新图
            for neigh, weight in edge_map.items():
                new_graph[i].append((neigh, weight))
        
        # 使用新图继续迭代（在实际Louvain算法中，这里应该递归调用）
        # 但根据题目要求，我们简化处理，只进行一次划分
        break
    
    # 收集结果
    community_result = defaultdict(set)
    for i in range(n):
        community_result[communities[i]].add(node_names[i])
    
    # 转换为向量并去重
    result = []
    for com_set in community_result.values():
        result.append(com_set)
    
    # 去重：如果两个社区的节点名称集相同，只保留一个
    unique_result = []
    for com_set in result:
        is_duplicate = False
        for existing_set in unique_result:
            if com_set == existing_set:
                is_duplicate = True
                break
        if not is_duplicate:
            unique_result.append(com_set)
    
    return unique_result

def cmp_community(a, b):
    return len(a) > len(b)

def main():
    in_name, out_suffix = input().strip().split()
    in_file = in_name + ".out"
    out_file = "louvain_output" + out_suffix + ".csv"
    
    with open(in_file, "r") as fin:
        n = int(fin.readline().strip())
        
        node_names = []
        graph = []
        
        for i in range(n):
            line = fin.readline().strip().split()
            name = line[0]
            m = int(line[1])
            node_names.append(name)
            edges = []
            for j in range(m):
                neighbor = int(line[2 + 2*j])
                weight = float(line[3 + 2*j])
                edges.append((neighbor, weight))
                
                # 对于无向图，添加反向边（如果输入没有提供的话）
                has_reverse = False
                for neigh, _ in graph[neighbor] if neighbor < len(graph) else []:
                    if neigh == i:
                        has_reverse = True
                        break
                if not has_reverse and neighbor < len(graph):
                    graph[neighbor].append((i, weight))
            graph.append(edges)
    
    # 运行 Louvain 算法
    communities = louvain_algorithm(node_names, graph)
    
    communities.sort(key=lambda x: len(x), reverse=True)
    
    # 输出结果
    with open(out_file, "w") as fout:
        countnum = 0
        for com_set in communities:
            if len(com_set) < 4:
                break
            countnum += 1
        print("Detected " + str(countnum) + " communities:")
        for com_set in communities:
            if len(com_set) < 4:
                break  # 只输出大小不少于4的社区
            print("Community of size " + str(len(com_set)) + ": ")
            for name in com_set:
                print(name + " ", end="")
            print()
            print("-------------------------------------------")
        
        fout.write("community_id,members\n")
        count = 0
        for com_set in communities:
            if len(com_set) < 4:
                break  # 只输出大小不少于4的社区
            count += 1
            fout.write("L" + str(count) + ",")
            for name in com_set:
                fout.write(name + " ")
            fout.write("\n")

if __name__ == "__main__":
    main()