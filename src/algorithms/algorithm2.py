import sys
import math
from collections import defaultdict

class Node:
    def __init__(self, name):
        self.name = name
        self.s = []      # 邻接节点索引
        self.m = 0       # 邻居数量
        self.ex = 0      # 标记是否被处理

sp = []                 # 物种列表
T = 0
n = 0
c = [0] * 10005
s1 = ""
s2 = ""
add = []                # 存储分组
ack = [0] * 1000005
spot = 0
edge = 0
degree = [0] * 105
maxd = 0
global count
count = 0
ec = [0] * 100005
k1=1    #0 for loose, 1 for very loose/tight
k2=3  #2.3 for very loose, 3 for loose/tight

def cmp_node(a, b):
    if a.m != b.m:
        return a.m > b.m
    return a.name < b.name

def cmp_pair(a, b):
    if a[1] != b[1]:
        return a[1] > b[1]
    return a[0] < b[0]

def in2():
    global T, sp
    with open("input.in", "r") as inFile:
        T = int(inFile.readline().strip())
        for _ in range(T):
            s1, s2 = inFile.readline().split()
            xi = -1
            yi = -1
            for i, node in enumerate(sp):
                if node.name == s1:
                    xi = i
                if node.name == s2:
                    yi = i
            if xi == -1:
                xi = len(sp)
                sp.append(Node(s1))
            if yi == -1:
                yi = len(sp)
                sp.append(Node(s2))
            sp[xi].s.append(yi)
            sp[xi].m += 1
            sp[yi].s.append(xi)
            sp[yi].m += 1

def dfs(u, dep, mp, group):
    c[u] = dep
    res = []
    temp = []
    renum = []
    thisgroup = group
    for v in sp[u].s:
        if sp[v].ex == -1:
            continue
        if c[v] == 0 or c[v] > dep + 1:
            renum.append(v)
            if v not in mp:
                mp[v] = 1
            else:
                mp[v] += 1
    for v, cnt in mp.items():
        if c[v] == 0 or c[v] > dep + 1:
            res.append((v, cnt))
    res.sort(key=lambda x: (-x[1], x[0]))
    ck = 0
    for v, cnt in res:
        if cnt >= dep - max(0, int((dep - k1) / k2)):
            ck += 1
            temp.append(v)
            mp.pop(v)
            thisgroup.append(v)
            dfs(v, dep + 1, mp, thisgroup)
            thisgroup.pop()
            mp[v] = cnt
        else:
            break
    for v in renum:
        mp[v] -= 1
        if mp[v] == 0:
            mp.pop(v)
    for v in temp:
        c[v] = 0
    if ck == 0:
        thisgroup.sort()
        add.append(thisgroup[:])

def merge(a, b):
    add[b].extend(add[a])
    add[b].sort()
    add[a].clear()
    if not add[b]:
        return
    add[a].append(add[b][0])
    for i in range(1, len(add[b])):
        if add[b][i] != add[b][i - 1]:
            add[a].append(add[b][i])
    add[b].clear()

def compare(a, b):
    if len(add[a]) < len(add[b]):
        a, b = b, a
    if len(add[a]) == 0:
        ack[b] = 1
        ack[a] = 1
        return
    if len(add[b]) == 0:
        ack[b] = 1
        return
    xi = yi = 0
    simi = 0
    ck = 0
    while xi < len(add[a]) and yi < len(add[b]):
        while xi < len(add[a]) and add[a][xi] < add[b][yi]:
            xi += 1
        if xi >= len(add[a]):
            ck = 1
            break
        if add[a][xi] == add[b][yi]:
            xi += 1
            yi += 1
            simi += 1
        else:
            ck = 1
            break
    if ck == 0:
        ack[b] = 1
    if simi / len(add[b]) >= 0.7:
        ack[b] = 1
        merge(a, b)

def iniout():
    with open("iniout.out", "w") as outFile:
        outFile.write(f"{len(sp)}\n")
        for node in sp:
            outFile.write(f"{node.name} {node.m}\n")
            for v in node.s:
                if sp[v].ex == -1:
                    continue
                outFile.write(f"{v} 1 ")
            outFile.write("\n")


def out():
    global spot, edge, maxd
    with open("out2(tight).out", "w") as outFile:
        outFile.write(f"{spot}\n")
        for i in range(n, len(sp)):
            node = sp[i]
            outFile.write(f"{node.name} {node.m}\n")
            for v in node.s:
                if sp[v].ex == -1:
                    continue
                outFile.write(f"{v - n} 1 ")
            outFile.write("\n")

def main():
    global n, spot, edge, maxd
    in2()
    iniout()
    print(f"Original total spots: {len(sp)}")
    sp.sort(key=lambda x: (-x.m, x.name))
    n = len(sp)
    for i in range(n):
        print(f"Processing species {i}: {sp[i].name}")
        for j in range(len(add)):
            ack[j] = 0
        for v in sp[i].s:
            if sp[v].ex == -1:
                continue
            ec[v] = 1
        add.clear()
        for j in range(len(sp)):
            c[j] = 0
        dfs(i, 1, {}, [])
        sp[i].ex = -1
        print(f"    Found {len(add)} groups.")
        for j in range(len(add)):
            if ack[j] == 1:
                continue
            for k in range(j + 1, len(add)):
                if ack[k] == 1:
                    continue
                compare(j, k)
        addsize = 0
        for j in range(len(add)):
            if ack[j] == 1:
                continue
            addsize += 1
            # print(f"      Real group {addsize}: ", end="")
            # for v in add[j]:
            #     print(f"{v},{sp[v].name} ", end="")
            # print()
            x = len(sp)
            sp.append(Node(sp[i].name))
            for v in add[j]:
                if ec[v] == 0:
                    continue
                sp[x].s.append(v)
                sp[x].m += 1
                sp[v].s.append(x)
                sp[v].m += 1
        print(f"    Found {addsize} real groups.")
        # print("-------------------------------------------")
        for v in sp[i].s:
            ec[v] = 0
    for node in sp:
        if node.ex == -1:
            continue
        spot += 1
        # print(f"{node.name} {node.m}")
        d = 0
        for v in node.s:
            if sp[v].ex == -1:
                continue
            d += 1
            # print(f"({v},{sp[v].name}) ", end="")
        # print()
        node.m = d
        edge += d
        maxd = max(maxd, d)
        degree[d] += 1
    edge //= 2
    print(f"Total spots: {spot}, Total relationships: {edge}")
    print("Degree distribution:")
    for i in range(1, maxd + 1):
        if degree[i] > 0:
            print(f"Degree {i}: {degree[i]} spots")

    out()

if __name__ == "__main__":
    main()