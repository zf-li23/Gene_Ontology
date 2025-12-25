import sys
import math

class Node:
    def __init__(self, name):
        self.name = name
        self.s = []  # list of int indices
        self.m = 0
        self.ex = 0  # initially 0, set to -1 later

sp = []
T = 0
n = 0
c = []
add = []
ack = []
ec = []
spot = 0
edge = 0
degree = [0] * 105
maxd = 0
k1=1    #0 for loose, 1 for very loose/tight
k2=2.3  #2.3 for very loose, 3 for loose/tight

def cmp_node(a, b):
    if a.m != b.m:
        return a.m > b.m
    return a.name < b.name

def cmp_pair(a, b):
    if a[1] != b[1]:
        return a[1] > b[1]
    return a[0] < b[0]

def in2():
    global sp, T
    with open("input.in", "r") as f:
        T = int(f.readline().strip())
        for _ in range(T):
            line = f.readline().strip().split()
            s1, s2 = line[0], line[1]
            xi = -1
            yi = -1
            for i in range(len(sp)):
                if sp[i].name == s1:
                    xi = i
                if sp[i].name == s2:
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
    global c, add
    c[u] = dep
    res = []
    temp = []
    for v in sp[u].s:
        if sp[v].ex == -1:
            continue
        if c[v] == 0 or c[v] > dep + 1:
            if v not in mp:
                mp[v] = 1
            else:
                mp[v] += 1
    for k, v in mp.items():
        if c[k] == 0 or c[k] > dep + 1:
            res.append((k, v))
    res.sort(key=lambda x: (-x[1], x[0]))
    ck = 0
    for r in res:
        if r[1] >= dep - max(0, int((dep - k1) / k2)):
            ck += 1
            temp.append(r[0])
            del mp[r[0]]
            group.append(r[0])
            dfs(r[0], dep + 1, mp, group)
            group.pop()
            mp[r[0]] = r[1]
        else:
            break
    for t in temp:
        c[t] = 0
    if ck == 0:
        group.sort()
        add.append(group[:])

def merge(a, b):
    global add
    for i in add[a]:
        add[b].append(i)
    add[b].sort()
    add[a] = []
    if not add[b]:
        return
    add[a].append(add[b][0])
    for i in range(1, len(add[b])):
        if add[b][i] != add[b][i - 1]:
            add[a].append(add[b][i])
    add[b] = []

def compare(a, b):
    global add, ack
    simi = 0
    if len(add[a]) < len(add[b]):
        a, b = b, a
    xi = 0
    yi = 0
    ck = 0
    while xi < len(add[a]) and yi < len(add[b]):
        while xi < len(add[a]):
            if add[a][xi] >= add[b][yi]:
                break
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
    with open("iniout.out", "w") as f:
        f.write(str(len(sp)) + '\n')
        for i in range(len(sp)):
            f.write(sp[i].name + ' ' + str(sp[i].m) + '\n')
            for j in sp[i].s:
                if sp[j].ex == -1:
                    continue
                f.write(str(j) + " 1 ")
            f.write('\n')

def out():
    global spot, n
    with open("out2(veryloose).out", "w") as f:
        f.write(str(spot) + '\n')
        for i in range(n, len(sp)):
            f.write(sp[i].name + ' ' + str(sp[i].m) + '\n')
            for j in sp[i].s:
                if sp[j].ex == -1:
                    continue
                f.write(str(j - n) + " 1 ")
            f.write('\n')

def main():
    global sp, n, c, add, ack, ec, spot, edge, degree, maxd
    in2()
    iniout()
    print("Original total spots: " + str(len(sp)))
    sp.sort(key=lambda x: (-x.m, x.name))
    n = len(sp)
    c = [0] * len(sp)
    ack = [0] * 1000005
    ec = [0] * len(sp)
    for i in range(n):
        print("Processing species " + str(i) + ": " + sp[i].name)
        for j in range(len(add)):
            ack[j] = 0
        for j in sp[i].s:
            if sp[j].ex == -1:
                continue
            ec[j] = 1
        add = []
        for j in range(len(sp)):
            c[j] = 0
        dfs(i, 1, {}, [])
        sp[i].ex = -1
        print("    Found " + str(len(add)) + " groups.")
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
            print("      Real group " + str(addsize) + ": ", end="")
            for k in add[j]:
                print(str(k) + "," + sp[k].name + " ", end="")
            print()
            x = len(sp)
            sp.append(Node(sp[i].name))
            for k in add[j]:
                v = k
                if ec[v] == 0:
                    continue
                sp[x].s.append(v)
                sp[x].m += 1
                sp[v].s.append(x)
                sp[v].m += 1
        print("    Found " + str(addsize) + " real groups.")
        print("-------------------------------------------")
        for j in sp[i].s:
            ec[j] = 0
    spot = 0
    edge = 0
    degree = [0] * 105
    maxd = 0
    for i in range(len(sp)):
        if sp[i].ex == -1:
            continue
        spot += 1
        print(sp[i].name + ' ' + str(sp[i].m))
        d = 0
        for j in sp[i].s:
            if sp[j].ex == -1:
                continue
            d += 1
            print("(" + str(j) + "," + sp[j].name + ") ", end="")
        sp[i].m = d
        edge += d
        maxd = max(maxd, d)
        degree[d] += 1
        print()
    edge //= 2
    print("Total spots: " + str(spot) + ", Total relationships: " + str(edge))
    print("Degree distribution:")
    for i in range(1, maxd + 1):
        if degree[i] > 0:
            print("Degree " + str(i) + ": " + str(degree[i]) + " spots")

    out()

if __name__ == "__main__":
    main()