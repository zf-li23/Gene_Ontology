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

def out0.5():
    with open("out_div0.5.out", "w") as f:
        f.write(str(len(sp)) + '\n')
        for i in range(len(sp)):
            f.write(sp[i].name + ' ' + str(sp[i].m) + '\n')
            for j in sp[i].s:
                if sp[j].ex == -1:
                    continue
                f.write(str(j) + " " + str(1 / math.sqrt(sp[i].m)) + " ")
            f.write('\n')

def out2():
    with open("out_div2.out", "w") as f:
        f.write(str(len(sp)) + '\n')
        for i in range(len(sp)):
            f.write(sp[i].name + ' ' + str(sp[i].m) + '\n')
            for j in sp[i].s:
                if sp[j].ex == -1:
                    continue
                f.write(str(j) + " " + str(1 / (sp[i].m*sp[i].m)) + " ")
            f.write('\n')

def out1():
    with open("out_div1.out", "w") as f:
        f.write(str(len(sp)) + '\n')
        for i in range(len(sp)):
            f.write(sp[i].name + ' ' + str(sp[i].m) + '\n')
            for j in sp[i].s:
                if sp[j].ex == -1:
                    continue
                f.write(str(j) + " " + str(1 / sp[i].m) + " ")
            f.write('\n')

def outAverage():
    with open("out_divAverage.out", "w") as f:
        f.write(str(len(sp)) + '\n')
        for i in range(len(sp)):
            f.write(sp[i].name + ' ' + str(sp[i].m) + '\n')
            for j in sp[i].s:
                if sp[j].ex == -1:
                    continue
                f.write(str(j) + " " + str(1 / (math.sqrt(sp[i].m)*math.sqrt(sp[j].m))) + " ")
            f.write('\n')

def main():
    global sp, n, c, add, ack, ec, spot, edge, degree, maxd
    in2()
    iniout()
    out0.5()
    out2()
    out1()
    outAverage()

if __name__ == "__main__":
    main()