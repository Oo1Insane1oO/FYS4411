from math import floor, ceil

def L(p,n):
    if n < 0:
        return 0
    return sum([L(i,n-i) for i in range(p)])

def N(n):
    return 1 + sum([N(i) for i in range(int(ceil(n/2.)))]) + L(int(ceil(n/2.-1)),n)

print N(20)
