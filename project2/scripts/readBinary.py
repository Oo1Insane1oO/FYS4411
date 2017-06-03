import numpy as np
import sys
import os

def reader(directory):
    alpha = 0
    beta = 0
    total = 0
    potential = 0
    kinetic = 0
    c = 0
    j = 0
    for filename in os.listdir(directory):
        data = np.fromfile(directory+filename, dtype=np.float64, count=-1, sep="")

        alpha += data[-2]
        beta += data[-1]
        i = 0
        while i<len(data[:-2]):
            total += data[i]
            potential += data[i+1]
            kinetic += data[i+2]
            i += 3
            j += 1

        c += 1
    return alpha/c, beta/c, total/j, potential/j, kinetic/j

directory = sys.argv[1];
alpha, beta, total, potential, kinetic = reader(directory)

print total
print potential
print kinetic
print potential - kinetic
print alpha
print beta
