import numpy as np
import sys
import os

alpha = 0
beta = 0
total = 0
potential = 0
kinetic = 0
count = 0
j = 0
directory = sys.argv[1];
for filename in os.listdir(directory):
    data = np.fromfile(directory+filename, dtype=np.float64, count=-1, sep="");

    alpha += data[-2]
    beta += data[-1]
    i = 0
    while i<len(data[:-2]):
        total += data[i]
        potential += data[i+1]
        kinetic += data[i+2]
        i += 3
        j += 1

    count += 1

print total/j
print potential/j
print kinetic/j
print alpha/count
print beta/count
