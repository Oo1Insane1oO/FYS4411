import numpy as np
import sys

filename = sys.argv[1];
data = np.fromfile(filename, dtype=np.float64, count=-1, sep="");

newData = np.zeros((len(data),6))

for i in range(len(data)/6):
    newData[i] = data[i:i+6]
