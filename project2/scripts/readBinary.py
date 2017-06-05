import numpy as np
import sys
import os
# import matplotlib.pyplot as plt

def reader(directory):
    alpha = 0
    beta = 0
    total = 0
    potential = 0
    kinetic = 0
    c = 0
    j = 0
    totalData = []
    for filename in os.listdir(directory):
        data = np.fromfile(directory+filename, dtype=np.float64, count=-1, sep="")

        if str(data[-2])=="nan" or str(data[-1])=="nan":
            continue;
        alpha += data[-2]
        beta += data[-1]
        i = 0
        while i<len(data[:-2]):
            if str(data[i])=="nan" or str(data[i+1])=="nan" or str(data[i+2])=="nan":
                continue;
            total += data[i]
            totalData.append(data[i])
            potential += data[i+1]
            kinetic += data[i+2]
            j += 1
            i += 3
        # end while
        c += 1
    # and for filename
    
    return alpha/c, beta/c, total/j, potential/j, kinetic/j, np.array(totalData)
# end function reader

def blocking(data, numBlocks):
#     numBlocks = numBlocks if (len(data)/numBlocks)%2==0 else
#     int(np.floor(numBlocks))
    minBlockSize = numBlocks/2
    maxBlockSize = len(data)/minBlockSize
    blockStep = int((maxBlockSize-minBlockSize+1)/numBlocks)
    blockMeans = np.zeros(numBlocks)
    for i in range(numBlocks):
        blockSize = int(minBlockSize + i*blockStep)
        tmpBlockNum = len(data)/blockSize
        tmpMean = np.zeros(tmpBlockNum)
        for j in range(len(data)/blockSize):
            tmpMean[j] = np.mean(data[j:(j+1)*blockSize])
        blockMeans[i] = np.sqrt(np.var(tmpMean)) / tmpBlockNum
#         blockMeans[i] = np.sqrt(np.var(data[i:(i+1)*blockSize])) / tmpBlockNum
    # end fori
    return blockMeans
# end function blocking

directory = sys.argv[1];
blockSize = int(sys.argv[2]);
alpha, beta, total, potential, kinetic, totalData = reader(directory)
# means = blocking(totalData, blockSize)

print total
print potential
print kinetic
print potential - kinetic
print alpha
print beta

# plt.plot(np.linspace(0,len(means)/blockSize,len(means)),means)
# plt.show()
