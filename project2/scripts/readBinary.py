import numpy as np
import sys
import os
import matplotlib.pyplot as plt

def addDecimal(x):
    return x + np.modf(x-int(x))[0]
# end function addDecimal

def averageFlatArea(A, eps=1e-10):
    idx = np.where(np.abs(np.gradient(A, eps))<eps)
    if not np.size(idx) and eps<=1:
        return averageFlatArea(A, addDecimal(eps))
    elif np.size(idx)==0 and eps>1:
        return np.array([-2,-1])
    else:
        return idx
    # end ifeifelse
# end function averateFlatArea

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

#         if str(data[-2])=="nan" or str(data[-1])=="nan":
#             continue;
        alpha += data[-2]
        beta += data[-1]
        i = 0
        while i<len(data[:-2]):
#             if str(data[i])=="nan" or str(data[i+1])=="nan" or str(data[i+2])=="nan":
#                 continue;
            totalData.append(data[i])
            total += data[i]
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
    numBlocks = numBlocks if (len(data)/numBlocks)%2==0 else int(np.int(numBlocks))
    minBlockSize = numBlocks / 2
    maxBlockSize = len(data) / minBlockSize
    blockStep = int((maxBlockSize - minBlockSize)/numBlocks)
    blockMeans = np.zeros(numBlocks)
    blockSizes = np.zeros(numBlocks, np.int64)
    for i in xrange(numBlocks):
        blockSizes[i] = maxBlockSize - i*blockStep
#         blockSizes[i] = minBlockSize + i*blockStep
        tmpBlockNum = int(np.ceil(len(data)/float(blockSizes[i])))
        tmpMean = np.zeros(tmpBlockNum)
        for j in xrange(tmpBlockNum):
            tmpMean[j] = np.mean(data[j*blockSizes[i]:(j+1)*blockSizes[i]]) / numBlocks
        blockMeans[i] = np.var(tmpMean) / tmpBlockNum
    # end fori
    return blockMeans, blockSizes
# end function blocking

def findDivNum(N, D):
    num = 0
    while N>=D:
        N /= D
        num += 1
    # end while
    return num
# end function findDivNum

def findFileName(string):
    oStart = 0
    oEnd = -1
    NE = ""
    for i,s in enumerate(string):
        if s=="w":
            oStart = i
            for j in xrange(len(string[i:])):
                if string[i+j]=="/":
                    oEnd = i+j
                    NE = string[i+j+1:-1]
                    break
                # end if
            # end forj
            break
        # end if
    # end foris
    return string[oStart:oEnd], NE
# end function findFileName

directory = sys.argv[1];
blockSize = int(sys.argv[2]);
oname, Nname = findFileName(directory)
fname = oname + "_" + Nname
if oname[-1] == "0":
    omega = float(oname[2])
else:
    omega = float(oname[1]+"."+oname[2:])
NE = float(Nname[1:])

alpha, beta, total, potential, kinetic, totalData = reader(directory)
means, meansx = blocking(totalData, blockSize)

flatidx = averageFlatArea(means)
meanMean = np.mean(means[flatidx])
print NE, omega, total, meanMean, potential, kinetic, alpha, beta

# maskOff = np.arange(0,np.size(means),int(np.floor(np.size(means)/10)))
plt.plot(meansx,means, 'b-')
# plt.fill_between(meansx, means-meanMean, means+meanMean, facecolor='r',
#         alpha=0.2)
# plt.errorbar(meansx[maskOff], means[maskOff], xerr=0, yerr=meanMean, fmt='o',
#         capthick=0.5, markersize=2.)
plt.xlabel("Blocksize")
plt.ylabel("std. dev.")
# plt.xscale("log", nonposx='clip')
plt.yscale("log", nonposy='clip')
plt.savefig("text/figures/Blocksize_" + fname + "_WOJ.pdf")
# plt.show()
