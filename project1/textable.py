import numpy as np
import os, sys

def findNumIdx(l, start, end):
    start = l.index(start)
    return int(l[start+1:start+l[start:].index(end)])
# end function findNumIdx

def findInfoIdx(num):
    if num == 2:
        return 0
    elif num == 6:
        return 1
    elif num == 12:
        return 2
    elif num == 20:
        return 3
# end function findInfoIdx

def findWriteIdx(num):
    if num == 0:
        return 2
    elif num == 1:
        return 6
    elif num == 2:
        return 12
    elif num == 3:
        return 20
# end function findWriteIdx

def findMidIdx(num, start):
    return num - start
# end function findMidIdx

def reader(directory,start,end):
    info = np.zeros((4,end-start+1,4))
    for fname in os.listdir(directory):
        if fname:
            tmpName = fname
            infoIdx = findInfoIdx(findNumIdx(tmpName, "e", "."))
            Eidx = findMidIdx(findNumIdx(tmpName, "E", "_"),start)
            with open(directory+fname, "r") as f:
                info[infoIdx][Eidx][0] = int(f.readline().split()[-1]) * 0.001
                for line in f.readlines():
                    sym = line.split()
                    if sym[0] == "Hartree-Fock":
                        info[infoIdx][Eidx][1] = int(sym[3])
                        info[infoIdx][Eidx][2] = float(sym[-1])
                    # end if
                    if sym[0] == "HF":
                        info[infoIdx][Eidx][3] = int(sym[-1]) * 0.001
                    # end if
                # end for line
            # end with open f
        # end if fname
    # and for fname
    return info
# end function reader

def trunc(fnum,truncPlace):
    strNum = str(fnum)
    dotIdx = strNum.index(".")

    numDecimals = len(strNum[dotIdx+1:])
    if numDecimals <= truncPlace:
        return strNum + "".join(["0" for i in range(truncPlace-numDecimals)])
    # end if

    truncIdx = dotIdx + truncPlace
    truncp1 = truncIdx + 1
    if truncp1 >= len(strNum):
        return strNum
    else:
        lastNum = int(strNum[truncp1])
        truncNum = int(strNum[truncIdx])
        if lastNum >= 5 and truncNum != 9:
            return strNum[:truncIdx] + str(int(strNum[truncIdx])+1)
        elif lastNum >= 5 and truncNum == 9:
            return strNum[:truncIdx] + "1"
        else:
            return strNum[:truncp1]
        # end ifelifelse
    # end ifelse
# end function trunc

directory = sys.argv[1]
start = int(sys.argv[2])
end = int(sys.argv[3])
infos = reader(directory,start,end)
 
decNum = 5
for g,i in enumerate(infos):
    idx = findWriteIdx(g) 
    with open("text/tables/" + directory[5:8] + "_e"+str(idx) + "Table.txt",
            "w") as wFile:
        """open file for writing"""
        for j,k in enumerate(i):
            wFile.write(str(j+start) + " & " + trunc(k[2],decNum) + " & " +
                    str(int(k[1])) + " \\\\\n")
#             wFile.write(str(j+start) + " & " + str(k[2]) + " & " +
#                     str(int(k[1])) + " \\\\\n")
        # end for jk
    # end with wFile
# end for gi
