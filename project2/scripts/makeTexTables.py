import numpy as np
import sys

fname = sys.argv[1]
outName = sys.argv[2]
rows = []
with open(fname, 'r') as dataFile:
    data = dataFile.readlines()
    for E in [2.,6.,12.,20.]:
        for w in [1,0.5,0.1,0.05,0.01]:
            for line in data:
                l = line.split()
                if float(l[0]) == E and float(l[1]) == w:
                    for i in xrange(len(l)):
                        if float(l[1])!=1.0 and i==0:
                            l[i] = " "
                        elif float(l[1])==1.0 and i==0:
                            l[i] = str(int(float(l[i])))
                        # end if
                        if i!=3 and i!=0:
                            l[i] = str(np.around(float(l[i]), decimals=5))
                        elif i==3:
                            l[i] = (l[i])[:7] + (l[i])[-4:]
                        # end if
                    if float(l[1])==0.01:
                        rows.append("".join([li+"&" for li in l])[:-1]+"\\\\\\hline")
                    else:
                        rows.append("".join([li+"&" for li in l])[:-1]+"\\\\")
                # end if
            # end for line
        # end forw
    # end forE
# end with open

with open(outName, 'w') as outFile:
    for line in rows:
        outFile.write(line)
    # end for
# end with open
