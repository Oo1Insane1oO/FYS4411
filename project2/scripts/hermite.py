##############################################################################
## Script for generating hermite polynomials upto order n (with H_0=1)      ##
##############################################################################

import math
import sys
import sympy as sp

def normal(n):
    """ calculate normal constant """
    return 1./(2**n*math.factorial(n)*math.sqrt(math.pi))
# end function normal

def hermite(x,n):
    """ Use recursive relation """
    if n<0:
        """ ignore negative idices """
        return 0
    elif n==0:
        """ first value """
        return 1
    else:
      return 2*(x*hermite(x,n-1)-(n-1)*hermite(x,n-2))
# end function hermite

def usage():
    """ print usage """
    print "USAGE: python hermite.py 'order' 'filename'\n"
# end function usage

def turnToCPP(n,H):
    """ turn hermite polynomial H of degree n into C++ template function """
    return ''' template<typename T> T H%(n)i(T x) {return %(expr)s;}''' % (
            {'n':n, 'expr':sp.printing.ccode(hermite(x,n))})
# end functtion turnToCPP

def appendToFile(codes, fname):
    """ append polynomial templates to file """
    with open(fname, "w+") as ofile:
        for c in codes:
            ofile.write(c+"\n")
        # end for c
        ofile.write("template<typename T> T H(T x, int n) {\n   switch(n) {\n")
        for i in range(len(codes)):
            ofile.write("       case %i: return H%i(x);\n" % (i,i))
        ofile.write("       default: return H0(x);\n")
        ofile.write("   }\n}")
    # end ofile
# end function appendToFile

if __name__ == "__main__":
    try:
        n = int(sys.argv[1])
        filename = sys.argv[2]
    except IndexError:
        usage()
        raise IndexError
    # end try-except

    x = sp.symbols('x')
    codes = [turnToCPP(i,hermite(x,i)) for i in range(n+1)]
    appendToFile(codes,filename)
# end ifmain
