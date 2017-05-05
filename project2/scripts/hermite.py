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
        ofile.write("#pragma GCC diagnostic push\n")
        ofile.write('#pragma GCC diagnostic ignored "-Wunused-parameter"\n')
        for c in codes:
            ofile.write(c+"\n")
        # end for c
        ofile.write("#pragma GCC diagnostic pop\n")
        ofile.write("template<typename T> T H(T x, int n) {\n") 
        ofile.write("   if (n > %i) {\n" % (len(codes)-1))
        ofile.write("       return -1;\n")
        ofile.write("   }\n")
        ofile.write("   switch(n) {\n")
        for i in range(len(codes)):
            # end if
            ofile.write("       case %i: return H%i(x);\n" % (i,i))
            # end if
        ofile.write("       default: return 0;\n")
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
