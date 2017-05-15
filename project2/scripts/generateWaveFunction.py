import sympy as sp

n = [
        [0,0],
        [1,0],
        [0,1],
        [1,1],
        [2,0],
        [0,2],
    ]

def H(x,n):
    """ Use recursive relation """
    if n<0:
        """ ignore negative idices """
        return 0
    elif n==0:
        """ first value """
        return 1
    else:
      return 2*(x*H(x,n-1)-(n-1)*H(x,n-2))
# end function hermite

def phi(alpha,omega,x,y,nx,ny):
    return (H(sp.sqrt(omega)*x,nx)*H(sp.sqrt(omega)*y,ny) *
            sp.exp(-alpha*omega/2*(x**2+y**2)))

def rij(x1,y1,x2,y2):
    return sp.sqrt((x1-x2)**2 + (y1-y2)**2)

def pj(i,j):
    if ((i%2 and j%2)):
        return 1./3
    elif (not (i%2) and not (j%2)):
        return 1./3
    else:
        return 1;

def psiT(alpha,beta,omega,r):
    PhiD = sp.Matrix([[phi(alpha,omega,r[i][0],r[i][1],n[j][0],n[j][1])
        for j in range(len(r)/2)] for i in range(len(r)/2)])
    PhiU = sp.Matrix([[phi(alpha,omega,r[i][0],r[i][1],n[j][0],n[j][1])
        for j in range(len(r)/2)] for i in range(len(r)/2,len(r))])
    g = 0
    for i in range(len(r)):
        for j in range(i+1,len(r)):
            g += pj(i,j) / (beta + 1/rij(r[i][0],r[i][1],r[j][0],r[j][1]))
    return PhiD.det() * PhiU.det() * sp.exp(g)

N = 6
r = [[sp.symbols('x%i' % i),sp.symbols('y%i' % i)] for i in range(1,N+1)]
rsyms = {}
for i in range(1,N+1):
    exec("{0} = sp.sqrt(r[%i-1][0]**2 + r[%i-1][1]**2)".format('r%i' % i) %
            (i,i))
    rsyms[sp.symbols("r%i" % i)] = eval("{0}".format('r%i' % i))

print rsyms

print psiT(1.,0.3,1.,r).subs(rsyms)
