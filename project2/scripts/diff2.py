import sympy as sy
import numpy as np

def makeDict(inputDict):
    returnDict = {}
    for i in inputDict:
        returnDict[i] = sy.symbols(inputDict[i][0], real=True,
                positive=True)**inputDict[i][1]
    return returnDict
# end function makeDict

def collectDenom(expr):
    """ collect input expression in terms of common denominators """
    commons = {}
    for arg in expr.args:
        if sy.denom(arg) in commons:
            commons[sy.denom(arg)] += sy.numer(arg)
        else:
            commons[sy.denom(arg)] = sy.numer(arg)
    return sum([sy.simplify(commons[denom])/denom for denom in commons])
# end function collectDenom

a, b, x1, x2, y1, y2 = sy.symbols('\\alpha \\beta x_1 x_2 y_1 y_2', real=True)
w, d = sy.symbols('\\omega a', real=True, positive=True)

r1sq = x1**2 + y1**2
r2sq = x2**2 + y2**2
dx12 = x1-x2
dy12 = y1-y2
r12sq = dx12**2 + dy12**2

replaceDict = makeDict({r1sq:['r_1',2],r2sq:['r_2',2],r12sq:['r_{12}',2]})

psiT = sy.simplify(sy.exp(a*w/2*(r1sq+r2sq)) *
        sy.exp(d*sy.sqrt(r12sq)/(1+b*sy.sqrt(r12sq))))

# r12 = sy.sqrt(r12sq)
# expr = 1 / (r12**2 * (b*r12+1)**4*(b*(x1-x2)**2+b*(y1-y2)**2+r12))

# print collectDenom(1/x1 + y1/x1 + y2/y1 + 0.5)
 
diff2PsiTx1 = sy.simplify(sy.simplify(sy.diff(psiT, x1, 2)).subs(replaceDict))
diff2PsiTy1 = sy.simplify(sy.simplify(sy.diff(psiT, y1, 2)).subs(replaceDict))
diff2PsiTx2 = sy.simplify(sy.simplify(sy.diff(psiT, x2, 2)).subs(replaceDict))
diff2PsiTy2 = sy.simplify(sy.simplify(sy.diff(psiT, y2, 2)).subs(replaceDict))

total = collectDenom(sy.simplify(0.5*(sy.simplify(-(diff2PsiTx1 + diff2PsiTy1 +
    diff2PsiTx2 + diff2PsiTy2)/psiT) + w**2*(r1sq + r2sq)) + 1/sy.sqrt(r12sq)))

print sy.latex(sy.simplify(sy.nsimplify(sy.simplify(sy.collect(sy.simplify(total.subs(replaceDict)), (b, b**2, b**3)).subs(replaceDict)).subs(replaceDict), rational=True)))
