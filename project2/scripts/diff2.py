import sympy as sy
import numpy as np

x1, x2, y1, y2 = sy.symbols('x_1 x_2 y_1 y_2', real=True)
a, b, w, d = sy.symbols('\\alpha \\beta \\omega a', real=True, positive=True)

r1sq = x1**2 + y1**2
r2sq = x2**2 + y2**2
dx12 = x1-x2
dy12 = y1-y2
r12sq = dx12**2 + dy12**2

replaceDict = {
        r1sq:sy.symbols('r_1', real=True, positive=True)**2,
        r2sq:sy.symbols('r_2', real=True, positive=True)**2,
        r12sq:sy.symbols('r_{12}', real=True, positive=True)**2}

psiT = sy.simplify(sy.exp(a*w/2*(r1sq+r2sq)) *
        sy.exp(d*sy.sqrt(r12sq)/(1+b*sy.sqrt(r12sq))))
 
diff2PsiTx1 = sy.simplify(sy.simplify(sy.diff(psiT, x1, 2)).subs(replaceDict))
diff2PsiTy1 = sy.simplify(sy.simplify(sy.diff(psiT, y1, 2)).subs(replaceDict))
diff2PsiTx2 = sy.simplify(sy.simplify(sy.diff(psiT, x2, 2)).subs(replaceDict))
diff2PsiTy2 = sy.simplify(sy.simplify(sy.diff(psiT, y2, 2)).subs(replaceDict))

total = sy.simplify(0.5*(sy.simplify(-(diff2PsiTx1 + diff2PsiTy1 + diff2PsiTx2
    + diff2PsiTy2)/psiT) + w**2*(r1sq + r2sq)) +
    1/sy.sqrt(r12sq)).collect(2*b**3*d/(b*sp.sqrt(r12sq) +
        1)**4).collect(2*b*a/(b*sy.sqrt(r12sq) + 1)**2)
 
print sy.latex(sy.simplify(total.subs(replaceDict)))
