import sympy as sp
import numpy as np

a, b, w, d, x1, x2, y1, y2 = sp.symbols('\\alpha \\beta \\omega a x1 x2 y1 y2')

r1 = sp.sqrt(x1*x1 + y1*y1)
r2 = sp.sqrt(x2*x2 + y2*y2)
r1sq = r1**2
r2sq = r2**2
dx12 = x1-x2
dy12 = y1-y2
r12 = sp.sqrt(dx12**2 + dy12**2)
r12sq = dx12**2 + dy12**2

replaceDict = {r1:sp.symbols('r_1'), 
                r2:sp.symbols('r_2'), 
                r12:sp.symbols('r_{12}'),
                r1sq:sp.symbols('r^2_1'),
                r2sq:sp.symbols('r^2_2'),
                r12sq:sp.symbols('r^2_{12}'),
                dx12:sp.symbols('x_{12}'),
                dy12:sp.symbols('y_{12}')}

psiT = sp.simplify(sp.exp(a*w/2*(r1sq+r2sq))*sp.exp(d*r12/(1+b*r12)))

diff2PsiTx1 = sp.simplify(sp.simplify(sp.diff(psiT, x1, 2)).subs(replaceDict))
diff2PsiTy1 = sp.simplify(sp.simplify(sp.diff(psiT, y1, 2)).subs(replaceDict))
diff2PsiTx2 = sp.simplify(sp.simplify(sp.diff(psiT, x2, 2)).subs(replaceDict))
diff2PsiTy2 = sp.simplify(sp.simplify(sp.diff(psiT, y2, 2)).subs(replaceDict))

print sp.latex(sp.simplify(sp.simplify(0.5*(sp.simplify(-(diff2PsiTx1 +
    diff2PsiTy1 + diff2PsiTx2 + diff2PsiTy2))/psiT + w**2*(r1sq + r2sq)) +
    1/r12).subs(replaceDict)))
