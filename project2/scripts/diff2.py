import sympy as sy
import numpy as np

x1, x2, y1, y2 = sy.symbols('x_1 x_2 y_1 y_2', real=True)
a, b, w, d = sy.symbols('\\alpha \\beta \\omega a', real=True, positive=True)

r1sq = x1*x1 + y1*y1
r2sq = x2*x2 + y2*y2
dx12 = x1-x2
dy12 = y1-y2
r12sq = dx12**2 + dy12**2

replaceDict = {
        r1sq:sy.symbols('r^2_1')**2,
        r2sq:sy.symbols('r^2_2')**2,
        dx12:sy.symbols('x_{12}'),
        dy12:sy.symbols('y_{12}'),
        r12sq:sy.symbols('r^2_{12}')**2}

psiT = sy.simplify(sy.exp(a*w/2*(r1sq+r2sq)) *
        sy.exp(d*sy.sqrt(r12sq)/(1+b*sy.sqrt(r12sq))))

diff2PsiTx1 = sy.simplify(sy.diff(psiT, x1, 2))
diff2PsiTy1 = sy.simplify(sy.diff(psiT, y1, 2))
diff2PsiTx2 = sy.simplify(sy.diff(psiT, x2, 2))
diff2PsiTy2 = sy.simplify(sy.diff(psiT, y2, 2))

print
sy.latex(sy.simplify(sy.simplify(0.5*(sy.simplify(-(diff2PsiTx1
    + diff2PsiTy1 + diff2PsiTx2 + diff2PsiTy2))/psiT + w**2*(r1sq + r2sq)) +
    1/sy.sqrt(r12sq)).subs(replaceDict)))
