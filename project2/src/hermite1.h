class HermitePolynomials
{
public:
    static double evaluate(int n, double x, double omega, double alpha);
    static double evaluateDerivative(int n, double x, double omega, double alpha);
    static double evaluateDoubleDerivative(int n, double x, double omega, double alpha);
    static double evaluateAlphaDerivative(int n, double x, double omega, double alpha);
};
