 template<typename T> T H0(T x) {return 1;}
 template<typename T> T H1(T x) {return 2*x;}
 template<typename T> T H2(T x) {return 4*pow(x, 2) - 2;}
 template<typename T> T H3(T x) {return 2*x*(4*pow(x, 2) - 2) - 8*x;}
 template<typename T> T H4(T x) {return -24*pow(x, 2) + 2*x*(2*x*(4*pow(x, 2) - 2) - 8*x) + 12;}
 template<typename T> T H5(T x) {return -16*x*(4*pow(x, 2) - 2) + 2*x*(-24*pow(x, 2) + 2*x*(2*x*(4*pow(x, 2) - 2) - 8*x) + 12) + 64*x;}
template<typename T> T H(T x, int n) {
   switch(n) {
       case 0: return H0(x);
       case 1: return H1(x);
       case 2: return H2(x);
       case 3: return H3(x);
       case 4: return H4(x);
       case 5: return H5(x);
   }
}