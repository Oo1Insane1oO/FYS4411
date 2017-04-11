 template<typename T> T H0(T x) {return 1;}
 template<typename T> T H1(T x) {return 2*x;}
 template<typename T> T H2(T x) {return 4*pow(x, 2) - 2;}
 template<typename T> T H3(T x) {return 2*x*(4*pow(x, 2) - 2) - 8*x;}
 template<typename T> T H4(T x) {return -24*pow(x, 2) + 2*x*(2*x*(4*pow(x, 2) - 2) - 8*x) + 12;}
 template<typename T> T H5(T x) {return -16*x*(4*pow(x, 2) - 2) + 2*x*(-24*pow(x, 2) + 2*x*(2*x*(4*pow(x, 2) - 2) - 8*x) + 12) + 64*x;}
 template<typename T> T H6(T x) {return 240*pow(x, 2) - 20*x*(2*x*(4*pow(x, 2) - 2) - 8*x) + 2*x*(-16*x*(4*pow(x, 2) - 2) + 2*x*(-24*pow(x, 2) + 2*x*(2*x*(4*pow(x, 2) - 2) - 8*x) + 12) + 64*x) - 120;}
 template<typename T> T H7(T x) {return 192*x*(4*pow(x, 2) - 2) - 24*x*(-24*pow(x, 2) + 2*x*(2*x*(4*pow(x, 2) - 2) - 8*x) + 12) + 2*x*(240*pow(x, 2) - 20*x*(2*x*(4*pow(x, 2) - 2) - 8*x) + 2*x*(-16*x*(4*pow(x, 2) - 2) + 2*x*(-24*pow(x, 2) + 2*x*(2*x*(4*pow(x, 2) - 2) - 8*x) + 12) + 64*x) - 120) - 768*x;}
 template<typename T> T H8(T x) {return -3360*pow(x, 2) + 280*x*(2*x*(4*pow(x, 2) - 2) - 8*x) - 28*x*(-16*x*(4*pow(x, 2) - 2) + 2*x*(-24*pow(x, 2) + 2*x*(2*x*(4*pow(x, 2) - 2) - 8*x) + 12) + 64*x) + 2*x*(192*x*(4*pow(x, 2) - 2) - 24*x*(-24*pow(x, 2) + 2*x*(2*x*(4*pow(x, 2) - 2) - 8*x) + 12) + 2*x*(240*pow(x, 2) - 20*x*(2*x*(4*pow(x, 2) - 2) - 8*x) + 2*x*(-16*x*(4*pow(x, 2) - 2) + 2*x*(-24*pow(x, 2) + 2*x*(2*x*(4*pow(x, 2) - 2) - 8*x) + 12) + 64*x) - 120) - 768*x) + 1680;}
 template<typename T> T H9(T x) {return -3072*x*(4*pow(x, 2) - 2) + 384*x*(-24*pow(x, 2) + 2*x*(2*x*(4*pow(x, 2) - 2) - 8*x) + 12) - 32*x*(240*pow(x, 2) - 20*x*(2*x*(4*pow(x, 2) - 2) - 8*x) + 2*x*(-16*x*(4*pow(x, 2) - 2) + 2*x*(-24*pow(x, 2) + 2*x*(2*x*(4*pow(x, 2) - 2) - 8*x) + 12) + 64*x) - 120) + 2*x*(-3360*pow(x, 2) + 280*x*(2*x*(4*pow(x, 2) - 2) - 8*x) - 28*x*(-16*x*(4*pow(x, 2) - 2) + 2*x*(-24*pow(x, 2) + 2*x*(2*x*(4*pow(x, 2) - 2) - 8*x) + 12) + 64*x) + 2*x*(192*x*(4*pow(x, 2) - 2) - 24*x*(-24*pow(x, 2) + 2*x*(2*x*(4*pow(x, 2) - 2) - 8*x) + 12) + 2*x*(240*pow(x, 2) - 20*x*(2*x*(4*pow(x, 2) - 2) - 8*x) + 2*x*(-16*x*(4*pow(x, 2) - 2) + 2*x*(-24*pow(x, 2) + 2*x*(2*x*(4*pow(x, 2) - 2) - 8*x) + 12) + 64*x) - 120) - 768*x) + 1680) + 12288*x;}
 template<typename T> T H10(T x) {return 60480*pow(x, 2) - 5040*x*(2*x*(4*pow(x, 2) - 2) - 8*x) + 504*x*(-16*x*(4*pow(x, 2) - 2) + 2*x*(-24*pow(x, 2) + 2*x*(2*x*(4*pow(x, 2) - 2) - 8*x) + 12) + 64*x) - 36*x*(192*x*(4*pow(x, 2) - 2) - 24*x*(-24*pow(x, 2) + 2*x*(2*x*(4*pow(x, 2) - 2) - 8*x) + 12) + 2*x*(240*pow(x, 2) - 20*x*(2*x*(4*pow(x, 2) - 2) - 8*x) + 2*x*(-16*x*(4*pow(x, 2) - 2) + 2*x*(-24*pow(x, 2) + 2*x*(2*x*(4*pow(x, 2) - 2) - 8*x) + 12) + 64*x) - 120) - 768*x) + 2*x*(-3072*x*(4*pow(x, 2) - 2) + 384*x*(-24*pow(x, 2) + 2*x*(2*x*(4*pow(x, 2) - 2) - 8*x) + 12) - 32*x*(240*pow(x, 2) - 20*x*(2*x*(4*pow(x, 2) - 2) - 8*x) + 2*x*(-16*x*(4*pow(x, 2) - 2) + 2*x*(-24*pow(x, 2) + 2*x*(2*x*(4*pow(x, 2) - 2) - 8*x) + 12) + 64*x) - 120) + 2*x*(-3360*pow(x, 2) + 280*x*(2*x*(4*pow(x, 2) - 2) - 8*x) - 28*x*(-16*x*(4*pow(x, 2) - 2) + 2*x*(-24*pow(x, 2) + 2*x*(2*x*(4*pow(x, 2) - 2) - 8*x) + 12) + 64*x) + 2*x*(192*x*(4*pow(x, 2) - 2) - 24*x*(-24*pow(x, 2) + 2*x*(2*x*(4*pow(x, 2) - 2) - 8*x) + 12) + 2*x*(240*pow(x, 2) - 20*x*(2*x*(4*pow(x, 2) - 2) - 8*x) + 2*x*(-16*x*(4*pow(x, 2) - 2) + 2*x*(-24*pow(x, 2) + 2*x*(2*x*(4*pow(x, 2) - 2) - 8*x) + 12) + 64*x) - 120) - 768*x) + 1680) + 12288*x) - 30240;}
template<typename T> T H(T x, int n) {
   switch(n) {
       case 0: return H0(x);
       case 1: return H1(x);
       case 2: return H2(x);
       case 3: return H3(x);
       case 4: return H4(x);
       case 5: return H5(x);
       case 6: return H6(x);
       case 7: return H7(x);
       case 8: return H8(x);
       case 9: return H9(x);
       case 10: return H10(x);
       default: return H0(x);
   }
}