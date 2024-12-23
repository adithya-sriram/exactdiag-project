// Second Quantized Operator Class
// An operator object holds a vector whose elements are themselves vectors. The outer vector holds the sum, and the inner vectors each hold the operator product strings
#include <vector>
#include <complex>

using namespace std;

namespace pyED {

    class SingleSiteOp {
        public:
            bool fermionic;
            // Supported operators: 
            // Spin 1/2 Ops: X(1) Y(2) Z(3) Sp(4) Sm(5)
            // (FUTURE) Fermionic Creation Ops: cD (6) c(7) n(8)
            // Bosonic Creation Ops: b (9) bD(10) nb(11)
            int site;
            int optype;

        SingleSiteOp(int n, int op);
        string to_string();

    };

    class OperatorTerm {
        public:
            vector<SingleSiteOp> terms;
            complex<double> amplitude;

    };

    class Operator {
        public:
            vector<OperatorTerm> ops;
        
        Operator();
        Operator operator+(Operator& o1);
        Operator operator-(Operator& o1);
        Operator operator*(Operator& o1);
        Operator operator*(complex<double> amp1);
        Operator& operator+=(Operator& o1);
        Operator& operator*=(Operator& o1);
        Operator& operator*=(complex<double> amp1);
        Operator power(int n);
        Operator dag();

        string to_string();
    };
    Operator Sx(int n);
    Operator Sy(int n);
    Operator Sz(int n);
    Operator Sp(int n);
    Operator Sm(int n);

    Operator b(int n);
    Operator bD(int n);

}