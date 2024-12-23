#include <vector>
#include <complex>
#include <map>
#include <string>
#include <tuple>

using namespace std;

// Lattice Class. Can create a lattice with L1 and/or L2 and orbitals. May also specify Bravais lattice vectors. Lattice vectors are specified using
// complex numbers. 


namespace pyED {

struct Lattice {

    public:
        int Nsites;
        int L1;
        int L2;
        int orbitals;
        complex<double> A1;
        complex<double> A2;
        map<tuple<int, int, int>, int> sites;
        map<int, tuple<double, double, int>> coords;
        vector<complex<double>> klist;
        Lattice();
        Lattice(int L1, int L2 = 1, int orbitals = 1, complex<double> a1 = complex<double>(1,0), complex<double> a2 = complex<double>(0,1));
        void drawLatticeASCII();
        
};


}