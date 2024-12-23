// Model class. Generates the statespace that we work in. Takes in a lattice as input and a list of symmetries, and then enumerates
// all the allowed states. 
// E.g. SpinHalf 
#include <vector>
#include <map>
#include <unordered_map>
#include "lattice.h"


using namespace std;

namespace pyED {
    
    class Model {
        public:
            int N;
            bool translation = false;
            vector<double> kx_vals;
            vector<double> ky_vals;
            int inversion = 0;
            vector<int> states;
            unordered_map<int, int> state_indices;
            int num_states;
            int modelType;

    };

    class SpinHalf : public Model {
        public:
            int Nup;
            bool SzT;
            Lattice lattice;
            SpinHalf(Lattice lattice, int Nup = -1, int inversion = 0, bool translation = false);
            void printQNs();
            
    };

    vector<int> getRepresentatives1D(int sites, int orbitals = 1);

    class RealSpaceBosons : public Model {
        public:
            Lattice lattice;
            RealSpaceBosons(Lattice lattice, int N, bool hardcore);
            bool hardcore;
            vector<vector<int>> FockSpaceStates;
            void printQNs();
    };

    class RealSpaceFermions : public Model {
        public:
            Lattice lattice;
            RealSpaceFermions(Lattice lattice, int N, bool number = false, bool parity = false);
            vector<vector<int>> FockSpaceStates;
            bool number;
            bool parity;
            void printQNs();
    };

}