#include <complex>
#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <algorithm>
#include "op.h"
#include "model.h"
#include "util.h"

using namespace std;
using Eigen::MatrixXcd;

namespace pyED {
    typedef Eigen::SparseMatrix<std::complex<double>> ComplexSparseMatrix;
    typedef vector<tuple<complex<double>, int>> arbitrarystate;


    arbitrarystate applyOpToSpinHalfState(Operator op, int rightstate, int sites) {
        arbitrarystate finalstate;
        for (int i = 0; i < op.ops.size(); i++) {
            complex<double> amp (1,0);
            amp *= op.ops[i].amplitude;
            int newstate = rightstate;
            //bitset<32> binary(newstate);
            for (int j = op.ops[i].terms.size() - 1; j >= 0; j--) {
                //bitset<32> binary(newstate);
                switch(op.ops[i].terms[j].optype) {
                    case(1) : 
                        newstate = flipBit(newstate, sites - op.ops[i].terms[j].site - 1);
                        //binary.flip(sites - op.ops[i].terms[j].site - 1);
                        //newstate = static_cast<int>(binary.to_ulong());
                        break;
                    case(2) :
                        //if (!binary.test(sites - 1 - op.ops[i].terms[j].site))
                        //    amp *= complex<double>(0,-1);
                        if (getBit(newstate, sites - 1 - op.ops[i].terms[j].site) == 0)
                            amp *= complex<double>(0,1);
                        else amp *= complex<double>(0,-1);
                        newstate = flipBit(newstate, sites - op.ops[i].terms[j].site - 1);
                        //binary.flip(sites - op.ops[i].terms[j].site - 1);
                        //newstate = static_cast<int>(binary.to_ulong());
                        break;
                    case(3) :
                        //if (binary.test(sites - 1 - op.ops[i].terms[j].site)) {
                        //    amp *= complex<double>(-1,0);
                        //}
                        if (getBit(newstate, sites - 1 - op.ops[i].terms[j].site) == 1)
                            amp *= complex<double>(-1,0);
                        break;
                    case(4) :
                        //if (!binary.test(sites - 1 - op.ops[i].terms[j].site)) {
                        //    binary.flip(sites - op.ops[i].terms[j].site - 1);
                        //    //newstate = static_cast<int>(binary.to_ulong());
                        //}
                        if (getBit(newstate, sites - 1 - op.ops[i].terms[j].site) == 0)
                            newstate = flipBit(newstate, sites - op.ops[i].terms[j].site - 1);
                        else amp *= complex<double>(0,0);
                        break;
                    case(5) :
                        //if (binary.test(sites - 1 - op.ops[i].terms[j].site)) {
                        //    binary.flip(sites - op.ops[i].terms[j].site - 1);
                            //newstate = static_cast<int>(binary.to_ulong());
                        //}
                        if (getBit(newstate, sites - 1 - op.ops[i].terms[j].site) == 1)
                            newstate = flipBit(newstate, sites - op.ops[i].terms[j].site - 1);
                        else amp *= complex<double>(0,0);
                        break;
                    }
                if (amp == complex<double>(0,0)) break;
            }
        if (amp != complex<double>(0,0)) {
            finalstate.push_back(tuple<complex<double>, int> (amp, newstate));
        }
        }
        return finalstate;
    }

    arbitrarystate applyOpToBosonState(Operator op, int rightstate, RealSpaceBosons model) {
        arbitrarystate finalstate;
        for (int i = 0; i < op.ops.size(); i++) {
            complex<double> amp (1,0);
            amp *= op.ops[i].amplitude;
            int newstate = rightstate;
            vector<int> FockState = model.FockSpaceStates[newstate];

            
            for (int j = op.ops[i].terms.size() - 1; j >= 0; j--) {
                switch(op.ops[i].terms[j].optype) {
                    case(9) : 
                        if (FockState[op.ops[i].terms[j].site] == 0)
                            amp *= 0;
                        else {
                            amp *= sqrt(FockState[op.ops[i].terms[j].site]);
                            FockState[op.ops[i].terms[j].site] -= 1;
                        }
                        break;
                    case(10) :
                        FockState[op.ops[i].terms[j].site] += 1;
                        amp *= sqrt(FockState[op.ops[i].terms[j].site]);
                        break;
                }
                if (amp == complex<double>(0,0)) break;
            }


            
            auto it = std::find(model.FockSpaceStates.begin(), model.FockSpaceStates.end(), FockState);
            if (it != model.FockSpaceStates.end()) {
                newstate = distance(model.FockSpaceStates.begin(), it);
                if (amp != complex<double>(0,0)) 
                    finalstate.push_back(tuple<complex<double>, int> (amp, newstate));
            }

        }
        return finalstate;
    }


    complex<double> getMatrixElementBasisStates(Operator op, int leftstate, int rightstate, int sites) {
        complex<double> matrixelement (0,0);
        arbitrarystate newstate = applyOpToSpinHalfState(op, rightstate, sites);
        for (int k = 0; k < newstate.size(); k++) {
            if (get<1>(newstate[k]) == leftstate) {
                matrixelement += get<0>(newstate[k]);
            }
        }
        return matrixelement;
    }
    
    /*
    complex<double> getInversionStateMatrixElement(Operator op, int state_i, int state_j, int sites, int inversion) {
        int state_i_flipped = Palindrome(state_i, sites);
        int state_j_flipped = Palindrome(state_j, sites);

        if (inversion == -1) {
            return 0.5 * (getMatrixElementBasisStates(op, state_i, state_j, sites) - 
            getMatrixElementBasisStates(op, state_i, state_j_flipped, sites) - getMatrixElementBasisStates(op, state_i_flipped, state_j, sites)
            + getMatrixElementBasisStates(op, state_i_flipped, state_j_flipped, sites));
        }
        else {
            if (state_i == state_i_flipped) {
                if (state_j == state_j_flipped) {
                    return getMatrixElementBasisStates(op, state_i, state_j, sites);
                }
                else {
                    return (1 / sqrt(2)) * (getMatrixElementBasisStates(op, state_i, state_j, sites) + getMatrixElementBasisStates(op, state_i, state_j_flipped, sites));
                }
            }
            else {
                if (state_j == state_j_flipped) {
                    return (1 / sqrt(2)) * (getMatrixElementBasisStates(op, state_i, state_j, sites) + getMatrixElementBasisStates(op, state_i_flipped, state_j, sites));
                }
                else {
                    return 0.5 * (getMatrixElementBasisStates(op, state_i, state_j, sites) + getMatrixElementBasisStates(op, state_i_flipped, state_j, sites) +
                    getMatrixElementBasisStates(op, state_i, state_j, sites) + getMatrixElementBasisStates(op, state_i, state_j_flipped, sites));
                }
            }

        }
    }
    */

    MatrixXcd constructDiagonalMatrix(Operator op, Model statemap) {
        int dim = statemap.num_states;
        MatrixXcd diag(dim,1);
        for (int i = 0; i < dim; i++) {
            diag(i,0) = getMatrixElementBasisStates(op, statemap.states[i], statemap.states[i], statemap.N);
        }
        return diag;
    }

    // Construct full dense matrix subject to possible spin sector (uses the lookup list)
    MatrixXcd constructDenseMatrix(Operator op, Model statemap) {
        int dim = statemap.num_states;        

        MatrixXcd op_mat = MatrixXcd::Zero(dim, dim);
        for (int i = 0; i < dim; i++) {
            arbitrarystate newstate = applyOpToSpinHalfState(op, statemap.states[i], statemap.N);
            for (int k = 0; k < newstate.size(); k++) {
                /*auto it = find(statemap.states.begin(), statemap.states.end(), get<1>(newstate[k]));
                if (it != statemap.states.end()) {
                    int j = distance(statemap.states.begin(), it);
                    op_mat(j,i) += get<0>(newstate[k]);
                } 
                */
                int j = statemap.state_indices[get<1>(newstate[k])];
                op_mat(j,i) += get<0>(newstate[k]);
            }
        }
        return op_mat;

        
    }

    // Construct dense matrix for a 1D translational space group symmetry
    MatrixXcd constructDenseMatrix1DTranslation(Operator op, Model statemap, double k) {
        vector<int> reducedStates = getCommensurateStates1D(statemap.states, statemap.N, k);
        int dim = reducedStates.size();
        MatrixXcd op_mat = MatrixXcd::Zero(dim, dim);
        // Loop through all representative states
        for (int j = 0; j < dim; j++) {
            int nbar = reducedStates[j];
            int nbar_eq = findEC1D(nbar, statemap.N);
            // Apply Hamiltonian to representative state
            arbitrarystate newstate = applyOpToSpinHalfState(op, nbar, statemap.N);
            // Loop through all the states generated by applying H to the representative state
            for (int b = 0; b < newstate.size(); b++) {
                tuple<int, int, int> repData = findRep1D(get<1>(newstate[b]), statemap.N);
                // Get the representative state, the distance away, and the equivalence class
                int mbar = get<0>(repData);
                int d  = get<1>(repData);
                int mbar_eq = get<2>(repData);
                
                
                // Find mbar in the list of representative states
                auto it = find(reducedStates.begin(), reducedStates.end(), mbar);
                if (it != reducedStates.end()) {
                    int i = distance(reducedStates.begin(), it);
                    op_mat(i,j) += get<0>(newstate[b]) * exp(complex<double>(0, -k * d)) * sqrt(nbar_eq) / sqrt(mbar_eq);
                }
                
                // int i = statemap.state_indices[get<1>(newstate[b])];
                // op_mat(i,j) += get<0>(newstate[b]) * exp(complex<double>(0, -k * d)) * sqrt(nbar_eq) / sqrt(mbar_eq);

            }
        }

        return op_mat;
    }

    // Construct full sparse matrix subject to possible spin sector (uses the lookup list)
    ComplexSparseMatrix constructSparseMatrix(Operator op, Model statemap) {
        int dim = statemap.num_states;
        ComplexSparseMatrix op_mat(dim, dim);
        op_mat.setZero(); 
        for (int i = 0; i < dim; i++) {
            arbitrarystate newstate = applyOpToSpinHalfState(op, statemap.states[i], statemap.N);
            for (int k = 0; k < newstate.size(); k++) {
                /*auto it = find(statemap.states.begin(), statemap.states.end(), get<1>(newstate[k]));
                if (it != statemap.states.end()) {
                    int j = distance(statemap.states.begin(), it);
                    op_mat(j,i) += get<0>(newstate[k]);
                } 
                */
                int j = statemap.state_indices[get<1>(newstate[k])];
                op_mat.coeffRef(j,i) += get<0>(newstate[k]);
            }
        }
        return op_mat;
        
    }

    // Construct dense matrix for a 1D inversion space group symmetry (uses the lookup list)
    MatrixXcd constructDenseMatrix1DInversion(Operator op, Model statemap) {
        int dim = statemap.num_states;
        MatrixXcd op_mat = MatrixXcd::Zero(dim, dim);
        for (int j = 0; j < dim; j++) {
            int nbar = statemap.states[j];
            int nbar_eq = 1;
            if (Palindrome(nbar, statemap.N) != nbar)
                nbar_eq = 2;
            
            arbitrarystate newstate = applyOpToSpinHalfState(op, nbar, statemap.N);
            for (int b = 0; b < newstate.size(); b++) {
                int mbar = get<1>(newstate[b]);
                int mbar_eq = 2;
                int flip = 0;
                if (Palindrome(mbar, statemap.N) < mbar) {
                    mbar = Palindrome(mbar, statemap.N);
                    flip = 1;
                }
                else if(Palindrome(mbar, statemap.N) == mbar)
                    mbar_eq = 1;


                auto it = statemap.state_indices.find(mbar);

                if (it != statemap.state_indices.end()) {
                    // Key "K" exists in the map
                    // Access the associated value using it->second
                    int i = it->second;
                    op_mat(i,j) += get<0>(newstate[b]) * pow(statemap.inversion, flip) * sqrt(nbar_eq) / sqrt(mbar_eq);

                } 
                /*
                // Find mbar in the list of representative states
                auto it = find(statemap.states.begin(), statemap.states.end(), mbar);
                if (it != statemap.states.end()) {
                    int i = distance(statemap.states.begin(), it);
                    op_mat(i,j) += get<0>(newstate[b]) * pow(statemap.inversion, flip) * sqrt(nbar_eq) / sqrt(mbar_eq);
                }
                */
                
            }

        }
        return op_mat;
    }

    // Construct dense matrix for a boson operator
    MatrixXcd constructDenseMatrixBoson(Operator op, RealSpaceBosons statemap) {
        int dim = statemap.num_states;        
        MatrixXcd op_mat = MatrixXcd::Zero(dim, dim);
        for (int i = 0; i < dim; i++) {
            arbitrarystate newstate = applyOpToBosonState(op, statemap.states[i], statemap);
            for (int k = 0; k < newstate.size(); k++) {
                /*auto it = find(statemap.states.begin(), statemap.states.end(), get<1>(newstate[k]));
                if (it != statemap.states.end()) {
                    int j = distance(statemap.states.begin(), it);
                    op_mat(j,i) += get<0>(newstate[k]);
                } 
                */
                int j = statemap.state_indices[get<1>(newstate[k])];
                op_mat(j,i) += get<0>(newstate[k]);
            }
        }
        return op_mat;

        
    }

    ComplexSparseMatrix constructSparseMatrixBoson(Operator op, RealSpaceBosons statemap) {
        int dim = statemap.num_states;
        ComplexSparseMatrix op_mat(dim, dim);
        op_mat.setZero(); 
        for (int i = 0; i < dim; i++) {
            arbitrarystate newstate = applyOpToBosonState(op, statemap.states[i], statemap);
            for (int k = 0; k < newstate.size(); k++) {
                /*auto it = find(statemap.states.begin(), statemap.states.end(), get<1>(newstate[k]));
                if (it != statemap.states.end()) {
                    int j = distance(statemap.states.begin(), it);
                    op_mat(j,i) += get<0>(newstate[k]);
                } 
                */
                int j = statemap.state_indices[get<1>(newstate[k])];
                op_mat.coeffRef(j,i) += get<0>(newstate[k]);
            }
        }
        return op_mat;
        
    }


}