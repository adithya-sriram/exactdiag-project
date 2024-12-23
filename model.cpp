#include <vector>
#include <iostream>
#include <string>
#include <algorithm>
#include <unordered_set>
#include "model.h"
#include "util.h"
#include <cmath>
#include <algorithm>
#include <numeric> // for std::accumulate



using namespace std;

namespace pyED {

    vector<int> getRepresentatives1D(int sites, int orbitals) {
        vector<int> representatives;
        std::unordered_set<int> usedStates;

        for (int i = 0; i < static_cast<int>(std::pow(2, sites)); i++) {
            if (usedStates.find(i) != usedStates.end()) {
                continue;
            } else {
                representatives.push_back(i);
                usedStates.insert(i);
                int newState = cycleBits(i, sites, orbitals);
                usedStates.insert(newState);
                while (newState != i) {
                    newState = cycleBits(newState, sites, orbitals);
                    usedStates.insert(newState);
                }
            }
        }
        return representatives;
    }

    SpinHalf::SpinHalf(Lattice lattice, int Nup, int inversion, bool translation) {
        this->lattice = lattice;
        this->N = lattice.Nsites;
        this-> Nup = Nup;
        this->modelType = 0;
        for (complex<double> kval : lattice.klist) {
            kx_vals.push_back(kval.real());
            ky_vals.push_back(kval.imag());
        }

        if (translation & (this->lattice.orbitals > 1))
            throw "Translation symmetry not yet implemented for >1 orbital!";

        // Translation
        if (translation) {
            // 1D
            if (this->lattice.L2 == 1) {
                if (Nup != -1) {
                    vector<int> candidates = getRepresentatives1D(this->N, this->lattice.orbitals);
                    int idx = 0;
                    for (int state : candidates) {
                        if (countSetBits(state) == Nup) {
                            this->states.push_back(state);
                            this->state_indices[state] = idx;
                            idx++;
                        }
                    }
                }
                else {
                    this->states = getRepresentatives1D(this->N, this->lattice.orbitals);
                    for (int idx = 0; idx < this->states.size(); idx++)
                        this->state_indices[this->states[idx]] = idx;
                }
                this->num_states = this->states.size();
            }
        }
        // Conserved Sz
        else if ((Nup != -1) & (inversion == 0)) {
            this->SzT = true;
            int idx = 0;
            for (int i = 0; i < pow(2,lattice.Nsites); i++) {
                if (countSetBits(i) == Nup) {
                    this->states.push_back(i);
                    this->state_indices[i] = idx;
                    idx++;
                }
            }
            this->num_states = this->states.size();
        }       
        // Inversion
        else if (((Nup == -1) & (inversion != 0))) {
            this->inversion = inversion;
            int idx = 0;
            for (int i = 0; i < pow(2, lattice.Nsites); i++) {
                int j = Palindrome(i, lattice.Nsites);
                if (i != j) {
                    if (std::find(this->states.begin(), this->states.end(), j) != this->states.end());
                    else {
                        this->states.push_back(i);
                        this->state_indices[i] = idx;
                        idx++;
                    } 
                }
                else if (inversion == 1) {
                    this->states.push_back(i);
                    this->state_indices[i] = idx;
                    idx++;
                }
            }
            this->num_states = this->states.size();
        }
        
        // Inversion and Conserved Sz
        else if((Nup != -1) & (inversion != 0)) {
            this->SzT = true;
            this->inversion = inversion;
            int idx = 0;
            for (int i = 0; i < pow(2,lattice.Nsites); i++) {
                if (countSetBits(i) == Nup) {
                    int j = Palindrome(i, lattice.Nsites);
                    if (i != j) {
                        if (std::find(this->states.begin(), this->states.end(), j) != this->states.end());
                        else{
                            this->states.push_back(i);
                            this->state_indices[i] = idx;
                            idx++;
                        } 
                    }
                    else if (inversion == 1) {
                        this->states.push_back(i);
                        this->state_indices[i] = idx;
                        idx++;
                    }
                }
            }
            this->num_states = this->states.size();
        }
        
        // Trivial
        else {
            this->num_states = pow(2,lattice.Nsites);
            for (int i = 0; i < this->num_states; i++) {
                this->states.push_back(i);
                this->state_indices[i] = i;
            }
        }
    } 

    void SpinHalf::printQNs() {
        if (this->Nup == -1)
            cout << "Spin Sector: All" << endl;
        else {
            cout << "Spin Sector: ";
            cout << this->Nup << endl;
        }
        cout << "Inversion: ";
        cout << this->inversion << endl;
        cout << "Translation: ";
        cout<< this->translation << endl;
        if (this->translation) {
            cout<< "Allowed kx: ";
            for (double kx : this->kx_vals) {
                cout << kx;
                cout << " ";
            }
            cout << endl;
            cout<< "Allowed kx: ";
            for (double ky : this->ky_vals) {
                cout << ky;
                cout << " ";
            }
        }
    }


    RealSpaceBosons::RealSpaceBosons(Lattice lattice, int N, bool hardcore) {
        this->lattice = lattice;
        this->N = N;
        this->modelType = 2;
        this->hardcore = hardcore;

        if (!hardcore) {
            vector<vector<int>> NumberStates = getPartitions(N);
            // Sort each vector in descending order
            for (vector<int>& vec : NumberStates) {
                sort(vec.begin(), vec.end(), greater<int>());
            }
            
            sort( NumberStates.begin(), NumberStates.end() );
            NumberStates.erase( unique( NumberStates.begin(), NumberStates.end() ), NumberStates.end() );



            for (int i = 0; i < NumberStates.size(); i++) {
                NumberStates[i] = padWithZeros(NumberStates[i], this->lattice.Nsites);
            }

            this->FockSpaceStates = generatePermutations(NumberStates[0]);
            
            for (int i = 0; i < NumberStates.size(); i++) {
                for (auto j : generatePermutations(NumberStates[i]))
                    this->FockSpaceStates.push_back(j);
            }

            sort(this->FockSpaceStates.begin(), this->FockSpaceStates.end(), [](const vector<int>& a, const vector<int>& b) {
                // Compare vectors in descending order
                for (int i = 0; i < min(a.size(), b.size()); ++i) {
                    if (a[i] != b[i]) {
                        return a[i] > b[i];
                    }
                }
                return a.size() > b.size();  // If vectors have the same elements up to this point, longer vector comes first
            });

            sort( this->FockSpaceStates.begin(), this->FockSpaceStates.end() );
            this->FockSpaceStates.erase( unique( this->FockSpaceStates.begin(), this->FockSpaceStates.end() ), this->FockSpaceStates.end() );

            reverse(this->FockSpaceStates.begin(), this->FockSpaceStates.end());

            

            this->num_states = this->FockSpaceStates.size();
            for (int i = 0; i < this->num_states; i++) {
                this->states.push_back(i);
                this->state_indices[i] = i;
            }
        }
        else {

            vector<vector<int>> NumberStates = getPartitions(N);
            // Sort each vector in descending order
            for (vector<int>& vec : NumberStates) {
                sort(vec.begin(), vec.end(), greater<int>());
            }
            
            sort( NumberStates.begin(), NumberStates.end() );
            NumberStates.erase( unique( NumberStates.begin(), NumberStates.end() ), NumberStates.end() );



            for (int i = 0; i < NumberStates.size(); i++) {
                NumberStates[i] = padWithZeros(NumberStates[i], this->lattice.Nsites);
            }


            this->FockSpaceStates = generatePermutations(NumberStates[0]);

            sort(this->FockSpaceStates.begin(), this->FockSpaceStates.end(), [](const vector<int>& a, const vector<int>& b) {
                // Compare vectors in descending order
                for (int i = 0; i < min(a.size(), b.size()); ++i) {
                    if (a[i] != b[i]) {
                        return a[i] > b[i];
                    }
                }
                return a.size() > b.size();  // If vectors have the same elements up to this point, longer vector comes first
            });

            sort( this->FockSpaceStates.begin(), this->FockSpaceStates.end() );
            this->FockSpaceStates.erase( unique( this->FockSpaceStates.begin(), this->FockSpaceStates.end() ), this->FockSpaceStates.end() );

            reverse(this->FockSpaceStates.begin(), this->FockSpaceStates.end());

            this->num_states = this->FockSpaceStates.size();
            for (int i = 0; i < this->num_states; i++) {
                this->states.push_back(i);
                this->state_indices[i] = i;
            }
            
        }

    }


    RealSpaceFermions::RealSpaceFermions(Lattice lattice, int N, bool number, bool parity) {
        this->lattice = lattice;
        this->N = N;
        this->modelType = 2;
        this->number = number;
        this->parity = parity;

        if (number) {
            vector<vector<int>> NumberStates = getPartitions(N);
            // Sort each vector in descending order
            for (vector<int>& vec : NumberStates) {
                sort(vec.begin(), vec.end(), greater<int>());
            }
            
            sort( NumberStates.begin(), NumberStates.end() );
            NumberStates.erase( unique( NumberStates.begin(), NumberStates.end() ), NumberStates.end() );



            for (int i = 0; i < NumberStates.size(); i++) {
                NumberStates[i] = padWithZeros(NumberStates[i], this->lattice.Nsites);
            }


            this->FockSpaceStates = generatePermutations(NumberStates[0]);

            sort(this->FockSpaceStates.begin(), this->FockSpaceStates.end(), [](const vector<int>& a, const vector<int>& b) {
                // Compare vectors in descending order
                for (int i = 0; i < min(a.size(), b.size()); ++i) {
                    if (a[i] != b[i]) {
                        return a[i] > b[i];
                    }
                }
                return a.size() > b.size();  // If vectors have the same elements up to this point, longer vector comes first
            });

            sort( this->FockSpaceStates.begin(), this->FockSpaceStates.end() );
            this->FockSpaceStates.erase( unique( this->FockSpaceStates.begin(), this->FockSpaceStates.end() ), this->FockSpaceStates.end() );

            reverse(this->FockSpaceStates.begin(), this->FockSpaceStates.end());

            this->num_states = this->FockSpaceStates.size();
            for (int i = 0; i < this->num_states; i++) {
                vector<int> numbers = this->FockSpaceStates[i];
                int sum = std::accumulate(numbers.begin(), numbers.end(), 0);
                if (sum == N) {
                    this->states.push_back(i);
                    this->state_indices[i] = i;
                }
            }



        }
        else if (parity) {
            
        }

    }

}