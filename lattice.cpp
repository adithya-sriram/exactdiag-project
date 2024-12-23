#include "lattice.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include <algorithm>
#include <string>
#include <stdlib.h>
#include <math.h>
#include <map>

using namespace std;

namespace pyED {

    Lattice::Lattice() {
        this->L1 = 1;
    }

    Lattice::Lattice(int L1, int L2, int orbitals, complex<double> a1, complex<double> a2) {
        this->L1 = L1;
        this->L2 = L2;
        this->orbitals = orbitals;
        this->A1 = a1;
        this->A2 = a2;
        this->Nsites = L1 * L2 * orbitals;
        int qubit_index = 0;

        for (int i = 0; i < L1; i++) {
            for (int j = 0; j < L2; j++) {
                for (int o = 0; o < orbitals; o++) {
                    double x = ((i * 1.0) * a1 + (j * 1.0) * a2).real();
                    double y = ((i * 1.0) * a1 + (j * 1.0) * a2).imag();

                    sites.insert({make_tuple(i,j,o),  qubit_index});
                    coords.insert({qubit_index, make_tuple(x, y, o)});
                    qubit_index += 1;
                }

                complex<double> k ( (2 * M_PI * i) / L1, (2 * M_PI * j) / L2);
                klist.push_back(k);
            }
        }
    
    }

    void Lattice::drawLatticeASCII() {
        if (this->coords.empty()) {
            std::cerr << "No qubit coordinates to draw the lattice." << std::endl;
            return;
        }

        // Find the minimum and maximum coordinates to determine lattice size
        double min_x = std::get<0>(this->coords.begin()->second);
        double max_x = min_x;
        double min_y = std::get<1>(this->coords.begin()->second);
        double max_y = min_y;

        for (const auto& [_, coord] : this->coords) {
            double x = std::get<0>(coord);
            double y = std::get<1>(coord);
            min_x = std::min(min_x, x);
            max_x = std::max(max_x, x);
            min_y = std::min(min_y, y);
            max_y = std::max(max_y, y);
        }

        // Calculate lattice dimensions
        int width = static_cast<int>(std::round(max_x - min_x) + 1);
        int height = static_cast<int>(std::round(max_y - min_y) + 1);

        // Create a 2D array to represent the lattice and initialize it with spaces
        std::vector<std::vector<char>> lattice(height, std::vector<char>(width, ' '));

        // Mark the qubit coordinates with '+'
        for (const auto& [_, coord] : coords) {
            double x = std::get<0>(coord) - min_x;
            double y = std::get<1>(coord) - min_y;
            lattice[static_cast<int>(std::round(y))][static_cast<int>(std::round(x))] = '+';
        }

        // Draw the lattice in ASCII
        for (int y = height - 1; y >= 0; --y) {
            for (int x = 0; x < width; ++x) {
                std::cout << lattice[y][x] << ' ';
            }
            std::cout << std::endl;
        }
    }

    //

}