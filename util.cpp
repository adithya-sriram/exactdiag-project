// Some useful helpful functions
#include <vector>
#include <complex>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <algorithm>
#include <cmath>


using namespace std;

namespace pyED {

    std::string prettyPrintComplex(const std::complex<double>& z) {
        double roundedReal = std::round(z.real() * 10) / 10.0;
        double roundedImag = std::round(z.imag() * 10) / 10.0;

        // Handle negative values correctly
        std::ostringstream realStream, imagStream;
        if (roundedReal >= 0) {
            realStream << std::fixed << std::setprecision(1) << roundedReal;
        } else {
            realStream << std::fixed << std::setprecision(1) << "-" << std::abs(roundedReal);
        }

        if (roundedImag >= 0) {
            imagStream << std::fixed << std::setprecision(1) << roundedImag;
        } else {
            imagStream << std::fixed << std::setprecision(1) << "-" << std::abs(roundedImag);
        }

        std::string roundedRealStr = realStream.str();
        std::string roundedImagStr = imagStream.str();

        return roundedRealStr + " + " + roundedImagStr + "i";
    }

    int Palindrome(int num, int numBits) {
        int reversed = 0;
        int original = num;
        
        for (int i = 0; i < numBits; ++i) {
            reversed = (reversed << 1) | (num & 1);
            num >>= 1;
        }
        return reversed;
    }

    int countSetBits(int num) {
        int count = 0;
        while (num) {
            num &= (num - 1);
            count++;
        }
        return count;
    }

    complex<double> chop(complex<double> val) {
        double val_real = val.real();
        double val_imag = val.imag();
        if (val_real < 1e-10) val_real = 0;
        if (val_imag < 1e-10) val_imag = 0;
        return complex<double> (val_real, val_imag);
    }

    int cycleBits(int num, int numBits, int shift) {
        int mask = (1 << numBits) - 1;
        return ((num << shift) | (num >> (numBits - shift))) & mask;
    }

    int findEC1D(int num, int sites) {
        int newstate = cycleBits(num, sites, 1);
        int cycleLength = 0;
        while (newstate !=  num) {
            cycleLength += 1;
            newstate = cycleBits(newstate, sites, 1);
        }
        return cycleLength + 1;
    }

    tuple<int, int, int> findRep1D(int num, int sites) {
        int rep = num;
        int newstate = num;
        int dist = 0;
        for (int i = 1; i < sites; i++) {
            newstate = cycleBits(newstate, sites, 1);
            if (newstate < rep) {
                rep = newstate;
                dist = i;
            }
        }
        return tuple<int, int, int> (rep, dist, findEC1D(num, sites));

    }

    vector<int> getCommensurateStates1D(vector<int> states, int sites, double k) {
        vector<int> commensurate_states;
        for (int state : states) {
            int state_eq = findEC1D(state, sites);
            double val = (k * state_eq) / (2 * M_PI);
            double intPart;
            double fracPart = modf(val, &intPart);
            double epsilon = 1e-10;

            if (fabs(fracPart) < epsilon) {
                commensurate_states.push_back(state);
            }

        }
        return commensurate_states;

    }

    void printLoadingBar(int progress, int total) {
        const int barWidth = 40;
        float ratio = static_cast<float>(progress) / total;
        int barCount = static_cast<int>(ratio * barWidth);

        std::cout << "[";
        for (int i = 0; i < barWidth; ++i) {
            if (i < barCount) {
                std::cout << "=";
            } else {
                std::cout << " ";
            }
        }
        std::cout << "] " << static_cast<int>(ratio * 100.0) << "%" << std::flush;
    }

    int flipBit(int num, int index) {
        // Create a mask with a 1 at the desired index
        int mask = 1 << index;
        
        // XOR the number with the mask to flip the bit
        return num ^ mask;
    }

    int getBit(int num, int index) {
        // Shift the bit at the desired index to the least significant position
        int shiftedBit = (num >> index) & 1;
        return shiftedBit;
    }

    
    vector<vector<int>> generatePartitions(int n, vector<int>& partition, int index) {
        vector<vector<int>> result;

        if (n == 0) {
            // Add the current partition to the result
            result.push_back(vector<int>(partition.begin(), partition.begin() + index));
            return result;
        }

        for (int i = 1; i <= n; ++i) {
            // Add the current integer to the partition
            partition[index] = i;

            // Recursively generate partitions for the remaining sum
            vector<vector<int>> subPartitions = generatePartitions(n - i, partition, index + 1);

            // Append subpartitions to the result
            result.insert(result.end(), subPartitions.begin(), subPartitions.end());
        }

        return result;
    };

    vector<vector<int>> getPartitions(int n) {
        vector<int> partition(n);
        return generatePartitions(n, partition, 0);
    };

    vector<int> padWithZeros(const vector<int>& input, int minLength) {
        vector<int> paddedVector(minLength, 0);  // Initialize with zeros

        // Copy the values from the input vector to the padded vector
        for (size_t i = 0; i < input.size() && i < static_cast<size_t>(minLength); ++i) {
            paddedVector[i] = input[i];
        }

        return paddedVector;
    }


    vector<vector<int>> generatePermutations(const vector<int>& input) {
        vector<vector<int>> permutations;

        // Sort the input vector to ensure permutations are generated in lexicographically sorted order
        vector<int> sortedInput = input;
        sort(sortedInput.begin(), sortedInput.end());

        do {
            permutations.push_back(sortedInput);
        } while (next_permutation(sortedInput.begin(), sortedInput.end()));

        return permutations;
    }



}
