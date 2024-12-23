// Some useful helpful functions
#include <vector>
#include <complex>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>

using namespace std;

namespace pyED {
    // Pretty prints a complex number
    std::string prettyPrintComplex(const std::complex<double>& z) ;
    // Return the integer of the flipped bitstring
    int Palindrome(int num, int numBits);

    // Sum the bits in the bit string of an int
    int countSetBits(int num);

    // Truncate a small number
    complex<double> chop(complex<double> val);

    // Cyclically permutes the bit string
    int cycleBits(int num, int numBits, int shift);

    // Find equivalence class of state for 1D translation
    int findEC1D(int num, int sites);

    // Given a computational basis state, returns the representative state, the distance to representative state, and the equivalence class numbmer
    tuple<int, int, int> findRep1D(int num, int sites);

    // Get representative states commensurate with a given k for 1D
    vector<int> getCommensurateStates1D(vector<int> states, int sites, double k);

    // Loading bar for for loops
    void printLoadingBar(int progress, int total);

    // Flip the bit at an index
    int flipBit(int num, int index);

    // Get the bit a specific index of a 
    int getBit(int num, int index);

    // Get all partitions of an int
    vector<vector<int>> generatePartitions(int n, vector<int>& partition, int index);
    vector<vector<int>> getPartitions(int n);

    // Pad a vector with 0s
    vector<int> padWithZeros(const vector<int>& input, int minLength);

    // Generate permutations of vector of ints
    vector<vector<int>> generatePermutations(const vector<int>& input);


}