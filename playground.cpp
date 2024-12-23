#include <iostream>
#include <vector>

using namespace std;

/*
// Get the binary representation of an integer
std::vector<int> getBinaryRepresentation(int integer, int totalBits) {
    std::vector<int> binaryRep(totalBits, 0);
    for (int i = totalBits - 1; i >= 0; i--) {
        binaryRep[i] = (integer >> i) & 1;
    }
    reverse(binaryRep.begin(), binaryRep.end());
    return binaryRep;
}

// Divide binary representation into segments
std::vector<int> divideBinaryIntoSegments(const std::vector<int>& binaryRep, int segmentSize) {
    std::vector<int> segments;
    int currentSegment = 0;

    for (int i = 0; i < binaryRep.size(); i++) {
        currentSegment = (currentSegment << 1) | binaryRep[i];

        if ((i + 1) % segmentSize == 0) {
            segments.push_back(currentSegment);
            currentSegment = 0;
        }
    }

    return segments;
}

// Function to cyclically permute the bits of an integer
int cycleBits(int num, int numBits, int shift) {
    int mask = (1 << numBits) - 1;
    return ((num << shift) | (num >> (numBits - shift))) & mask;
}



int main() {
    int integer = 18;
    int segmentSize = 4;
    int totalBits = segmentSize * 2; // Total bits is a multiple of segmentSize

    std::vector<int> binaryRep = getBinaryRepresentation(integer, totalBits);
    std::vector<int> segments = divideBinaryIntoSegments(binaryRep, segmentSize);
    std::vector<int> concatendatedBinary = getBinaryRepresentation(segments[0], segmentSize);

    for (int i = 1; i < segments.size(); i++) {
        vector<int> vector2 = getBinaryRepresentation(segments[i], segmentSize);
        concatendatedBinary.insert(concatendatedBinary.end(), vector2.begin(), vector2.end());
    }

    for (int i : concatendatedBinary)
        cout << i;

    return 0;
}

// Function to cyclically permute the bits of an integer
int cyclicallyPermuteBits(int num, int numBits, int shift) {
    int mask = (1 << numBits) - 1;
    return ((num << shift) | (num >> (numBits - shift))) & mask;
}

// Function to convert an integer to a cyclically permuted bit string of fixed length
int convertToCyclicBitString(int integer, int totalBits, int shift) {
    return cyclicallyPermuteBits(integer, totalBits, shift);
}

*/

int cycleBits(int num, int numBits, int shift) {
    int mask = (1 << numBits) - 1;
    return ((num << shift) | (num >> (numBits - shift))) & mask;
}

tuple<int, int> findRep1D(int num, int sites) {
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
    return tuple<int, int> (rep, dist);

}

int main() {
    int integer = 5;
    int sites = 3;
    tuple<int, int> rep = findRep1D(integer, sites);
    cout << get<0>(rep) << endl;
    cout << get<1>(rep);

    return 0;
}