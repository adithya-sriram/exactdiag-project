#include <complex>
#include <iostream>
#include <Eigen/Dense>
#include "op.h"
#include "model.h"

using namespace std;
using Eigen::MatrixXcd;

namespace pyED {
    typedef Eigen::SparseMatrix<std::complex<double>> ComplexSparseMatrix;
    typedef vector<tuple<complex<double>, int>> arbitrarystate;

    arbitrarystate applyOpToSpinHalfState(Operator op, int rightstate, int sites);
    arbitrarystate applyOpToBosonState(Operator op, int rightstate, RealSpaceBosons model);
    complex<double> getMatrixElementBasisStates(Operator op, int leftstate, int rightstate, int sites);
    complex<double> getInversionStateMatrixElement(Operator op, int state_i, int state_j, int sites, int inversion);
    MatrixXcd constructDiagonalMatrix(Operator op, Model statemap);
    MatrixXcd constructDenseMatrix(Operator op, Model statemap);
    ComplexSparseMatrix constructSparseMatrix(Operator op, Model statemap);
    MatrixXcd constructDenseMatrix1DTranslation(Operator op, Model statemap, double k);
    MatrixXcd constructDenseMatrix1DInversion(Operator op, Model statemap);
    MatrixXcd constructDenseMatrixBoson(Operator op, RealSpaceBosons statemap);
    ComplexSparseMatrix constructSparseMatrixBoson(Operator op, RealSpaceBosons statemap);
}