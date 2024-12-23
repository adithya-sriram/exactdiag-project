#include <Python.h>
#include <pybind11/pybind11.h>
#include "pybind11/complex.h"
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include "construct_matrix.h"
#include "util.h"

int add(int i, int j) {
    return i + j;
}

namespace py = pybind11;
using namespace std;
using namespace pyED;


PYBIND11_MODULE(pyED, m) {
    m.doc() = "pybind11 example plugin"; // optional module docstring
    m.def("add", &add, "A function that adds two numbers");

    py::class_<Lattice>(m, "Lattice", "Lattice object creates the real space that qubits (qudits) lie on.")
        .def(py::init<int, int, int, complex<double>, complex<double>>(), py::arg("L1"), py::arg("L2") = 1, py::arg("Orbitals") = 1, py::arg("a1") = complex<double> (1,0), py::arg("a2") = complex<double> (0,1))
        .def_readonly("sites", &Lattice::sites)
        .def_readonly("coords", &Lattice::coords)
        .def_readonly("a1", &Lattice::A1)
        .def_readonly("a2", &Lattice::A2)
        .def_readonly("L1", &Lattice::L1)
        .def_readonly("L2", &Lattice::L2)
        .def_readonly("orbitals", &Lattice::orbitals)
        .def_readonly("N", &Lattice::Nsites)
        .def_readonly("klist", &Lattice::klist)
        .def("draw", &Lattice::drawLatticeASCII);

    py::class_<Model>(m, "Model")
        .def_readonly("states", &Model::states)
        .def_readonly("state_indices", &Model::state_indices)
        .def_readonly("num_states", &Model::num_states)
        .def_readonly("N", &Model::N);;

    py::class_<SpinHalf, Model>(m, "SpinHalf", "Spin One Half system on an input lattice.")
        .def(py::init<Lattice, int, int, bool>(), py::arg("Lattice"), py::arg("Nup") = -1, py::arg("inversion") = 0, py::arg("translation") = false)
        .def("printQNs", &SpinHalf::printQNs);

    py::class_<RealSpaceBosons, Model>(m, "RealSpaceBosons", "Bosons on an input lattice.")
        .def(py::init<Lattice, int, bool>(), py::arg("Lattice"), py::arg("N"), py::arg("hardcore"))
        .def_readonly("FockSpaceStates", &RealSpaceBosons::FockSpaceStates);

    py::class_<RealSpaceFermions, Model>(m, "RealSpaceFermions", "Fermions on an input lattice.")
        .def(py::init<Lattice, int, bool, bool>(), py::arg("Lattice"), py::arg("N"), py::arg("number"), py::arg("parity"))
        .def_readonly("FockSpaceStates", &RealSpaceFermions::FockSpaceStates);

    py::class_<SingleSiteOp>(m, "Single Site Operator")
        .def_readonly("site", &SingleSiteOp::site);
    py::class_<OperatorTerm>(m, "Operator Term");
    py::class_<Operator>(m, "Operator", "Operator in Spin Chain or Second Quantization")
        .def(py::init<>())
        .def("__add__", &Operator::operator+)
        .def("__sub__", &Operator::operator-)
        .def("__mul__", [](Operator& self, Operator& o1) { return self.operator*(o1); })  // Operator * Operator
        .def("__mul__", [](Operator& self, std::complex<double> amp1) { return self.operator*(amp1); })  // Operator * complex<double>
        .def("__iadd__", &Operator::operator+=)
        .def("__imul__", [](Operator& self, Operator& o1) { return self.operator*(o1); })  // Operator * Operator
        .def("__imul__", [](Operator& self, std::complex<double> amp1) { return self.operator*(amp1); })  // Operator * complex<double>
        .def("__pow__", &Operator::power)
        //.def("dag", &Operator::dag)
        .def("print", &Operator::to_string);
    
    m.def("Sx", &Sx, "Sx Operator for Spin 1/2");
    m.def("Sy", &Sy, "Sy Operator for Spin 1/2");
    m.def("Sz", &Sz, "Sz Operator for Spin 1/2");
    m.def("Sp", &Sp, "Sp Operator for Spin 1/2");
    m.def("Sm", &Sm, "Sm Operator for Spin 1/2");

    m.def("b", &b, "Annihilation Operator for RealSpaceBosons");
    m.def("bD", &bD, "Creation Operator for RealSpaceBosons");
    

    m.def("constructDenseMatrix", &constructDenseMatrix, "Construct the Dense Matrix of an Operator", py::arg("Operator"), py::arg("Model"));
    m.def("constructSparseMatrix", &constructSparseMatrix, "Construct the Sparse Matrix of an Operator", py::arg("Operator"), py::arg("Model"));
    m.def("constructDiagonalMatrix", &constructDiagonalMatrix, "Construct the Diagonal Matrix of an Operator", py::arg("Operator"), py::arg("Model"));
    m.def("constructDenseMatrix1DTranslation", &constructDenseMatrix1DTranslation, "Construct the Dense Matrix of an Operator in a k sector (only for 1D!)", py::arg("Operator"), py::arg("Model"), py::arg("k"));
    m.def("constructDenseMatrix1DInversion", &constructDenseMatrix1DInversion, "Construct the Dense Matrix of an Inversion Symmetric Operator (only for 1D!)", py::arg("Operator"), py::arg("Model"));
    m.def("constructDenseMatrixBoson", &constructDenseMatrixBoson, "Construct the Dense Matrix of a Boson Operator", py::arg("Operator"), py::arg("RealSpaceBosons"));
    m.def("constructDenseMatrixBoson", &constructDenseMatrixBoson, "Construct the Dense Matrix of a Boson Operator", py::arg("Operator"), py::arg("RealSpaceBosons"));
    m.def("constructSparseMatrixBoson", &constructSparseMatrixBoson, "Construct the Sparse Matrix of a Boson Operator", py::arg("Operator"), py::arg("RealSpaceBosons"));


    
    m.def("findRep1D", &findRep1D, "When there is translational symmetry, given a state, find the representative state and the distance away");    
}
