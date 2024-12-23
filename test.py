import numpy as np
import unittest
import sys
sys.path.append("/Users/adithyasriram/Research Files/Khemani/repos")
sys.path.append('/Users/adithyasriram/Research Files/Khemani/repos/Exact Diagonalization Project/build')
from SpinLibrary import *
from pyED import *


class TestOperatorAlgebra(unittest.TestCase):

    def test_Sx(self):
        Ham = Sx(0)
        self.assertEqual(Ham.print(), "(1.0 + 0.0i)Sx(0)", "Sx creation failed")
    
    def test_Sy(self):
        Ham = Sy(0)
        self.assertEqual(Ham.print(), "(1.0 + 0.0i)Sy(0)", "Sy creation failed")

    def test_Sz(self):
        Ham = Sz(0)
        self.assertEqual(Ham.print(), "(1.0 + 0.0i)Sz(0)", "Sz creation failed")

    def test_add(self):
        Ham = Sx(0) + Sx(1)
        self.assertEqual(Ham.print(), "(1.0 + 0.0i)Sx(0) + (1.0 + 0.0i)Sx(1)", "Sz creation failed")
    
    def test_multiply(self):
        Ham = Sx(0) * Sx(1)
        self.assertEqual(Ham.print(), "(1.0 + 0.0i)Sx(1) * Sx(0)", "Sz creation failed")

class TestDenseMatrixConstruction(unittest.TestCase):

    def test_PauliX(self):
        # PyEd
        model = SpinHalf(Lattice(1))
        Ham = Sx(0)
        Ham_mat = constructDenseMatrix(Ham, model)

        # Spin Library
        s0_list,sx_list,sy_list,sz_list = gen_s0sxsysz(1)
        Ham_mat_py = (sx_list[0]).toarray()

        self.assertEqual(Ham_mat.tolist(), Ham_mat_py.tolist(), 'Pauli X Failed.')
    
    def test_PauliY(self):
        # PyEd
        model = SpinHalf(Lattice(1))
        Ham = Sy(0)
        Ham_mat = constructDenseMatrix(Ham, model)

        # Spin Library
        s0_list,sx_list,sy_list,sz_list = gen_s0sxsysz(1)
        Ham_mat_py = (sy_list[0]).toarray()

        self.assertEqual(Ham_mat.tolist(), Ham_mat_py.tolist(), 'Pauli Y Failed.')

    def test_PauliZ(self):
        # PyEd
        model = SpinHalf(Lattice(1))
        Ham = Sz(0)
        Ham_mat = constructDenseMatrix(Ham, model)

        # Spin Library
        s0_list,sx_list,sy_list,sz_list = gen_s0sxsysz(1)
        Ham_mat_py = (sz_list[0]).toarray()

        self.assertEqual(Ham_mat.tolist(), Ham_mat_py.tolist(), 'Pauli Z Failed.')

    def test_two_site_ising(self) :
        # PyEd
        model = SpinHalf(Lattice(2))
        Ham = Sz(0)
        Ham *= Sz(1)
        Ham_mat = constructDenseMatrix(Ham, model)

        # Spin Library
        s0_list,sx_list,sy_list,sz_list = gen_s0sxsysz(2)
        Ham_mat_py = (sz_list[0] * sz_list[1]).toarray()

        self.assertEqual(Ham_mat.tolist(), Ham_mat_py.tolist(), 'Two Site Ising Failed.')

    def test_eight_site_ising(self) :
        # PyEd
        model = SpinHalf(Lattice(8))
        Ham = Operator()
        for n in range(model.N - 1):
            Ham += Sz(n) * Sz(n+1) 
        Ham_mat = constructDenseMatrix(Ham, model)

        # Spin Library
        s0_list,sx_list,sy_list,sz_list = gen_s0sxsysz(8)
        Ham_mat_py = sz_list[0] * sz_list[1]
        for i in range(1, 7):
            Ham_mat_py += sz_list[i] * sz_list[i+1]
        Ham_mat_py = Ham_mat_py.toarray()
        
        self.assertTrue((Ham_mat == Ham_mat_py).all(), 'Eight Site Ising Failed')

    def test_eight_site_ising_Sz4(self) :
        # PyEd
        model = SpinHalf(Lattice(8), 4)
        Ham = Operator()
        for n in range(model.N - 1):
            Ham += Sz(n) * Sz(n+1) 
        Ham_mat = constructDenseMatrix(Ham, model)

        # Spin Library
        s0_list,sx_list,sy_list,sz_list = gen_s0sxsysz(8)
        Ham_mat_py = sz_list[0] * sz_list[1]
        for i in range(1, 7):
            Ham_mat_py += sz_list[i] * sz_list[i+1]
        Q = Sz_proj(8)[4]
        Ham_mat_py = Q * Ham_mat_py * Q.T
        Ham_mat_py = Ham_mat_py.toarray()

        self.assertEqual(Ham_mat.tolist(), Ham_mat_py.tolist(), 'Eight Site Ising with Sz = 4 Failed.')

    def test_eight_site_ising_Sz2(self) :
        # PyEd
        model = SpinHalf(Lattice(8), 2)
        Ham = Operator()
        for n in range(model.N - 1):
            Ham += Sz(n) * Sz(n+1) 
        Ham_mat = constructDenseMatrix(Ham, model)

        # Spin Library
        s0_list,sx_list,sy_list,sz_list = gen_s0sxsysz(8)
        Ham_mat_py = sz_list[0] * sz_list[1]
        for i in range(1, 7):
            Ham_mat_py += sz_list[i] * sz_list[i+1]
        Q = Sz_proj(8)[-3]
        Ham_mat_py = Q * Ham_mat_py * Q.T
        Ham_mat_py = Ham_mat_py.toarray()

        self.assertEqual(Ham_mat.tolist(), Ham_mat_py.tolist(), 'Eight Site Ising with Sz = 2 Failed.')

    def test_eight_site_TFIM(self):
        # PyEd
        model = SpinHalf(Lattice(8))
        Ham = Operator()
        for n in range(model.N - 1):
            Ham += Sz(n) * Sz(n+1) + Sx(n)
        Ham += Sx(model.N - 1)
        Ham_mat = constructDenseMatrix(Ham, model)

        # Spin Library
        s0_list,sx_list,sy_list,sz_list = gen_s0sxsysz(8)
        Ham_mat_py = sz_list[0] * sz_list[1] + sx_list[0]
        for i in range(1, 7):
            Ham_mat_py += sz_list[i] * sz_list[i+1] + sx_list[i]
        Ham_mat_py += sx_list[7]
        Ham_mat_py = Ham_mat_py.toarray()
        
        self.assertTrue((Ham_mat == Ham_mat_py).all(), 'Eight Site TFIM Failed')

    def test_eight_site_TFIM_YField(self):
        # PyEd
        model = SpinHalf(Lattice(8))
        Ham = Operator()
        for n in range(model.N - 1):
            Ham += Sz(n) * Sz(n+1) + Sy(n)
        Ham += Sy(model.N - 1)
        Ham_mat = constructDenseMatrix(Ham, model)

        # Spin Library
        s0_list,sx_list,sy_list,sz_list = gen_s0sxsysz(8)
        Ham_mat_py = sz_list[0] * sz_list[1] + sy_list[0]
        for i in range(1, 7):
            Ham_mat_py += sz_list[i] * sz_list[i+1] + sy_list[i]
        Ham_mat_py += sy_list[7]
        Ham_mat_py = Ham_mat_py.toarray()
        
        self.assertTrue((Ham_mat == Ham_mat_py).all(), 'Eight Site TFIM Y Field Failed')

class TestSparseMatrixConstruction(unittest.TestCase):

    def test_PauliX_sparse(self):
        # PyEd
        model = SpinHalf(Lattice(1))
        Ham = Sx(0)
        Ham_mat = constructSparseMatrix(Ham, model).toarray()

        # Spin Library
        s0_list,sx_list,sy_list,sz_list = gen_s0sxsysz(1)
        Ham_mat_py = (sx_list[0]).toarray()

        self.assertTrue((Ham_mat == Ham_mat_py).all(), 'Pauli X Sparse Failed.')
    
    def test_PauliY_sparse(self):
        # PyEd
        model = SpinHalf(Lattice(1))
        Ham = Sy(0)
        Ham_mat = constructSparseMatrix(Ham, model).toarray()

        # Spin Library
        s0_list,sx_list,sy_list,sz_list = gen_s0sxsysz(1)
        Ham_mat_py = (sy_list[0]).toarray()

        self.assertTrue((Ham_mat == Ham_mat_py).all(), 'Pauli Y Sparse Failed')

    def test_PauliZ_sparse(self):
        # PyEd
        model = SpinHalf(Lattice(1))
        Ham = Sz(0)
        Ham_mat = constructSparseMatrix(Ham, model)

        # Spin Library
        s0_list,sx_list,sy_list,sz_list = gen_s0sxsysz(1)
        Ham_mat_py = (sz_list[0]).toarray()


        self.assertTrue((Ham_mat == Ham_mat_py).all(), 'Pauli Z Sparse Failed.')

    def test_two_site_ising_sparse(self) :
        # PyEd
        model = SpinHalf(Lattice(2))
        Ham = Sz(0)
        Ham *= Sz(1)
        Ham_mat = constructSparseMatrix(Ham, model).toarray()

        # Spin Library
        s0_list,sx_list,sy_list,sz_list = gen_s0sxsysz(2)
        Ham_mat_py = (sz_list[0] * sz_list[1]).toarray()

        self.assertTrue((Ham_mat == Ham_mat_py).all(), 'Two Site Ising Sparse Failed.')

    def test_eight_site_ising_sparse(self) :
        # PyEd
        model = SpinHalf(Lattice(8))
        Ham = Operator()
        for n in range(model.N - 1):
            Ham += Sz(n) * Sz(n+1) 
        Ham_mat = constructSparseMatrix(Ham, model).toarray()

        # Spin Library
        s0_list,sx_list,sy_list,sz_list = gen_s0sxsysz(8)
        Ham_mat_py = sz_list[0] * sz_list[1]
        for i in range(1, 7):
            Ham_mat_py += sz_list[i] * sz_list[i+1]
        Ham_mat_py = Ham_mat_py.toarray()
        
        self.assertTrue((Ham_mat == Ham_mat_py).all(), 'Eight Site Ising Failed')

    def test_eight_site_ising_Sz4_sparse(self) :
        # PyEd
        model = SpinHalf(Lattice(8), 4)
        Ham = Operator()
        for n in range(model.N - 1):
            Ham += Sz(n) * Sz(n+1) 
        Ham_mat = constructSparseMatrix(Ham, model).toarray()

        # Spin Library
        s0_list,sx_list,sy_list,sz_list = gen_s0sxsysz(8)
        Ham_mat_py = sz_list[0] * sz_list[1]
        for i in range(1, 7):
            Ham_mat_py += sz_list[i] * sz_list[i+1]
        Q = Sz_proj(8)[4]
        Ham_mat_py = Q * Ham_mat_py * Q.T
        Ham_mat_py = Ham_mat_py.toarray()

        self.assertTrue((Ham_mat == Ham_mat_py).all(), 'Eight Site Ising with Sz = 4 Failed.')

    def test_eight_site_ising_Sz2_sparse(self) :
        # PyEd
        model = SpinHalf(Lattice(8), 2)
        Ham = Operator()
        for n in range(model.N - 1):
            Ham += Sz(n) * Sz(n+1) 
        Ham_mat = constructSparseMatrix(Ham, model).toarray()

        # Spin Library
        s0_list,sx_list,sy_list,sz_list = gen_s0sxsysz(8)
        Ham_mat_py = sz_list[0] * sz_list[1]
        for i in range(1, 7):
            Ham_mat_py += sz_list[i] * sz_list[i+1]
        Q = Sz_proj(8)[-3]
        Ham_mat_py = Q * Ham_mat_py * Q.T
        Ham_mat_py = Ham_mat_py.toarray()

        self.assertTrue((Ham_mat == Ham_mat_py).all(), 'Eight Site Ising with Sz = 2 Failed.')

    def test_eight_site_TFIM_sparse(self):
        # PyEd
        model = SpinHalf(Lattice(8))
        Ham = Operator()
        for n in range(model.N - 1):
            Ham += Sz(n) * Sz(n+1) + Sx(n)
        Ham += Sx(model.N - 1)
        Ham_mat = constructSparseMatrix(Ham, model).toarray()

        # Spin Library
        s0_list,sx_list,sy_list,sz_list = gen_s0sxsysz(8)
        Ham_mat_py = sz_list[0] * sz_list[1] + sx_list[0]
        for i in range(1, 7):
            Ham_mat_py += sz_list[i] * sz_list[i+1] + sx_list[i]
        Ham_mat_py += sx_list[7]
        Ham_mat_py = Ham_mat_py.toarray()
        
        self.assertTrue((Ham_mat == Ham_mat_py).all(), 'Eight Site TFIM Failed')

    def test_eight_site_TFIM_YField_sparse(self):
        # PyEd
        model = SpinHalf(Lattice(8))
        Ham = Operator()
        for n in range(model.N - 1):
            Ham += Sz(n) * Sz(n+1) + Sy(n)
        Ham += Sy(model.N - 1)
        Ham_mat = constructSparseMatrix(Ham, model).toarray()

        # Spin Library
        s0_list,sx_list,sy_list,sz_list = gen_s0sxsysz(8)
        Ham_mat_py = sz_list[0] * sz_list[1] + sy_list[0]
        for i in range(1, 7):
            Ham_mat_py += sz_list[i] * sz_list[i+1] + sy_list[i]
        Ham_mat_py += sy_list[7]
        Ham_mat_py = Ham_mat_py.toarray()
        
        self.assertTrue((Ham_mat == Ham_mat_py).all(), 'Eight Site TFIM Y Field Failed')
    
class TestInversionSymmetry1D(unittest.TestCase):

    def test_two_site_ising_invPlus(self) :
        # PyEd
        model = SpinHalf(Lattice(2), inversion = 1)
        Ham = Sz(0)
        Ham *= Sz(1)
        Ham_mat = constructDenseMatrix1DInversion(Ham, model)

        # Spin Library
        s0_list,sx_list,sy_list,sz_list = gen_s0sxsysz(2)
        Ham_mat_py = (sz_list[0] * sz_list[1]).toarray()
        InvProjs = gen_inv_proj(2, 1/2)
        Q = InvProjs[1]
        Ham_mat_py = Q @ Ham_mat_py @ Q.T

        self.assertTrue(np.allclose(Ham_mat, Ham_mat_py, rtol=1e-05, atol=1e-08))

    def test_two_site_ising_invMinus(self) :
        # PyEd
        model = SpinHalf(Lattice(2), inversion = -1)
        Ham = Sz(0)
        Ham *= Sz(1)
        Ham_mat = constructDenseMatrix1DInversion(Ham, model)

        # Spin Library
        s0_list,sx_list,sy_list,sz_list = gen_s0sxsysz(2)
        Ham_mat_py = (sz_list[0] * sz_list[1]).toarray()
        InvProjs = gen_inv_proj(2, 1/2)
        Q = InvProjs[0]
        Ham_mat_py = Q @ Ham_mat_py @ Q.T

        self.assertTrue(np.allclose(Ham_mat, Ham_mat_py, rtol=1e-05, atol=1e-08))

    def test_eight_site_TFIM_invPlus(self):
        # PyEd
        model = SpinHalf(Lattice(8), inversion = 1)
        Ham = Operator()
        for n in range(model.N - 1):
            Ham += Sz(n) * Sz(n+1) + Sx(n)
        Ham += Sx(model.N - 1)
        Ham_mat = constructDenseMatrix1DInversion(Ham, model)

        # Spin Library
        s0_list,sx_list,sy_list,sz_list = gen_s0sxsysz(8)
        Ham_mat_py = sz_list[0] * sz_list[1] + sx_list[0]
        for i in range(1, 7):
            Ham_mat_py += sz_list[i] * sz_list[i+1] + sx_list[i]
        Ham_mat_py += sx_list[7]
        Ham_mat_py = Ham_mat_py.toarray()
        InvProjs = gen_inv_proj(8, 1/2)
        Q = InvProjs[1]
        Ham_mat_py = Q @ Ham_mat_py @ Q.T
        
        self.assertTrue(np.allclose(Ham_mat, Ham_mat_py, rtol=1e-05, atol=1e-08), 'Eight Site TFIM +1 Inversion Failed')

    def test_eight_site_TFIM_invMinus(self):
        # PyEd
        model = SpinHalf(Lattice(8), inversion = -1)
        Ham = Operator()
        for n in range(model.N - 1):
            Ham += Sz(n) * Sz(n+1) + Sx(n)
        Ham += Sx(model.N - 1)
        Ham_mat = constructDenseMatrix1DInversion(Ham, model)

        # Spin Library
        s0_list,sx_list,sy_list,sz_list = gen_s0sxsysz(8)
        Ham_mat_py = sz_list[0] * sz_list[1] + sx_list[0]
        for i in range(1, 7):
            Ham_mat_py += sz_list[i] * sz_list[i+1] + sx_list[i]
        Ham_mat_py += sx_list[7]
        Ham_mat_py = Ham_mat_py.toarray()
        InvProjs = gen_inv_proj(8, 1/2)
        Q = InvProjs[0]
        Ham_mat_py = Q @ Ham_mat_py @ Q.T
        
        self.assertTrue(np.allclose(Ham_mat, Ham_mat_py, rtol=1e-05, atol=1e-08), 'Eight Site TFIM -1 Inversion Failed')

    def test_eight_site_TFIM_YField_invPlus(self):
        # PyEd
        model = SpinHalf(Lattice(8), inversion = 1)
        Ham = Operator()
        for n in range(model.N - 1):
            Ham += Sz(n) * Sz(n+1) + Sy(n)
        Ham += Sy(model.N - 1)
        Ham_mat = constructDenseMatrix1DInversion(Ham, model)

        # Spin Library
        s0_list,sx_list,sy_list,sz_list = gen_s0sxsysz(8)
        Ham_mat_py = sz_list[0] * sz_list[1] + sy_list[0]
        for i in range(1, 7):
            Ham_mat_py += sz_list[i] * sz_list[i+1] + sy_list[i]
        Ham_mat_py += sy_list[7]
        Ham_mat_py = Ham_mat_py.toarray()
        InvProjs = gen_inv_proj(8, 1/2)
        Q = InvProjs[1]
        Ham_mat_py = Q @ Ham_mat_py @ Q.T
        
        self.assertTrue(np.allclose(Ham_mat, Ham_mat_py, rtol=1e-05, atol=1e-08), 'Eight Site TFIM Y Field +1 Inversion Failed')

    def test_eight_site_TFIM_YField_invMinus(self):
        # PyEd
        model = SpinHalf(Lattice(8), inversion = -1)
        Ham = Operator()
        for n in range(model.N - 1):
            Ham += Sz(n) * Sz(n+1) + Sy(n)
        Ham += Sy(model.N - 1)
        Ham_mat = constructDenseMatrix1DInversion(Ham, model)

        # Spin Library
        s0_list,sx_list,sy_list,sz_list = gen_s0sxsysz(8)
        Ham_mat_py = sz_list[0] * sz_list[1] + sy_list[0]
        for i in range(1, 7):
            Ham_mat_py += sz_list[i] * sz_list[i+1] + sy_list[i]
        Ham_mat_py += sy_list[7]
        Ham_mat_py = Ham_mat_py.toarray()
        InvProjs = gen_inv_proj(8, 1/2)
        Q = InvProjs[0]
        Ham_mat_py = Q @ Ham_mat_py @ Q.T
        
        self.assertTrue(np.allclose(Ham_mat, Ham_mat_py, rtol=1e-05, atol=1e-08), 'Eight Site TFIM Y Field -1 Inversion Failed')

    def test_10_site_MFIM_invPlus(self):
        # PyEd
        model = SpinHalf(Lattice(10), inversion = 1)
        Ham = Operator()
        for n in range(model.N - 1):
            Ham += Sz(n) * Sz(n+1) + Sx(n) + Sz(n)
        Ham += Sx(model.N - 1) + Sz(model.N - 1)
        Ham_mat = constructDenseMatrix1DInversion(Ham, model)

        # Spin Library
        s0_list,sx_list,sy_list,sz_list = gen_s0sxsysz(10)
        Ham_mat_py = sz_list[0] * sz_list[1] + sx_list[0] + sz_list[0]
        for i in range(1, 9):
            Ham_mat_py += sz_list[i] * sz_list[i+1] + sz_list[i] + sx_list[i]
        Ham_mat_py += sx_list[9] + sz_list[9]
        Ham_mat_py = Ham_mat_py.toarray()
        InvProjs = gen_inv_proj(10, 1/2)
        Q = InvProjs[1]
        Ham_mat_py = Q @ Ham_mat_py @ Q.T
        
        self.assertTrue(np.allclose(Ham_mat, Ham_mat_py, rtol=1e-05, atol=1e-08), '10 Site MFIM Inversion +1 Failed')

    def test_10_site_MFIM_invMinus(self):
        # PyEd
        model = SpinHalf(Lattice(10), inversion = -1)
        Ham = Operator()
        for n in range(model.N - 1):
            Ham += Sz(n) * Sz(n+1) + Sx(n) + Sz(n)
        Ham += Sx(model.N - 1) + Sz(model.N - 1)
        Ham_mat = constructDenseMatrix1DInversion(Ham, model)

        # Spin Library
        s0_list,sx_list,sy_list,sz_list = gen_s0sxsysz(10)
        Ham_mat_py = sz_list[0] * sz_list[1] + sx_list[0] + sz_list[0]
        for i in range(1, 9):
            Ham_mat_py += sz_list[i] * sz_list[i+1] + sz_list[i] + sx_list[i]
        Ham_mat_py += sx_list[9] + sz_list[9]
        Ham_mat_py = Ham_mat_py.toarray()
        InvProjs = gen_inv_proj(10, 1/2)
        Q = InvProjs[0]
        Ham_mat_py = Q @ Ham_mat_py @ Q.T
        
        self.assertTrue(np.allclose(Ham_mat, Ham_mat_py, rtol=1e-05, atol=1e-08), '10 Site MFIM Inversion -1 Failed')

# Tests compare the list of eigenvalues across all k sectors to the eigenvalues obtained from Dense matrix diagonalization
class TestTranslationSymmetry1D(unittest.TestCase):

    def test_eight_site_ising_1DTR(self) :
        # PyEd
        lattice = Lattice(8)
        model = SpinHalf(lattice, translation = True)
        Ham = Operator()
        for n in range(model.N):
            Ham += Sz(n) * Sz((n+1) % model.N)

        eigvals = []
        for k in lattice.klist:
            Ham_mat = constructDenseMatrix1DTranslation(Ham, model, k.real)
            eigvals = eigvals + np.linalg.eigvals(Ham_mat).real.tolist()

        eigvals.sort()
        full_eigvals = np.linalg.eigvalsh(constructDenseMatrix(Ham,SpinHalf(Lattice(8))))
        
        self.assertTrue(np.allclose(eigvals, full_eigvals, rtol=1e-05, atol=1e-08), 'Eight Site Ising 1D TR Failed')

    def test_eight_site_TFIM_1DTR(self) :
        # PyEd
        lattice = Lattice(8)
        model = SpinHalf(lattice, translation = True)
        Ham = Operator()
        for n in range(model.N):
            Ham += Sz(n) * Sz((n+1) % model.N) + Sx(n)

        eigvals = []
        for k in lattice.klist:
            Ham_mat = constructDenseMatrix1DTranslation(Ham, model, k.real)
            eigvals = eigvals + np.linalg.eigvals(Ham_mat).real.tolist()

        eigvals.sort()
        full_eigvals = np.linalg.eigvalsh(constructDenseMatrix(Ham,SpinHalf(Lattice(8))))
        
        self.assertTrue(np.allclose(eigvals, full_eigvals, rtol=1e-05, atol=1e-08), 'Eight Site TFIM 1D TR Failed')

    def test_eight_site_MFIM_1DTR(self) :
        # PyEd
        lattice = Lattice(8)
        model = SpinHalf(lattice, translation = True)
        Ham = Operator()
        for n in range(model.N):
            Ham += Sz(n) * Sz((n+1) % model.N) + Sx(n) + Sz(n)

        eigvals = []
        for k in lattice.klist:
            Ham_mat = constructDenseMatrix1DTranslation(Ham, model, k.real)
            eigvals = eigvals + np.linalg.eigvals(Ham_mat).real.tolist()

        eigvals.sort()
        full_eigvals = np.linalg.eigvalsh(constructDenseMatrix(Ham,SpinHalf(Lattice(8))))
        
        self.assertTrue(np.allclose(eigvals, full_eigvals, rtol=1e-05, atol=1e-08), 'Eight Site MFIM 1D TR Failed')

    def test_eight_site_MFIM_1DTR_anisotropic(self) :
        # PyEd
        lattice = Lattice(8)
        model = SpinHalf(lattice, translation = True)
        Ham = Operator()
        for n in range(model.N):
            Ham += Sz(n) * Sz((n+1) % model.N) + Sx(n) * 0.5 + Sz(n) * 0.7

        eigvals = []
        for k in lattice.klist:
            Ham_mat = constructDenseMatrix1DTranslation(Ham, model, k.real)
            eigvals = eigvals + np.linalg.eigvals(Ham_mat).real.tolist()

        eigvals.sort()
        full_eigvals = np.linalg.eigvalsh(constructDenseMatrix(Ham,SpinHalf(Lattice(8))))
        
        self.assertTrue(np.allclose(eigvals, full_eigvals, rtol=1e-05, atol=1e-08), 'Eight Site MFIM Anisotropic 1D TR Failed')
        


if __name__ == '__main__':
    unittest.main()

