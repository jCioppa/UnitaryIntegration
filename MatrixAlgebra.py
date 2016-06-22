from Helpers import *
from sage.all import *

class MatrixAlgebra():

    def __init__(self, n_):
        self.n = n_

    def e(self, i,j):
        return matrix(QQ, self.n, self.n, lambda k,l: delta(k+1,l+1,i,j))


class TensorAlgebra():

    def __init__(self, n_, d_):
        self.n = n_
        self.d = d_

    def e(self, I,J):
        MatAlg = MatrixAlgebra(self.n)
        m_list = [MatAlg.e(I[i], J[i]) for i in range(0,self.d)]
        return tensorProduct(m_list)


