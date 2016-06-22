from Helpers import *
from MatrixAlgebra import *
from sage.all import *

class Perm():

    def __init__(self, p):
        self.size = len(p) 
        self.element = Permutation(p)

    def inverse(self):
        return self.element.inverse()

    def embed(self, n):
        p = self.element.inverse()
        inds = indices(self.size, n)
        Mats = MatrixAlgebra(n)
        alltensors = []
        for list in inds:  
            ind = indexSet(list)
            tensors = [Mats.e(ind(p[j-1]), ind(j)) for j in range(1,self.size+1)]  
            alltensors.append(tensors)         
        all_tensors = [tensorProduct(list) for list in alltensors]
        return sum(all_tensors)


class SymGroup():

    def __init__(self, d):
        self.size = d

    def elements(self):
        return  [Perm(Permutation(p)) for p in SymmetricGroup(self.size).list()]

    def embed(self, p, n):
        q = Perm(p)
        return q.embed(n)

class CS():

    def __init__(self, d):
        self.dimension = fact(d)
        self.size = d
        self.algebra = SymmetricGroupAlgebra(QQ,d)


    def permutationBasis(self):
        return self.algebra.basis().list()


    def e(self, T, S = -1):
        if S == -1:
            return self.algebra.epsilon_ik(T,T)
        else:
            return self.algebra.epsilon_ik(T,S)


    def orthogonalBasis(self, shape = -1):

        if shape == -1:
            return [self.e(T[0], T[1]) for T in CartesianProduct(StandardTableaux(self.size), StandardTableaux(self.size)) if T[0].shape() == T[1].shape()]
        else:
            return [self.e(T[0], T[1]) for T in CartesianProduct(StandardTableaux(self.size), StandardTableaux(self.size)) if T[0].shape() == T[1].shape() == shape]

    def yjm(self, k):
        return self.algebra.jucys_murphy(k)


    def idempotents(self, shape = -1):
        if shape == -1:
            return [self.e(T) for T in StandardTableaux(self.size).list()]
        else:
            return [self.e(T) for T in StandardTableaux(self.size).list() if T.shape() == shape]


    def coefficients(self, x):
        return [[p, x.coefficient(p)] for p in Permutations(self.size).list()]    


    def embed(self, x, n):
        Perms = SymGroup(self.size) 
        return sum([x[1]*Perms.embed(x[0], n) for x in self.coefficients(x)]) 

