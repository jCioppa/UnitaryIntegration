# Unitary Integral Calculator -- 
#
# Timothy Josel Cioppa, Benoit Collins
# 2012 - University of Ottawa
#
# This code can be implemented in Sage and calculates unitary integrals.
#
# --------------------------------------------------------------------------------------------------------
# Functions: 
# 
# TabContent(T,i): 							returns the ith content of the standard tableau T
# ContentVector(T): 						returns the content vector of T
# AxialDistance(i): 						returns the ith axial distance when evaluated on a tableau T
# YoungOrthForm(i,T):						returns Young's Orthogonal Form of s_i = (i,i+1) acting on v_T; NOTE: this is buggy and needs to be fixed
# TabClasses(n):							returns a dictionary of Tableau classes of size n
# MatrixElement(i, n): 						returns the matrix unit in CS_n corresponding so s_i (BUGGY)
# MatrixElementStats(i, n): 				displays matrix statistics about the above matrix
# QS(n): 									the symmetric group algebra Q[S_n]
# QSBasis(n): 								returns a set of basis elements of QS(n) corresponding to permutations
# PairedQSBasis(n): 						returns (perm, basis) pairs, since basis elements of QS(n) are NOT permutation objects
# TabSize(T): 								returns the size of the tableau T
# bar(T): 									returns the tableau with the box [n] removed
# E(T, n): 									central idempotent E_T in CS_n
# trans_perm(T, S, n): 						returns the transition permutation from T to S in S_n as an element of QS(n)
# MatrixUnit(T,S,n): 						returns the matrix unit E_{T,S} in QS(n)
# delta(i, j, k, l): 						dirac delta function
# MUnit(i,j,n): 							matrix unit E_{i,j} in M_n
# InnerProduct(A, B, n): 					inner product in M_n
# tensorProduct(matrixList): 				returns the tensor product of a list of matrices
# indexSet(list): 							returns an index function with values in the list given 
# indices(d, n): 							returns all indices of length d taking values in 1,..,n
# EmSd(p, d, n): 							transforms a permutation p in S_d into an element of M_n^{otimes d}
# GenericMUnit(I,J,d,n):	 				returns the matrix unit E_{i_1, j_1}\otimes ... where I,J are index sets
# basecoefficientpairs(x, n): 				returns pairs of (base, coefficient) for an element x in QS(n)
# algToTensor(p,d,n): 						converts an element p in QS(d) into an element of M_n^\otimes(d)
# UnitaryIntegral(I,J,K,L,d,n): 			given index sets I,J,K,L, of length d, taking values in 1,...,n, this returns the corresponding unitary integral (very slowly)
# --------------------------------------------------------------------------------------------------------

# returns ith content of T
def TabContent(T,i):
    for y in range(len(T)):
        for x in range(len(T[y])):
            if T[y][x] == i: return x-y
    return 'undefined'
        
# returns the content vector of T                               
def ContentVector(T):
    return map(lambda n: TabContent(T,n), range(1, T.size()+1)) 
    
# returns ith axial distance as a function to be evaluated on tableaux     
def AxialDistance(i):
        return lambda T: TabContent(T,i+1) - TabContent(T,i) 


def YoungOrthForm(i, T):       
    if TabContent(T,i+1) == TabContent(T,i)+1:    
        return matrix([[1]])
    elif TabContent(T,i+1) == TabContent(T,i)-1:    
        return matrix([[-1]])
    else:
        return matrix([[1/AxialDistance(i)(T), sqrt(1-1/(AxialDistance(i)(T))^2)],[sqrt(1-  1/(AxialDistance(i)(T))^2), -1/AxialDistance(i)(T)]])

def TabClasses(n): 
    tabs = StandardTableaux(n).list()
    dict = {}
    for T in tabs:
        dict[T.shape()] = T
    return dict    

def MatrixElement(i, n):
    tabs = TabClasses(n).values()
    mats = map(lambda T: YoungOrthForm(i, T), tabs)
    matt = block_diagonal_matrix(mats)
    return matt

def MatrixElementStats(i, n): 
    M = MatrixElement(i,n)
    print 'determinant: ', M.determinant()
    print 'trace: ', M.trace()
    print 'spectrum: ', M.eigenvalues()

def QS(n): return SymmetricGroupAlgebra(QQ, n)
def QSBasis(n): return [QS(n)(p) for p in Permutations(n)]
def PairedQSBasis(n): return [(p, QS(n)(p)) for p in Permutations(n)] 

# ith content of the tableau T
def TabContent(T,i):
    for y in range(len(T)):
        for x in range(len(T[y])):
            if T[y][x] == i: return x-y
    return -999

# size of a tableau T
def TabSize(T):
    S = StandardTableau(copy(T))
    return S.size()

# removes the final box from a tableau T
def bar(T):
    temp = []
    for i in range(len(T)):
        if TabSize(T) in T[i]:            # changed from T.size() to avoid problems
            final = copy(T[i])
            if len(final) >= 2:
                final.pop()
                temp.append(final)  
        else:
            temp.append(copy(T[i]))
    return temp    

# central idempotents corresponding to the tableau T in CS_n
def E(T, n):
    return SymmetricGroupAlgebra(QQ,n).epsilon_ik(T,T)

def E(T,S,d):
    return SymmetricGroupAlgebra(QQ,d).epsilon_ik(T,S)



# permutation tranforming T into S
def trans_perm(T, S, n):
    if not (T.shape() == S.shape()): 
        return 'incompatible shapes'
    else: 
        dict = {1:1}
        for i in range(len(T)):
            for j in range(len(T[i])):
                dict[T[i][j]] = S[i][j]
    perm = []
    for i in dict.keys():
        perm.append(dict[i])
    
    return QS(n)(perm)
    
# matrix units in the symmetric group algebra (without the coefficients, but it doesn't matter)    
def MatrixUnit(T,S,n):
    return SymmetricGroupAlgebra(QQ,n).epsilon_ik(T,S)
   # if not (T.shape() == S.shape()): 
#    return 'incompatible types'
#    else:
#        return E(T,n)*trans_perm(S,T,n)*E(S,n)


# kronecker delta
def delta(i, j, k, l):
    if ((i == k) and (j == l)): return 1
    else: return 0
    
    
# matrix unit E_ij    
def MUnit(i,j,n):
    return matrix(QQ, n, n, lambda k,l: delta(k+1,l+1,i,j))


# matrix inner product
def InnerProduct(A, B, n):
    sum = 0
    for i in range(0,n):
        for j in range(0,n):
            sum+=A[i][j]*B[i][j]
    return sum
    
    
    
# tensor product of a list of matrices    
def tensorProduct(matrixList): 
    d = len(matrixList)   
    if d == 1:
        return matrixList[0]
    if d == 2:
        return matrixList[0].tensor_product(matrixList[1])
    else:    
        return matrixList[0].tensor_product(tensorProduct([matrixList[i] for i in range(1, d)]))


def indexSet(list):
    return lambda k: list[k-1]

#return a list of all index sets [i_1,...,i_d] ranging from 1 to n
def indices(d, n):
    if d<1: 
        return []
    if d == 1:
        return [[i] for i in range(1, n+1)]
    else:
        ret = []
        singletons = indices(1,n)
        temp = indices(d-1,n)
        for list in temp:
            for sing in singletons:
                cpylist = copy(list)
                cpylist.append(sing[0])
                ret.append(cpylist)
        return ret


#embedds S_d into M_n^times(d)
def EmSd(q, d, n):
    p = q.inverse()
    inds = indices(d,n)
    alltensors = []
    for list in inds:  
        it = indexSet(list)
        tensors = [MUnit(it(p[j-1]), it(j), n) for j in range(1,d+1)]  
        alltensors.append(tensors)         
    merp = [tensorProduct(list) for list in alltensors]
    return sum(merp)

#takes a pair of list if indices, like ((i_1,...,i_d), (j_1,...,j_d)) and produces the tensor of matrix units in #M_n(C)
def GenericMUnit(I,J,d,n):
    return tensorProduct([MUnit(I[j], J[j], n) for j in range(d)])

# takes an element x of the group algebra QS(n) and returns a list of (base, coefficient) pairs
def basecoefficientpairs(x, d):
    return [[p, x.coefficient(p)] for p in Permutations(d).list()]    
    
# returns the inner product of elements in the symmetric group algebra S_d    
def SymInnerProduct(x,y,d):
    return sum([x.coefficient(p)*y.coefficient(p) for p in Permutations(d).list()])     

# take an element p of the group algebra QS(d) and return its corresponding element in M_n \otimes ... \otimes #M_n
def isComp(S,d):
    return not (MatrixUnit(S[0], S[1], d) == 'incompatible types')

# embeds QS(d) into M_n \otimes .... \otimes M_n
def algToTensor(p,d,n):
   return sum([x[1]*EmSd(x[0], d, n) for x in basecoefficientpairs(p, d)]) 

def tabPairs(d):
    return [S for S in CartesianProduct(StandardTableaux(d).list(), StandardTableaux(d).list()).list() if S[0].shape() == S[1].shape()] 

def asdf(d,n):
    set = [algToTensor(MatrixUnit(S[0], S[1], d),d,n) for S in tabPairs(d)]     
    return set

# THE FINAL FUNCTION!!! Takes indices I, J, K, L of length d and produces the unitary integral!
def UnitaryIntegral(I,J,K,L,d,n):    
    B = GenericMUnit(J,L,d,n)
    C = GenericMUnit(I,K,d,n)
    set = asdf(d,n)
    return sum([(1/InnerProduct(A,A,n^d))*InnerProduct(B,A,n^d)*InnerProduct(A,C,n^d) for A in set if not InnerProduct(A,A,n^d) == 0])

