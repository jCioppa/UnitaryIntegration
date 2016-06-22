from Helpers import *
from Diagram import *
from SymmetricGroupAlgebra import *

# idempotent E_T in C[S_d] or E_TS if S is given 
# NOT NEEDED
def E(d,T,S = -1):
    Alg = CS(d)
    if S == -1:
        return Alg.e(T,T)
    else:
        return Alg.e(T,S)

# orthogonal basis of C[S_d]
# if we specify a diagram's shape lamb, we only return those basis elements corresponding to that diagram
# NOT NEEDED
def orthogonalBasis(d, lamb = -1):
    if lamb == -1:
        return [E(d,T[0], T[1])  for T in CartesianProduct(StandardTableaux(d).list(), StandardTableaux(d).list()) if T[0].shape() == T[1].shape()]
    else:
        return [E(d,T[0], T[1])  for T in CartesianProduct(StandardTableaux(d).list(), StandardTableaux(d).list()) if T[0].shape() == T[1].shape() == lamb]


#return all idempotents, and restrict to a specific diagram if lamb is given 
# NOT NEEDD
def idempotents(d, lamb = -1):
    if lamb == -1:
        return [E(d,T) for T in StandardTableaux(d).list()]
    else:
        return [E(d,T) for T in StandardTableaux(d).list() if T.shape() == lamb]

#matrix unit e_I,J in M_n, or e_IJ in M_n^d if d is given
# NOT NEEDED
def e(I,J,n, d = -1):
    if d == -1:
        return matrix(QQ, n, n, lambda k,l: delta(k+1,l+1,I,J))
    else:
        m_list = [e(I[i], J[i], n) for i in range(0,d)]
        return tensorProduct(m_list)


def EmSd(q,d,n):
    p = q.inverse()
    inds = indices(d,n)
    alltensors = []
    for list in inds:  
        ind = indexSet(list)
        tensors = [e(ind(p[j-1]), ind(j), n) for j in range(1,d+1)]  
        alltensors.append(tensors)         
    merp = [tensorProduct(list) for list in alltensors]
    return sum(merp)

# takes an element x of the group algebra QS(n) and returns a list of (base, coefficient) pairs
def basecoefficientpairs(x, d):
    return [[p, x.coefficient(p)] for p in Permutations(d).list()]    


# returns the inner product of elements in the symmetric group algebra S_d    
def SymInnerProduct(x,y,d):
    return sum([x.coefficient(p) * y.coefficient(p) for p in Permutations(d).list()])     


# embeds QS(d) into M_n \otimes .... \otimes M_n
def algToTensor(p,d,n):
   return sum([x[1]*EmSd(x[0], d, n) for x in basecoefficientpairs(p, d)]) 


# NOT NEEDED
def tabPairs(d, lamb = -1):
    if lamb == -1:
        return [S for S in CartesianProduct(StandardTableaux(d).list(), StandardTableaux(d).list()).list() if S[0].shape() == S[1].shape()] 
    else:
        return [S for S in CartesianProduct(StandardTableaux(d).list(), StandardTableaux(d).list()).list() if S[0].shape() == S[1].shape() == lamb] 



def MatrixBasis(d,n, lamb = -1):
    Return = []
    Pairs = tabPairs(d,lamb)
    for pair in Pairs:
        if len(pair[0]) <= n:
            Return.append(algToTensor(E(d,pair[0], pair[1]), d, n))
    return Return


def Integral_TS(I,J,K,L,d,n, T, S):

    if (len(T) > n):
        return 0

    E_TS = algToTensor(E(d,T,S),d,n)

    return (1/InnerProduct(E_TS,E_TS,n**d)) * InnerProduct(e(J,L,n,d),E_TS,n**d) * InnerProduct(E_TS,e(I,K,n,d), n**d)

# calculate integral
# if we specify lambda, only sum up for that diagram
def Integral(I,J,K,L,d,n, lamb = -1):

    sett = MatrixBasis(d,n, lamb)
    sum = 0

    for A in sett:
        if InnerProduct(A,A,n**d) != 0:
            sum += (1/InnerProduct(A,A,n**d)) * InnerProduct(e(J,L,n,d),A,n**d) * InnerProduct(A,e(I,K,n,d), n**d) 

    return sum


def SUNTest0(T,S):

    I = [1,1,2,2]
    J = [1,2,1,2]
    K = [2,1,2,1]
    L = [1,2,2,1]
   
    I1 = Integral_TS(I,J,J,J,4,2,T,S)
    I2 = Integral_TS(I,J,J,K,4,2,T,S)
    I3 = Integral_TS(I,J,J,L,4,2,T,S)

    return "[" + str(I1) + "," + str(I2) + "," + str(I3) + "] => " + str(I1 + I2 - 2 * I3)

def SunTest1():

    T1 = [[1,2,3,4]]
    
    T2 = [[1,2,3], [4]]
    T3 = [[1,2,4], [3]]
    T4 = [[1,3,4], [2]]

    T5 = [[1,2], [3,4]]
    T6 = [[1,3], [2,4]]

    T7 = [[1,2], [3], [4]]
    T8 = [[1,3], [2], [4]]
    T9 = [[1,4], [2], [3]]

    T10 = [[1], [2], [3], [4]]

    with open("out2.txt", "a") as myfile:
        myfile.write("(T1, T1, " + str(SUNTest0(T1, T1)) + ")\n")
        myfile.write("----------------\n")
        myfile.write("(T2, T2, " + str(SUNTest0(T2, T2)) + ")\n")
        myfile.write("(T3, T3, " + str(SUNTest0(T3, T3)) + ")\n")
        myfile.write("(T4, T4, " + str(SUNTest0(T4, T4)) + ")\n")
        myfile.write("(T2, T3, " + str(SUNTest0(T2, T3)) + ")\n")
        myfile.write("(T3, T2, " + str(SUNTest0(T3, T2)) + ")\n")
        myfile.write("(T2, T4, " + str(SUNTest0(T2, T4)) + ")\n")
        myfile.write("(T4, T2, " + str(SUNTest0(T4, T2)) + ")\n")
        myfile.write("(T3, T4, " + str(SUNTest0(T3, T4)) + ")\n")
        myfile.write("(T4, T3, " + str(SUNTest0(T4, T3)) + ")\n")
        myfile.write("----------------\n")
        myfile.write("(T5, T5, " + str(SUNTest0(T5, T5)) + ")\n")
        myfile.write("(T6, T6, " + str(SUNTest0(T6, T6)) + ")\n")
        myfile.write("(T5, T6, " + str(SUNTest0(T5, T6)) + ")\n")
        myfile.write("(T6, T5, " + str(SUNTest0(T6, T5)) + ")\n")
        myfile.write("----------------\n")
        myfile.write("(T7, T7, " + str(SUNTest0(T7, T7)) + ")\n")
        myfile.write("(T8, T8, " + str(SUNTest0(T8, T8)) + ")\n")
        myfile.write("(T9, T9, " + str(SUNTest0(T9, T9)) + ")\n")
        myfile.write("(T7, T8, " + str(SUNTest0(T7, T8)) + ")\n")
        myfile.write("(T8, T7, " + str(SUNTest0(T8, T7)) + ")\n")
        myfile.write("(T7, T9, " + str(SUNTest0(T7, T9)) + ")\n")
        myfile.write("(T9, T7, " + str(SUNTest0(T9, T7)) + ")\n")
        myfile.write("(T8, T9, " + str(SUNTest0(T8, T9)) + ")\n")
        myfile.write("(T9, T8, " + str(SUNTest0(T9, T8)) + ")\n")
        myfile.write("----------------\n")
        myfile.write("(T10, T10, " + str(SUNTest0(T10, T10)) + ")\n")

def sun_test(T = [[1,3], [2,4]],S = [[1,2], [3,4]]):

    I = [1,2,1,2]
    J = [1,1,2,2]

    I0 = [1,2,1,2] # (1 x 1)
    I1 = [1,2,2,1] # (1 x \sigma)
    I2 = [2,1,1,2] # (\sigma x 1)
    I3 = [2,1,2,1] # (\sigma x \sigma)

    Int0 = Integral_TS(I,J,I,I0,4,2, T,S)
    Int1 = Integral_TS(I,J,I,I1,4,2, T,S)
    Int2 = Integral_TS(I,J,I,I2,4,2, T,S)
    Int3 = Integral_TS(I,J,I,I3,4,2, T,S)

    sum =  Int0 - Int1 - Int2 + Int3
    return ([Int0, Int1, Int2, Int3], sum)


def SunTest4():

    T1 = [[1,2,3,4]]
    
    T2 = [[1,2,3], [4]]
    T3 = [[1,2,4], [3]]
    T4 = [[1,3,4], [2]]

    T5 = [[1,2], [3,4]]
    T6 = [[1,3], [2,4]]

    T7 = [[1,2], [3], [4]]
    T8 = [[1,3], [2], [4]]
    T9 = [[1,4], [2], [3]]

    T10 = [[1], [2], [3], [4]]

    with open("out3.txt", "a") as myfile:
        myfile.write("(T1, T1, " + str(sun_test(T1, T1)) + ")\n")
        myfile.write("----------------\n")
        myfile.write("(T2, T2, " + str(sun_test(T2, T2)) + ")\n")
        myfile.write("(T2, T3, " + str(sun_test(T2, T3)) + ")\n")
        myfile.write("(T2, T4, " + str(sun_test(T2, T4)) + ")\n")
        myfile.write("(T3, T2, " + str(sun_test(T3, T2)) + ")\n")
        myfile.write("(T3, T3, " + str(sun_test(T3, T3)) + ")\n")
        myfile.write("(T3, T4, " + str(sun_test(T3, T4)) + ")\n")
        myfile.write("(T4, T2, " + str(sun_test(T4, T2)) + ")\n")
        myfile.write("(T4, T3, " + str(sun_test(T4, T3)) + ")\n")
        myfile.write("(T4, T4, " + str(sun_test(T4, T4)) + ")\n")
        myfile.write("----------------\n")
        myfile.write("(T5, T5, " + str(sun_test(T5, T5)) + ")\n")
        myfile.write("(T5, T6, " + str(sun_test(T5, T6)) + ")\n")
        myfile.write("(T6, T5, " + str(sun_test(T6, T5)) + ")\n")
        myfile.write("(T6, T6, " + str(sun_test(T6, T6)) + ")\n")
        myfile.write("----------------\n")
        myfile.write("(T7, T7, " + str(sun_test(T7, T7)) + ")\n")
        myfile.write("(T7, T8, " + str(sun_test(T7, T8)) + ")\n")
        myfile.write("(T7, T9, " + str(sun_test(T7, T9)) + ")\n")
        myfile.write("(T8, T7, " + str(sun_test(T8, T7)) + ")\n")
        myfile.write("(T8, T8, " + str(sun_test(T8, T8)) + ")\n")
        myfile.write("(T8, T9, " + str(sun_test(T8, T9)) + ")\n")
        myfile.write("(T9, T7, " + str(sun_test(T9, T7)) + ")\n")
        myfile.write("(T9, T8, " + str(sun_test(T9, T8)) + ")\n")
        myfile.write("(T9, T9, " + str(sun_test(T9, T9)) + ")\n")
        myfile.write("----------------\n")
        myfile.write("(T10, T10, " + str(sun_test(T10, T10)) + ")\n")



#answer is -1/6
def SUNTest():

    I = [1,1,2,2]
    J = [1,2,1,2]
    K = [2,1,2,1]
    L = [1,2,2,1]

    T = [[1,2,3],[4]]
    S = [[1,2,3],[4]]

    I1 = Integral_TS(I,J,J,J,4,2, T,S)
    I2 = Integral_TS(I,J,J,K,4,2, T,S)
    I3 = Integral_TS(I,J,J,L,4,2, T,S)

    print I1
    print I2
    print I3

    return I1 + I2 - 2 * I3


def PermIndices():
    I = [1,2,3,1,2,3,1,2,3]
    inds = []
    for p in CartesianProduct(Permutations(3).list(), Permutations(3).list(), Permutations(3).list()):
        p1 = p[0]
        p2 = p[1]
        p3 = p[2]
        inds.append(((p1.signature()*p2.signature()*p3.signature()), [p1(1), p1(2), p1(3), p2(1), p2(2), p2(3), p3(1), p3(2), p3(3)]))
    return inds


def SunTest2() :

    I = [1,2,3,1,2,3,1,2,3]
    J = [1,1,1,2,2,2,3,3,3]
    inds = PermIndices()

    integral = 0

    for ind in inds:
        s = ind[0]
        K = ind[1]
        integral += (s * Integral_TS(I,J,I,K, 9, 3, [[1,2,3,4,5,6,7,8],[9]], [[1,2,3,4,5,6,7,9],[8]]))
    
    return integral


def mat_coeff(x,d,n,I,J):
# returns coefficient <x, e_IJ>

        e_IJ = e(I,J,n,d)
        inn = InnerProduct(e_IJ, x, d**n) 
        return (I,J,inn)

def mat_coeff_test():
 
    d = 4
    n = 2
    I = [1,2,1,2]
    J = [1,2,1,2]

    S = [[1,2],[3,4]]
    T = [[1,3],[2,4]]

    x = algToTensor(E(d,T,S),d,n)

    return mat_coeff(x,d,n,I,J)

def Coefftest(x, d, n, I = []):
# returns coefficient triples (I,J, <x, e_IJ>)

    inds = CartesianProduct(indices(d,n), indices(d,n))
    ret = []
    for ind in inds:
        if I != []:
            if ind[0] == I:
                inn = mat_coeff(x,d,n,ind[0], ind[1])
                if inn[2] != 0 :
                    ret.append(inn)
        else: 
            inn = mat_coeff(x,d,n,ind[0], ind[1])
            if inn[2] != 0 :
                ret.append(inn)

    return ret

def bar(T,d):
    tBar = []
    for row in T:
        newRow = []
        for entry in row:
            if entry != d:
                newRow.append(entry)
        if newRow != []:
            tBar.append(newRow)
    return StandardTableaux(d-1)(tBar)

def P(T,d):
# returns the product P(T,d) = \prod_{T!=S} ....
    T = StandardTableaux(d)(T)
    prod = 1
    tabs = StandardTableaux(d).list()
    for S in tabs:
        if (T != S) and (bar(T,d) == bar(S,d)):
            X_d = SymmetricGroupAlgebra(QQ,d).jucys_murphy(d)
            aS = S.content(d)
            aT = T.content(d)
            term = (-X_d + aS) / (aS - aT)
            return term
            prod *= term
    return prod

def ets_coeff():
   
    d = 4
    n = 2

    I = [1,2,1,2]
    J = [1,1,2,2]

    T = [[1,2,3,4]]
    S = [[1,3],[2,4]]
    R = [[1,2],[3,4]]
    R = [[1,2,3],[4]]

    p = P(R,d)
    p = SymmetricGroupAlgebra(QQ,d).jucys_murphy(d);
    x = algToTensor(p**1,d,n)

#   interesting: x = yjm(4)
#   if you sum up all the coefficients for e_J,\sigmaI of x you get :

#   for x^0 you get 1
#   for x^1 you get 2
#   for x^2 you get 6
#   for x^3 you get 18
#   for x^4 you get 54
#   for x^5 you get 162
#   for x^k you get .... obvious pattern though

    return Coefftest(x, d, n,J) 

def test():

    T = [[1,3],[2,4]]
    S = [[1,2],[3,4]]
    d = 4
    n = 2
    P = algToTensor(e_hat(T) + e_hat(S), d,n)
    E_TS = algToTensor(E(d,T,S),d,n)

    return Coefftest(E_TS)


def SymmetrizerTest():

    T = [[1,2],[3,4]]
    S = [[1,3],[2,4]]

    E_TS = E(4, T,S)
    E_ST = E(4,S,T)
    E_T = E(4,T,T)
    E_S = E(4,S,S)
    P = e_hat(T) + e_hat(S)
    return  E_T * P - E_T

def ProjectionTest():

    cs = CS(4)
    A = cs.algebra

    I = [1,2,1,2]
    J = [1,1,2,2]

    I0 = [1,2,1,2] # (1 x 1)
    I1 = [1,2,2,1] # (1 x \sigma)
    I2 = [2,1,1,2] # (\sigma x 1)
    I3 = [2,1,2,1] # (\sigma x \sigma)

    P = e_hat([[1,3],[2,4]]) 
    x = algToTensor(P,4,2)

    eII = e(I,I,4,2)
    eJI0 = e(J,I0,4,2) 
    eJI1 = e(J,I1,4,2) 
    eJI2 = e(J,I2,4,2) 
    eJI3 = e(J,I3,4,2) 

    return x
    return InnerProduct(x, eII, 16)

# return sigma(I) for an index set I
def permute(index, perm, d):
    ret = []
    for k in range(1,d+1):
        ret.append(index[perm(k)-1])
    return ret

#return content vector of T
def contents(T,d):
    T = StandardTableaux(d)(T)
    cont = []
    for i in range(1,d+1):
        cont.append(T.content(i))
    return cont

def lamb_content(lamb, d):
    ret = []
    for i in range(0,len(lamb)):
        for j in range(-i, -i + lamb[i]):
            ret.append(j)
    return ret

#returns <E_TS, E_TS>
def SquaredNorm(lamb,d,n):
    dim = StandardTableaux(lamb).cardinality()
    ret = dim / fact(d)
    for c in lamb_content(lamb,d):
        ret *= (n + c)
    return ret

# inner product <x, e_IJ>, for x in C[S_d]
def ip(x,d,n,I,J):
    perms = Permutations(d).list()
    ret = 0
    for sigma in perms:
        if permute(I,sigma, d) == J:
            ret += x.coefficient(sigma)
    return ret

def new_integral(I,J,K,L, d, n):
    base = orthogonalBasis(d)
    sum = 0
    for E_TS in base:
        sum +=  ip(E_TS,d,n,J,L) * ip(E_TS,d,n,I,K)
    return sum
    
# return a dictionary of coefficient => [permutation] for an element of the symmetric group algebra
def MatUnitTest():

    J = [1,1,2,2]
    I = [1,2,1,2]

    T = [[1,2],[3,4]]
    S = [[1,3],[2,4]]

    dd = 4

    x = E(dd,S,T)
    d = {}

    for p in Permutations(dd).list():
        if x.coefficient(p) not in d.keys():
            d[x.coefficient(p)] = []
        d[x.coefficient(p)].append((p, I == permute(J, p, dd)));

    return d

