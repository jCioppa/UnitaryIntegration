def fact(n):
    if n == 1:
        return 1
    else:
        return n * fact(n-1)

def delta(a,b,c,d):
    if (a == c) and (b == d):
        return 1
    else:
        return 0

def trace(A,n):
    sum = 0
    for i in range(0,n):
        sum += A[i][i]*A[i][i]
    return sum

def normalized_trace(A,n):
    tr = trace(A,n)
    return tr / n

#inner product of matrices A,B
def InnerProduct(A, B, n):
    sum = 0
    for i in range(0,n):
        for j in range(0,n):
            sum+=A[i][j]*B[i][j]
    return sum

#tensor product of a list of matrices
def tensorProduct(matrixList): 
    d = len(matrixList)   
    if d == 1:
        return matrixList[0]
    if d == 2:
        return matrixList[0].tensor_product(matrixList[1])
    else:    
        return matrixList[0].tensor_product(tensorProduct([matrixList[i] for i in range(1, d)]))


#given a list [i_1, i_2, ..., i_d] returns the function sigma(k) = i_k
def indexSet(list):
    return lambda k: list[k-1]


#return all elements of [1..n]^d as a list of lists
def indices(d,n):
    if d == 0:
        return []
    if d == 1:
        return [[i] for i in range(1,n+1)]
    else:
        ret = []
        for i in range(1,n+1):
            for x in indices(d-1,n):
                ret.append(x + [i])
        return ret
 
