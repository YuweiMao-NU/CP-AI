# Crystal Plasticity Code.

import math,sys



############ Calculate A(n,n)*B(n,n)############

def matmul(A,B):
    C = [[0. for ii in range(len(A))] for jj in range(len(A))]
    for i in range(len(A)):
        for j in range(len(A)):
            for k in range(len(A)):
                C[i][j] += A[i][k]*B[k][j]
    return C

