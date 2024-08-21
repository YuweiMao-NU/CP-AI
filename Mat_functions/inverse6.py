# Crystal Plasticity Code.

import math,sys


      
################# inverse nbyn matrix ######################
def matinv(a,n):
    ipivot = [0 for j in range(n)]
    pivot = [0 for j in range(n)]
    index = [[0 for i in range(2)] for j in range(n)]
    for i in range(n):
        amax=0.0
        for j in range(n):
            if ipivot[j] != 1:
                for k in range(n):
                    if ipivot[k] <= 1:
                        if ipivot[k] != 1:
                            if abs(amax) <= abs(a[j][k]):
                                irow = j
                                icolum = k
                                amax = a[j][k]
                            
        ipivot[icolum] += 1
        if irow != icolum:
            for l in range(n):
                swap = a[irow][l]
                a[irow][l] = a[icolum][l]
                a[icolum][l] = swap
                
        index[i][0] = irow
        index[i][1] = icolum
        pivot[i] = a[icolum][icolum]
        a[icolum][icolum] = 1.0
        
        for l in range(n):
            a[icolum][l] = a[icolum][l]/pivot[i]
        for l1 in range(n):
            if l1 != icolum:
                t = a[l1][icolum]
                a[l1][icolum] = 0.0
                for l in range(n):
                    a[l1][l] = a[l1][l]-a[icolum][l]*t

    for i in range(n):
        l = n-i-1
        if index[l][0] != index[l][1]:
            jrow = index[l][0]
            jcolum = index[l][1]
            for k in range(n):
                swap = a[k][jrow]
                a[k][jrow] = a[k][jcolum]
                a[k][jcolum] = swap

    return a


################# end of inverse nbyn matrix ######################

