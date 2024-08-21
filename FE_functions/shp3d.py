# Crystal Plasticity Code.

import math,sys
from Mat_functions import inverse3



##### Get global shape functions and their derivatives at a given gauss point#####
############################ N_i(X,t) and gradN_i(X,t) ###########################
def SHP3D(ss,tt,zz,XL=[[]]):
    si = [-1,1,1,-1,-1,1,1,-1]
    ti = [-1,-1,1,1,-1,-1,1,1]
    zi = [-1,-1,-1,-1,1,1,1,1]
    SHPL = [[0. for ii in range(8)] for jj in range(4)]
#----------- Local shape functions with derivatives with repect to si,ti and zi
    for i in range(8):
        N1 = 1.0+ss*si[i]
        N2 = 1.0+tt*ti[i]
        N3 = 1.0+zz*zi[i]
        SHPL[3][i] = N1*N2*N3/8.0
        SHPL[0][i] = si[i]*N2*N3/8.0
        SHPL[1][i] = ti[i]*N1*N3/8.0
        SHPL[2][i] = zi[i]*N1*N2/8.0
#--------- Get global values in refrence configuration

    XS = [[0. for ii in range(3)] for jj in range(3)]
    SX = [[0. for ii in range(3)] for jj in range(3)]
 

    for i in range(3):
        for j in range(3):
            for k in range(8):
            	XS[i][j]+=XL[k][j]*SHPL[i][k]
             
    SXD = inverse3.matinv3(XS)
    SX = SXD[0] ; xsJ = SXD[1]
    
    SHPG = [[0. for ii in range(8)] for jj in range(3)]
    for i in range(8):
        T0 = SX[0][0]*SHPL[0][i]+SX[0][1]*SHPL[1][i]+SX[0][2]*SHPL[2][i]
        T1 = SX[1][0]*SHPL[0][i]+SX[1][1]*SHPL[1][i]+SX[1][2]*SHPL[2][i]
        T2 = SX[2][0]*SHPL[0][i]+SX[2][1]*SHPL[1][i]+SX[2][2]*SHPL[2][i]
        SHPG[0][i] = T0
        SHPG[1][i] = T1
        SHPG[2][i] = T2
                
    return SHPG,xsJ
############### End of calculation of shape functions and drivatives ##################

