# Crystal Plasticity Code.

import math,sys


########### Calculate deformation gradient (dx/dX=sum(x*gradN_i(X,t))) ################
def calc_F(cord,SHP):
    
    dxt_dx0 = 0.0 ; dxt_dy0 = 0.0 ; dxt_dz0 = 0.0
    dyt_dx0 = 0.0 ; dyt_dy0 = 0.0 ; dyt_dz0 = 0.0
    dzt_dx0 = 0.0 ; dzt_dy0 = 0.0 ; dzt_dz0 = 0.0

    for I in range(8):
        dxt_dx0 += SHP[0][I]*cord[I][0]
        dxt_dy0 += SHP[1][I]*cord[I][0]
        dxt_dz0 += SHP[2][I]*cord[I][0]

        dyt_dx0 += SHP[0][I]*cord[I][1]
        dyt_dy0 += SHP[1][I]*cord[I][1]
        dyt_dz0 += SHP[2][I]*cord[I][1]

        dzt_dx0 += SHP[0][I]*cord[I][2]
        dzt_dy0 += SHP[1][I]*cord[I][2]
        dzt_dz0 += SHP[2][I]*cord[I][2]

        FG = [[0. for ii in range(3)] for jj in range(3)]
        FG[0][0] = dxt_dx0 ; FG[0][1] = dxt_dy0 ; FG[0][2] = dxt_dz0
        FG[1][0] = dyt_dx0 ; FG[1][1] = dyt_dy0 ; FG[1][2] = dyt_dz0
        FG[2][0] = dzt_dx0 ; FG[2][1] = dzt_dy0 ; FG[2][2] = dzt_dz0


    return FG
#################### End of calculation of deformation gradient ###################

