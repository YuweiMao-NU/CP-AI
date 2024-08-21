# Crystal Plasticity Code.

import math,sys
from FE_functions import gauss
from FE_functions import shp3d






################### calculation of stiffness matrix ################
def STIF3D(ngauss,nnode,xyz,strs,dep):
    vol_gp = [0.0 for i in range(ngauss)]
    gaussw = gauss.pgauss()
             
    ES = [[0.0 for i in range(24)] for j in range(24)]
    ESGEOM = [[0.0 for i in range(24)] for j in range(24)]
    for ig in range(ngauss):

        SHP = shp3d.SHP3D(gaussw[0][ig],gaussw[1][ig],gaussw[2][ig],xyz)
        xsj = SHP[1]*gaussw[3][ig]

        vol_gp[ig] = xsj

        B = [[0.0 for i in range(24)] for j in range(6)]
        for i in range(nnode):
            ii = i*3
            B[0][ii+0] = SHP[0][0][i] ; B[0][ii+1] = 0.0          ; B[0][ii+2] = 0.0
            B[1][ii+0] = 0.0          ; B[1][ii+1] = SHP[0][1][i] ; B[1][ii+2] = 0.0
            B[2][ii+0] = 0.0          ; B[2][ii+1] = 0.0          ; B[2][ii+2] = SHP[0][2][i]
            B[3][ii+0] = SHP[0][1][i] ; B[3][ii+1] = SHP[0][0][i] ; B[3][ii+2] = 0.0
            B[4][ii+0] = 0.0          ; B[4][ii+1] = SHP[0][2][i] ; B[4][ii+2] = SHP[0][1][i]
            B[5][ii+0] = SHP[0][2][i] ; B[5][ii+1] = 0.0          ; B[5][ii+2] = SHP[0][0][i]

        DB = [[0.0 for i in range(24)] for j in range(6)]
        for i in range(6):
            for j in range(24):
                DB[i][j] = 0.0
                for k in range(6):
                    DB[i][j] += dep[ig][i][k]*B[k][j]

        for i in range(24):
            for j in range(24):
                for k in range(6):
                    ES[i][j] += B[k][i]*DB[k][j]*xsj




        BNL = [[0.0 for i in range(24)] for j in range(9)]
        DBNL = [[0.0 for i in range(24)] for j in range(9)]
        stressgeo = [[0.0 for i in range(9)] for j in range(9)]

        for i in range(8):
            icol1 = i*3+0
            icol2 = icol1+1
            icol3 = icol2+1

            BNL[0][icol1] = SHP[0][0][i]
            BNL[1][icol1] = SHP[0][1][i]
            BNL[2][icol1] = SHP[0][2][i]
            BNL[3][icol2] = SHP[0][0][i]
            BNL[4][icol2] = SHP[0][1][i]
            BNL[5][icol2] = SHP[0][2][i]
            BNL[6][icol3] = SHP[0][0][i]
            BNL[7][icol3] = SHP[0][1][i]
            BNL[8][icol3] = SHP[0][2][i]

        for i in range(0,9,3):
            stressgeo[i][i] = strs[ig][0]
            stressgeo[i][i+1] = strs[ig][3]
            stressgeo[i][i+2] = strs[ig][5]
            stressgeo[i+1][i] = strs[ig][3]
            stressgeo[i+1][i+1] = strs[ig][1]
            stressgeo[i+1][i+2] = strs[ig][4]
            stressgeo[i+2][i] = strs[ig][5]
            stressgeo[i+2][i+1] = strs[ig][4]
            stressgeo[i+2][i+2] = strs[ig][2]


        for i in range(9):
            for j in range(24):
                for k in range(9):
                    DBNL[i][j] += stressgeo[i][k]*BNL[k][j]

        for i in range(24):
            for j in range(24):
                for k in range(9):
                    ESGEOM[i][j] += BNL[k][i]*DBNL[k][j]*xsj


        for i in range(24):
            for j in range(24):
                    ES[i][j] += ESGEOM[i][j] 


    return ES,vol_gp

################# end of element stiffness matrix calculation #########

      
