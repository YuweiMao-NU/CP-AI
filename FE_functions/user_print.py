# Crystal Plasticity Code.

import math,sys


      


########### post processing ################
def user_print(tnode,nelem,ngauss,gdisp,garray_tau,vol_g):
    
    stressvol = [0.0 for i in range(6)]
    strainvol = [0.0 for i in range(6)]

    voltot = 0.0

    for ielem in range(nelem):
        for ig in range(ngauss):
            TMP1 = garray_tau[ielem][ig].defg
            TMPt = list(zip(*TMP1))
            TMP2 = [[0.0 for i in range(3)] for j in range(3)]
            TMP2[0][0] = TMP2[1][1] = TMP2[2][2] = 1.0
            TMP3 = [[0.0 for i in range(3)] for j in range(3)]
            for i in range(3):
                for j in range(3):
                    for k in range(3):
                        TMP3[i][j] += TMPt[i][k]*TMP1[k][j]

            for i in range(3):
                for j in range(3):
                    TMP3[i][j] = 0.5*(TMP3[i][j]-TMP2[i][j])

            strain = [0.0 for i in range(6)]
            strain[0] = TMP3[0][0] ; strain[1] = TMP3[1][1] ; strain[2] = TMP3[2][2] 
            strain[3] = TMP3[0][1] ; strain[4] = TMP3[1][2] ; strain[5] = TMP3[0][2] 

            strainvol[0] += strain[0]*vol_g[ielem][ig]
            strainvol[1] += strain[1]*vol_g[ielem][ig]
            strainvol[2] += strain[2]*vol_g[ielem][ig]

    for ielem in range(nelem):
        for ig in range(ngauss):
            stressvol[0] += garray_tau[ielem][ig].cauchy[0][0]*vol_g[ielem][ig]
            stressvol[1] += garray_tau[ielem][ig].cauchy[1][1]*vol_g[ielem][ig]
            stressvol[2] += garray_tau[ielem][ig].cauchy[2][2]*vol_g[ielem][ig]
            voltot += vol_g[ielem][ig]

    stressvol[0] = stressvol[0]/voltot
    stressvol[1] = stressvol[1]/voltot
    stressvol[2] = stressvol[2]/voltot

    strainvol[0] = strainvol[0]/voltot
    strainvol[1] = strainvol[1]/voltot
    strainvol[2] = strainvol[2]/voltot


    return stressvol,strainvol
    
