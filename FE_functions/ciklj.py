# Crystal Plasticity Code.

import math,sys






############## Elastic Material 2nd and 4th order tensor ###################

def cmat(c11,c12,c44,qrot):

    cij = [[0. for i in range(6)] for j in range(6)] 
    
#---------------- cij for cubic materials 
    cij[0][0] = cij[1][1] = cij[2][2] = c11
    cij[0][1] = cij[0][2] = cij[1][0] = cij[1][2] = cij[2][0] = cij[2][1] = c12
    cij[3][3] = cij[4][4] = cij[5][5] = c44


    cijkl = [[[[0. for i in range(3)] for j in range(3)] for k in range(3)] for l in range(3)] 

#----------------- Get cijkl  
    for i in range(3):
        for j in range(3):
            cijkl[i][i][j][j] = cij[i][j]

    for i in range(3):
        cijkl[i][i][0][1] = cij[i][3]
        cijkl[i][i][1][0] = cij[i][3]
        cijkl[i][i][1][2] = cij[i][4]
        cijkl[i][i][2][1] = cij[i][4]
        cijkl[i][i][0][2] = cij[i][5]
        cijkl[i][i][2][0] = cij[i][5]
        
        cijkl[0][1][i][i] = cijkl[1][0][i][i] = cij[3][i]
        cijkl[1][2][i][i] = cijkl[2][1][i][i] = cij[4][i]
        cijkl[0][2][i][i] = cijkl[2][0][i][i] = cij[5][i]
    
    cijkl[0][1][0][1] = cijkl[0][1][1][0] = cijkl[1][0][0][1] = cijkl[1][0][1][0] = cij[3][3]
    cijkl[0][1][1][2] = cijkl[0][1][2][1] = cijkl[1][0][1][2] = cijkl[1][0][2][1] = cij[3][4]
    cijkl[0][1][0][2] = cijkl[0][1][2][0] = cijkl[1][0][0][2] = cijkl[1][0][2][0] = cij[3][5]
    cijkl[1][2][0][1] = cijkl[1][2][1][0] = cijkl[2][1][0][1] = cijkl[2][1][1][0] = cij[4][3]
    cijkl[1][2][1][2] = cijkl[1][2][2][1] = cijkl[2][1][1][2] = cijkl[2][1][2][1] = cij[4][4]
    cijkl[1][2][0][2] = cijkl[1][2][2][0] = cijkl[2][1][0][2] = cijkl[2][1][2][0] = cij[4][5]
    cijkl[0][2][0][1] = cijkl[0][2][1][0] = cijkl[2][0][0][1] = cijkl[2][0][1][0] = cij[5][3]
    cijkl[0][2][1][2] = cijkl[0][2][2][1] = cijkl[2][0][1][2] = cijkl[2][0][2][1] = cij[5][4]
    cijkl[0][2][0][2] = cijkl[0][2][2][0] = cijkl[2][0][0][2] = cijkl[2][0][2][0] = cij[5][5]


#----------------- C_MAT for Crytal with respect to orientation 
    C_mat = [[[[0. for i in range(3)] for j in range(3)] for k in range(3)] for l in range(3)] 
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    
                    for m in range(3):
                        for n in range(3):
                            for p in range(3):
                                for q in range(3):
                                    C_mat[i][j][k][l] += qrot[i][m]*qrot[j][n]*qrot[k][p]*qrot[l][q]*cijkl[m][n][p][q]


    return C_mat
########################## End of calculation C_mat ##########################

