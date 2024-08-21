# Crystal Plasticity Code.

import math,sys


######### Inflate 6*6 matrix to 4rt order tenso ################
def inflate_ten(as_ten):
    as_mat = [[[[0.0 for i in range(3)] for j in range(3)] for k in range(3)] for l in range(3)]
    for i in range(3):
        for j in range(3):
            as_mat[i][i][j][j] = as_ten[i][j]

    for i in range(3):
        as_mat[i][i][0][1] = as_mat[i][i][1][0] = 0.5*as_ten[i][3]
        as_mat[i][i][1][2] = as_mat[i][i][2][1] = 0.5*as_ten[i][4]
        as_mat[i][i][0][2] = as_mat[i][i][2][0] = 0.5*as_ten[i][5]

    for j in range(3):
        as_mat[0][1][j][j] = as_mat[1][0][j][j] = as_ten[3][j]
        as_mat[1][2][j][j] = as_mat[2][1][j][j] = as_ten[4][j]
        as_mat[0][2][j][j] = as_mat[2][0][j][j] = as_ten[5][j]

    as_mat[0][1][0][1] = as_mat[0][1][1][0] = as_mat[1][0][0][1] = as_mat[1][0][1][0] = 0.5*as_ten[3][3]
    as_mat[0][1][1][2] = as_mat[0][1][2][1] = as_mat[1][0][1][2] = as_mat[1][0][2][1] = 0.5*as_ten[3][4]
    as_mat[0][1][0][2] = as_mat[0][1][2][0] = as_mat[1][0][0][2] = as_mat[1][0][2][0] = 0.5*as_ten[3][5]

    as_mat[1][2][0][1] = as_mat[1][2][1][0] = as_mat[2][1][0][1] = as_mat[2][1][1][0] = 0.5*as_ten[4][3]
    as_mat[1][2][1][2] = as_mat[1][2][2][1] = as_mat[2][1][1][2] = as_mat[2][1][2][1] = 0.5*as_ten[4][4]
    as_mat[1][2][0][2] = as_mat[1][2][2][0] = as_mat[2][1][0][2] = as_mat[2][1][2][0] = 0.5*as_ten[4][5]

    as_mat[0][2][0][1] = as_mat[0][2][1][0] = as_mat[2][0][0][1] = as_mat[2][0][1][0] = 0.5*as_ten[5][3]
    as_mat[0][2][1][2] = as_mat[0][2][2][1] = as_mat[2][0][1][2] = as_mat[2][0][2][1] = 0.5*as_ten[5][4]
    as_mat[0][2][0][2] = as_mat[0][2][2][0] = as_mat[2][0][0][2] = as_mat[2][0][2][0] = 0.5*as_ten[5][5]

    return as_mat
############ end of Inflate of 6*6 matrix to 4rt order tensor ###########

