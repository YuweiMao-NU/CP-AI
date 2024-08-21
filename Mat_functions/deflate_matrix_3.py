# Crystal Plasticity Code.

import math,sys


#################  Reduce 4th order 3*3*3*3 tensor to 2nd 6*6 matrix ########

def reduce_mat(as_mat): 
       as_redu = [[0.0 for i in range(6)] for j in range(6)]

       for i in range(3):
           for j in range(3):
               as_redu[i][j] = as_mat[i][i][j][j]

       for i in range(3):
           as_redu[i][3] = as_mat[i][i][0][1]+as_mat[i][i][1][0]
           as_redu[i][4] = as_mat[i][i][1][2]+as_mat[i][i][2][1]
           as_redu[i][5] = as_mat[i][i][0][2]+as_mat[i][i][2][0]

       for j in range(3):
           as_redu[3][j] = as_mat[0][1][j][j]
           as_redu[4][j] = as_mat[1][2][j][j]
           as_redu[5][j] = as_mat[0][2][j][j]

       as_redu[3][3] = as_mat[0][1][0][1]+as_mat[0][1][1][0]
       as_redu[3][4] = as_mat[0][1][1][2]+as_mat[0][1][2][1]
       as_redu[3][5] = as_mat[0][1][0][2]+as_mat[0][1][2][0]

       as_redu[4][3] = as_mat[1][2][0][1]+as_mat[1][2][1][0]
       as_redu[4][4] = as_mat[1][2][1][2]+as_mat[1][2][2][1]
       as_redu[4][5] = as_mat[1][2][0][2]+as_mat[1][2][2][0]

       as_redu[5][3] = as_mat[0][2][0][1]+as_mat[0][2][1][0]
       as_redu[5][4] = as_mat[0][2][1][2]+as_mat[0][2][2][1]
       as_redu[5][5] = as_mat[0][2][0][2]+as_mat[0][2][2][0]

       return as_redu

    
#################  end of Reduce 4th order 3*3*3*3 tensor to 2nd 6*6 matrix ########

#################  Reduce 4th order 3*3*3*3 W_mat tensor to 2nd 6*6 matrix ########

def reduce_wmat(as_mat): 
       as_redu = [[0.0 for i in range(6)] for j in range(6)]

       for i in range(3):
           for j in range(3):
               as_redu[i][j] = as_mat[i][i][j][j]

       for i in range(3):
           as_redu[i][3] = (as_mat[i][i][0][1]+as_mat[i][i][1][0])/2.0
           as_redu[i][4] = (as_mat[i][i][1][2]+as_mat[i][i][2][1])/2.0
           as_redu[i][5] = (as_mat[i][i][0][2]+as_mat[i][i][2][0])/2.0

       for j in range(3):
           as_redu[3][j] = as_mat[0][1][j][j]
           as_redu[4][j] = as_mat[1][2][j][j]
           as_redu[5][j] = as_mat[0][2][j][j]

       as_redu[3][3] = (as_mat[0][1][0][1]+as_mat[0][1][1][0])/2.0
       as_redu[3][4] = (as_mat[0][1][1][2]+as_mat[0][1][2][1])/2.0
       as_redu[3][5] = (as_mat[0][1][0][2]+as_mat[0][1][2][0])/2.0

       as_redu[4][3] = (as_mat[1][2][0][1]+as_mat[1][2][1][0])/2.0
       as_redu[4][4] = (as_mat[1][2][1][2]+as_mat[1][2][2][1])/2.0
       as_redu[4][5] = (as_mat[1][2][0][2]+as_mat[1][2][2][0])/2.0

       as_redu[5][3] = (as_mat[0][2][0][1]+as_mat[0][2][1][0])/2.0
       as_redu[5][4] = (as_mat[0][2][1][2]+as_mat[0][2][2][1])/2.0
       as_redu[5][5] = (as_mat[0][2][0][2]+as_mat[0][2][2][0])/2.0

       return as_redu

    
#################  end of Reduce 4th order 3*3*3*3 tensor to 2nd 6*6 matrix ########

