# Crystal Plasticity Code.

import math,sys
from Mat_functions import deflate_matrix_2




################# Calculation of GN vector ###################

def calc_GNvec(n_slip,S_star,S_trial,dgam,C_alpha):
    GN_vec = [0.0 for i in range(6)]
    TMP1 = [[0.0 for i in range(3)] for j in range(3)]

    for i  in range(3):
        for j in range(3):
                TMP1[i][j] = 0.0
                for k in range(n_slip):
                    TMP1[i][j] += dgam[k]*C_alpha[k][i][j]

    for i in range(3):
        for j in range(3):
            TMP1[i][j] += S_star[i][j]-S_trial[i][j]

    GN_vec = deflate_matrix_2.reduce_vec(TMP1)

    return GN_vec
    

################# end of Calculation of GN vector ###################

