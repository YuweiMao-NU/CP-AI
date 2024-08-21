# Crystal Plasticity Code.

import math,sys
from Mat_functions import matrix_mlutip
from Mat_functions import inverse3




########################## Calculate matrix [A] ##############################
def calc_A(F,Fp):
    F_tau = [[0. for ii in range(3)] for jj in range(3)]
    FT_tau = [[0. for ii in range(3)] for jj in range(3)]
    FpI_t = [[0. for ii in range(3)] for jj in range(3)]
    FpIT_t = [[0. for ii in range(3)] for jj in range(3)]
    TMP1 = [[0. for ii in range(3)] for jj in range(3)]
    TMP2 = [[0. for ii in range(3)] for jj in range(3)]
    A_mat = [[0. for ii in range(3)] for jj in range(3)]

    FpID = inverse3.matinv3(Fp)
    FpI_t = FpID[0] ; FpIT_t = list(zip(*FpI_t))
    
    for i in range(3):
        for j in range(3):
            F_tau[i][j] = F[i][j]
            FT_tau[i][j] = F[j][i]

    TMP1 = matrix_mlutip.matmul(F_tau,FpI_t)
    TMP2 = matrix_mlutip.matmul(FT_tau,TMP1)
    A_mat = matrix_mlutip.matmul(FpIT_t,TMP2)

    return A_mat
########################## End of A_mat ##################################


################ Calculate matrix [B_alpha] ##############################
def calc_B(n_slip,A_mat,schmid):
    B_alpha = [[[0. for i in range(3)] for j in range(3)] for k in range(n_slip)]
    TMP1 = [[0. for ii in range(3)] for jj in range(3)]
    TMP1T = [[0. for ii in range(3)] for jj in range(3)]
    TMP20 = [[0. for ii in range(3)] for jj in range(3)]
    TMP21 = [[0. for ii in range(3)] for jj in range(3)]
    TMP2 = [[0. for ii in range(3)] for jj in range(3)]
    for k in range(n_slip):
        for j in range(3):
            for i in range(3):
                TMP1[i][j] = schmid[k][i][j]
                TMP1T[i][j] = schmid[k][j][i]
                

        TMP20 = matrix_mlutip.matmul(A_mat,TMP1)
        TMP21= matrix_mlutip.matmul(TMP1T,A_mat)

        for it in range(3):
            for jt in range(3):
                TMP2[it][jt] = TMP20[it][jt]+TMP21[it][jt]
        
        for i in range(3):
            for j in range(3):
                B_alpha[k][i][j]=TMP2[i][j]
        
    return B_alpha
###################### End of B_alpha #####################################
    
############# Calculate the trial elastic stress S_trial ##################
def stress_trial(C_mat,A_mat):
    delta_kron = [[0. for ii in range(3)] for jj in range(3)]
    S_trial = [[0. for ii in range(3)] for jj in range(3)]
    
    delta_kron[0][0] = delta_kron[1][1]= delta_kron[2][2] = 1.0
	
    for i in range(3):
        for j in range(3):
            for m in range(3):
                for n in range(3):
                    S_trial[i][j] += C_mat[i][j][m][n]*0.5*(A_mat[m][n]-delta_kron[m][n])

    return S_trial
########################## End of S_trial ####################################

################ Calculate the matrix C_ALPHA for each slip system ###########
def calc_C_alpha(n_slip,C_mat,B_alpha):
    C_alpha = [[[0. for i in range(3)] for j in range(3)] for k in range(n_slip)]    

    for k in range(n_slip):
        for i in range(3):
            for j in range(3):
                for m in range(3):
                    for n in range(3):
                        C_alpha[k][i][j] += C_mat[i][j][m][n]*0.5*B_alpha[k][m][n]

    return C_alpha


####################### End of C_alpha ####################################

