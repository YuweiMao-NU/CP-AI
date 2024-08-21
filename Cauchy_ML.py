from Mat_functions import inverse3
from Mat_functions import matrix_mlutip

def Cauchy_ML(F_tau, Fp_tau, C_mat):
    Fe_tau = [[0.0 for i in range(3)] for j in range(3)]
    Fp_tauI = [[0.0 for i in range(3)] for j in range(3)]
    FptauD = inverse3.matinv3(Fp_tau)
    Fp_tauI = FptauD[0]

    for i in range(3):
        for j in range(3):
            for k in range(3):
                Fe_tau[i][j] += F_tau[i][k] * Fp_tauI[k][j]

    E_e = [[0.0 for i in range(3)] for j in range(3)]
    delta_kron = [[0. for ii in range(3)] for jj in range(3)]

    delta_kron[0][0] = delta_kron[1][1] = delta_kron[2][2] = 1.0

    for i in range(3):
        for j in range(3):
            E_e[i][j] = 0.0
            for k in range(3):
                E_e[i][j] += Fe_tau[k][i] * Fe_tau[k][j]
            E_e[i][j] = 0.5 * (E_e[i][j] - delta_kron[i][j])

    S_Cauchy = [[0.0 for i in range(3)] for j in range(3)]
    S_Second = [[0.0 for i in range(3)] for j in range(3)]

    for i in range(3):
        for j in range(3):
            S_Second[i][j] = 0.0
            for m in range(3):
                for n in range(3):
                    S_Second[i][j] += C_mat[i][j][m][n] * E_e[m][n]

    det_F_tau = (F_tau[0][0] * F_tau[1][1] * F_tau[2][2]
                 - F_tau[0][0] * F_tau[1][2] * F_tau[2][1] - F_tau[1][0] * F_tau[0][1] * F_tau[2][2]
                 + F_tau[1][0] * F_tau[0][2] * F_tau[2][1] + F_tau[2][0] * F_tau[0][1] * F_tau[1][2]
                 - F_tau[2][0] * F_tau[0][2] * F_tau[1][1])

    for i in range(3):
        for j in range(3):
            S_Second[i][j] = S_Second[i][j] / det_F_tau

    TPMcs = [[0.0 for i in range(3)] for j in range(3)]
    TPMcs = matrix_mlutip.matmul(S_Second, list(zip(*Fe_tau)))

    S_Cauchy = matrix_mlutip.matmul(Fe_tau, TPMcs)

    for i in range(3):
        for j in range(3):
            S_Cauchy[i][j] = S_Cauchy[i][j] / det_F_tau

    return S_Cauchy, S_Second