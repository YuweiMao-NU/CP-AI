import math,sys

from Mat_functions import deflate_matrix_W
from Mat_functions import matrix_mlutip
from Mat_functions import inverse3
from Mat_functions import inflate_matrix_6

# ---------------------- Polar Decomposition of Rlative Gradient Deformation -------------------

def polardecomp(F_t, F_tau):
    # -------------- Calculation of relative deformation gradient F_REL = F_tau*(F_t)^-1 -----------
    F_REL = [[0.0 for i in range(3)] for j in range(3)]
    F_tI = [[0. for i in range(3)] for j in range(3)]

    F_tID = inverse3.matinv3(F_t)
    F_tI = F_tID[0]

    for i in range(3):
        for j in range(3):
            F_REL[i][j] = 0.0
            for k in range(3):
                F_REL[i][j] += F_tau[i][k] * F_tI[k][j]

    # ----------------------- square root of a positive matrix U=sqrt(C) ------------

    R = [[0. for i in range(3)] for j in range(3)]
    C = [[0. for i in range(3)] for j in range(3)]
    Csquare = [[0. for i in range(3)] for j in range(3)]
    Iden = [[0. for i in range(3)] for j in range(3)]
    U = [[0. for i in range(3)] for j in range(3)]
    invU = [[0. for i in range(3)] for j in range(3)]

    Iden[0][0] = Iden[1][1] = Iden[2][2] = 1.0

    F_REL_T = list(zip(*F_REL))
    # ----------------------- C = FTF --------------------------------------
    for i in range(3):
        for j in range(3):
            for k in range(3):
                C[i][j] += F_REL_T[i][k] * F_REL[k][j]

    # ---------------------- C^2 ------------------------------------------
    for i in range(3):
        for j in range(3):
            for k in range(3):
                Csquare[i][j] += C[i][k] * C[k][j]

    I_C = C[0][0] + C[1][1] + C[2][2]

    # --------------------- Invarients of C -------------------------------
    I_Csquare = Csquare[0][0] + Csquare[1][1] + Csquare[2][2]
    II_C = 0.5 * (I_C ** 2 - I_Csquare)
    III_C = C[0][0] * (C[1][1] * C[2][2] - C[1][2] * C[2][1]) - \
            C[0][1] * (C[1][0] * C[2][2] - C[2][0] * C[1][2]) + \
            C[0][2] * (C[1][0] * C[2][1] - C[1][1] * C[2][0])

    k = I_C ** 2 - 3.0 * II_C

    l = I_C * (I_C ** 2 - 4.5 * II_C) + 13.5 * III_C

    ######## added for initialization ##########

    if abs(k) < 1.0e-10:
        k = 1.0e-10

    ###########################################

    t1 = l / (k ** 1.5)

    if t1 < -1.0:
        t1 = -1.0
    elif t1 > 1.0:
        t1 = 1.0

    teta = (1.0 / 3.0) * math.acos(t1)

    alpha = 2.0 * math.pi / 3.0

    # --------------------- Eigenvalues of U -----------------------------
    lamda1 = math.sqrt((1.0 / 3.0) * (I_C + 2.0 * math.sqrt(k) * math.cos(teta)))
    lamda2 = math.sqrt((1.0 / 3.0) * (I_C + 2.0 * math.sqrt(k) * math.cos(alpha + teta)))
    lamda3 = math.sqrt((1.0 / 3.0) * (I_C + 2.0 * math.sqrt(k) * math.cos(2.0 * alpha + teta)))

    # --------------------- Invarients of U ---------------------------
    I_U = lamda1 + lamda2 + lamda3
    II_U = lamda1 * lamda2 + lamda1 * lamda3 + lamda2 * lamda3
    III_U = lamda1 * lamda2 * lamda3

    # --------------------- U and Inverse of U ------------------------

    for i in range(3):
        for j in range(3):
            U[i][j] = (1.0 / (I_U * II_U - III_U)) * (I_U * III_U * Iden[i][j] + \
                                                      (I_U ** 2 - II_U) * C[i][j] - Csquare[i][j])
            invU[i][j] = (1.0 / III_U) * (II_U * Iden[i][j] - \
                                          I_U * U[i][j] + C[i][j])
    # ----------------------- R = FU^-1 -----------------------------
    for i in range(3):
        for j in range(3):
            for k in range(3):
                R[i][j] += F_REL[i][k] * invU[k][j]

    # ---------------------- End of Polar Decomposition ---------------

    return U, R

def reduce_mat(as_mat):
    as_redu = [[0.0 for i in range(6)] for j in range(6)]

    for i in range(3):
        for j in range(3):
            as_redu[i][j] = as_mat[i][i][j][j]

    for i in range(3):
        as_redu[i][3] = as_mat[i][i][0][1] + as_mat[i][i][1][0]
        as_redu[i][4] = as_mat[i][i][1][2] + as_mat[i][i][2][1]
        as_redu[i][5] = as_mat[i][i][0][2] + as_mat[i][i][2][0]

    for j in range(3):
        as_redu[3][j] = as_mat[0][1][j][j]
        as_redu[4][j] = as_mat[1][2][j][j]
        as_redu[5][j] = as_mat[0][2][j][j]

    as_redu[3][3] = as_mat[0][1][0][1] + as_mat[0][1][1][0]
    as_redu[3][4] = as_mat[0][1][1][2] + as_mat[0][1][2][1]
    as_redu[3][5] = as_mat[0][1][0][2] + as_mat[0][1][2][0]

    as_redu[4][3] = as_mat[1][2][0][1] + as_mat[1][2][1][0]
    as_redu[4][4] = as_mat[1][2][1][2] + as_mat[1][2][2][1]
    as_redu[4][5] = as_mat[1][2][0][2] + as_mat[1][2][2][0]

    as_redu[5][3] = as_mat[0][2][0][1] + as_mat[0][2][1][0]
    as_redu[5][4] = as_mat[0][2][1][2] + as_mat[0][2][2][1]
    as_redu[5][5] = as_mat[0][2][0][2] + as_mat[0][2][2][0]

    return as_redu

# ------------ Forth order Kronecker Delta ---------------------
def delta_kron4():
    delta_kron = [[0. for ii in range(3)] for jj in range(3)]
    delta_kron[0][0] = delta_kron[1][1] = delta_kron[2][2] = 1.0

    delta_kron4d = [[[[0. for i in range(3)] for j in range(3)] for k in range(3)] for l in range(3)]

    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    delta_kron4d[i][j][k][l] += delta_kron[i][k] * delta_kron[j][l] + delta_kron[i][l] * delta_kron[j][
                        k]

    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    delta_kron4d[i][j][k][l] = 0.5 * delta_kron4d[i][j][k][l]

    return delta_kron4d


# ---------------- Calculation of elasto plastic modulus or elastic Jacobian ---------

def W_mat(F_t,F_tau,Fp_tau,Fp,C_mat,schmid,S_star,dgam,dgam_dta, C_alpha):
    n_slip = len(schmid)

    # ---------Elastic stretch and rotation------------
    Poldecomp = polardecomp(F_t, F_tau)
    Strh_el = Poldecomp[0]
    Rot_el = Poldecomp[1]
    # ------------------------------------------------
    # Fe_tau = self._Fe_tau_cache
    Fe_tau = [[0.0 for i in range(3)] for j in range(3)]
    Fp_tauI = [[0.0 for i in range(3)] for j in range(3)]
    FptauD = inverse3.matinv3(Fp_tau)
    Fp_tauI = FptauD[0]

    for i in range(3):
        for j in range(3):
            for k in range(3):
                Fe_tau[i][j] += F_tau[i][k] * Fp_tauI[k][j]
    #-----------------------------------------------------

    Fp_t_inv = [[0.0 for i in range(3)] for j in range(3)]
    Fp_tID = inverse3.matinv3(Fp)
    Fp_t_inv = Fp_tID[0]
    Fe_t = [[0.0 for i in range(3)] for j in range(3)]
    for i in range(3):
        for j in range(3):
            Fe_t[i][j] = 0.0
            for k in range(3):
                Fe_t[i][j] += F_t[i][k] * Fp_t_inv[k][j]

    # ------ Calculation of L_mat-------
    L_mat = [[[[0. for i in range(3)] for j in range(3)] for k in range(3)] for l in range(3)]

    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    L_mat[i][j][k][l] = 0.0
                    for m in range(3):
                        L_mat[i][j][k][l] += Fe_t[k][i] * Strh_el[l][m] * Fe_t[m][j] + Fe_t[m][i] * Strh_el[m][k] * \
                                             Fe_t[l][j]

    # --------- Calculation of D_mat ---------
    D_mat = [[[[0. for i in range(3)] for j in range(3)] for k in range(3)] for l in range(3)]

    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    D_mat[i][j][k][l] = 0.0
                    for m in range(3):
                        for n in range(3):
                            D_mat[i][j][k][l] += 0.5 * C_mat[i][j][m][n] * L_mat[m][n][k][l]

    # -------------- Calculation of G_alpha----------------------
    G_alpha = [[[[[0. for i in range(3)] for j in range(3)] for k in range(3)] for l in range(3)] for alpha in
               range(len(schmid))]

    for alpha in range(len(schmid)):
        for m in range(3):
            for n in range(3):
                for k in range(3):
                    for l in range(3):
                        G_alpha[alpha][m][n][k][l] = 0.0
                        for p in range(3):
                            G_alpha[alpha][m][n][k][l] += L_mat[m][p][k][l] * schmid[alpha][p][n] + \
                                                          schmid[alpha][m][p] * L_mat[p][n][k][l]
    # --------- Calculation of J_alpha ---------
    J_alpha = [[[[[0. for i in range(3)] for j in range(3)] for k in range(3)] for l in range(3)] for alpha in
               range(len(schmid))]

    for alpha in range(len(schmid)):
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    for l in range(3):
                        J_alpha[alpha][i][j][k][l] = 0.0
                        for m in range(3):
                            for n in range(3):
                                J_alpha[alpha][i][j][k][l] += 0.5 * C_mat[i][j][m][n] * G_alpha[alpha][m][n][k][l]

    # ------------- Calculation of Q_mat -------------------
    #--------replace-----------
    # RJ_reduced = self._RJ_reduced_cache
    GT_mat = [[[0. for i in range(3)] for j in range(3)] for k in range(n_slip)]

    for k in range(n_slip):
        for i in range(3):
            for j in range(3):
                GT_mat[k][i][j] = 0.5 * (schmid[k][i][j] + schmid[k][j][i]) * dgam_dta[k]

    RJ_mat = [[[[0. for i in range(3)] for j in range(3)] for k in range(3)] for l in range(3)]
    for i in range(3):
        for j in range(3):
            for m in range(3):
                for n in range(3):
                    for k in range(n_slip):
                        RJ_mat[i][j][m][n] += C_alpha[k][i][j] * GT_mat[k][m][n]

    delta4 = delta_kron4()
    for i in range(3):
        for j in range(3):
            for m in range(3):
                for n in range(3):
                    RJ_mat[i][j][m][n] += delta4[i][j][m][n]

    RJ_reduced = reduce_mat(RJ_mat)


    # ----------------------------

    K_inv = inflate_matrix_6.inflate_ten(RJ_reduced)

    # dgam = self._dgam_cache

    TMP_4d = [[[[0. for i in range(3)] for j in range(3)] for k in range(3)] for l in range(3)]
    for m in range(3):
        for n in range(3):
            for k in range(3):
                for l in range(3):
                    TMP_4d[m][n][k][l] = 0.0
                    for alpha in range(len(schmid)):
                        TMP_4d[m][n][k][l] += dgam[alpha] * J_alpha[alpha][m][n][k][l]

    Q_mat = [[[[0. for i in range(3)] for j in range(3)] for k in range(3)] for l in range(3)]
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    Q_mat[i][j][k][l] = 0.0
                    for m in range(3):
                        for n in range(3):
                            Q_mat[i][j][k][l] += K_inv[i][j][m][n] * (D_mat[m][n][k][l] - TMP_4d[m][n][k][l])

    # -------------- Calculation of R_alpha -----------------
    R_alpha = [[[0. for i in range(3)] for j in range(3)] for k in range(len(schmid))]
    # GT_mat = self._GT_mat_chache
    for alpha in range(len(schmid)):
        for i in range(3):
            for j in range(3):
                R_alpha[alpha][i][j] = 0.0
                for k in range(3):
                    for l in range(3):
                        R_alpha[alpha][i][j] += GT_mat[alpha][k][l] * Q_mat[k][l][i][j]
    # ------------- Calculation of S_mat---------------------
    TMP_4d = [[[[0. for i in range(3)] for j in range(3)] for k in range(3)] for l in range(3)]
    for k in range(3):
        for l in range(3):
            for p in range(3):
                for j in range(3):
                    TMP_4d[k][l][p][j] = 0.0
                    for alpha in range(len(schmid)):
                        TMP_4d[k][l][p][j] += R_alpha[alpha][k][l] * schmid[alpha][p][j]

    TMP_2d = [[0. for i in range(3)] for j in range(3)]
    for p in range(3):
        for j in range(3):
            TMP_2d[p][j] = 0.0
            for alpha in range(len(schmid)):
                TMP_2d[p][j] += dgam[alpha] * schmid[alpha][p][j]

    S_mat = [[[[0. for i in range(3)] for j in range(3)] for k in range(3)] for l in range(3)]
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    S_mat[i][j][k][l] = Rot_el[i][k] * Fe_t[l][j]

    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    for p in range(3):
                        S_mat[i][j][k][l] -= Rot_el[i][k] * Fe_t[l][p] * TMP_2d[p][j]

    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    for m in range(3):
                        for n in range(3):
                            for p in range(3):
                                S_mat[i][j][k][l] -= Rot_el[i][m] * Strh_el[m][n] * Fe_t[n][p] * TMP_4d[k][l][p][j]

    # ------------------ Calculation of inverse of Fe ----------------
    Fe_tau_inv = [[0.0 for i in range(3)] for j in range(3)]
    Fe_tauID = inverse3.matinv3(Fe_tau)
    Fe_tau_inv = Fe_tauID[0]
    det_Fe_tau = Fe_tauID[1]

    # S_star = self._S_star_cache

    W_mat = [[[[0. for i in range(3)] for j in range(3)] for k in range(3)] for l in range(3)]
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    for m in range(3):
                        for n in range(3):
                            W_mat[i][j][k][l] += S_mat[i][m][k][l] * S_star[m][n] * Fe_tau[j][n]
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    for m in range(3):
                        for n in range(3):
                            W_mat[i][j][k][l] += Fe_tau[i][m] * Q_mat[m][n][k][l] * Fe_tau[j][n]
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    for m in range(3):
                        for n in range(3):
                            W_mat[i][j][k][l] += Fe_tau[i][m] * S_star[m][n] * S_mat[j][n][k][l]
    TMP_sf = [[0.0 for i in range(3)] for j in range(3)]
    for k in range(3):
        for l in range(3):
            for p in range(3):
                for q in range(3):
                    TMP_sf[k][l] += S_mat[p][q][k][l] * Fe_tau_inv[q][p]
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    for m in range(3):
                        for n in range(3):
                            W_mat[i][j][k][l] -= Fe_tau[i][m] * S_star[m][n] * Fe_tau[j][n] * TMP_sf[k][l]
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    W_mat[i][j][k][l] = W_mat[i][j][k][l] / det_Fe_tau

    Dep = deflate_matrix_W.reduce_wmat(W_mat)

    return Dep