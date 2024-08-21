# Crystal Plasticity Code.

import math,sys

from FE_functions import res_shear_stess
from FE_functions import calc_GN

from Mat_functions import deflate_matrix_3
from Mat_functions import deflate_matrix_W
from Mat_functions import deflate_matrix_2
from Mat_functions import inflate_vec
from Mat_functions import inflate_matrix_6
from Mat_functions import inverse6
from Mat_functions import matrix_mlutip
from Mat_functions import inverse3


from constitutive_model_PL import *





################ User Material ###########

class UMAT:
    def __init__(self,S_trial,C_mat,Fp,C_alpha,schmid,F_t,F_tau,res,dgam,dgam_dta,props,S_star0):
        self.S_trial = S_trial
        self.C_mat = C_mat
        self.Fp = Fp
        self.C_alpha = C_alpha
        self.schmid = schmid
        self.F_t = F_t
        self.F_tau = F_tau
        self.res = res
        self.dgam = dgam
        self.dgam_dta = dgam_dta
        self.props = props
        self.S_star0 = S_star0


#------------ Forth order Kronecker Delta ---------------------
    def delta_kron4(self):
        delta_kron = [[0. for ii in range(3)] for jj in range(3)]
        delta_kron[0][0] = delta_kron[1][1]= delta_kron[2][2] = 1.0

        delta_kron4d = [[[[0. for i in range(3)] for j in range(3)] for k in range(3)] for l in range(3)] 

        for i in range(3):
            for j in range(3):
                for k in range(3):
                    for l in range(3):
                        delta_kron4d[i][j][k][l] += delta_kron[i][k]*delta_kron[j][l]+delta_kron[i][l]*delta_kron[j][k]

        for i in range(3):
            for j in range(3):
                for k in range(3):
                    for l in range(3):
                        delta_kron4d[i][j][k][l] = 0.5*delta_kron4d[i][j][k][l]
        
        return delta_kron4d
#---------------------------------------------------------------

#------------------ Iteration process for time integration of Fp --------        

    def itr(self):
        n_slip = len(self.schmid)

        tau_alpha = [0.0 for k in range(n_slip)]
        tau_alpha = res_shear_stess.resolvedshear(n_slip,self.S_trial,self.schmid)

        res_ssd = self.res
        dgam = self.dgam
        dgam_dta = self.dgam_dta
        dgam = [0.0 for k in range(n_slip)]
        dgam_dta = [0.0 for k in range(n_slip)]
       
        
#---------------------- Update state variables ----------------
        prop = self.props
        s1 = update_statev(tau_alpha,res_ssd,dgam,dgam_dta,prop)
        u1 = s1.RES()
        dgam = u1[0]
        dgam_dta = u1[1]

        

        TMP1 = [[0.0 for it in range(3)] for jt in range(3)]
        for it in range(3):
            for jt in range(3):
                for kt in range(n_slip):
                    TMP1[it][jt] += dgam[kt]*self.C_alpha[kt][it][jt]


#-------------------- Calculation of first value of 2nd PKS-----------------
        S_star = [[0.0 for it in range(3)] for jt in range(3)]
        for it in range(3):
            for jt in range(3):
                S_star[it][jt] = self.S_trial[it][jt]-TMP1[it][jt]


#        print "Start of iteration loop for time integration"
#----------------- iteration loop for S_star -------------------
        niter = 20
        iitr = 0
        ratio_norm = 1.0e10
        while (iitr <= 20 and abs(ratio_norm) >= 0.001):

 #------------ Get the resolved shear stress and update state variables with S_star-----------------
            tau_alpha = [0.0 for k in range(n_slip)]
            tau_alpha = res_shear_stess.resolvedshear(n_slip,S_star,self.schmid)

            
            s2 = update_statev(tau_alpha,res_ssd,dgam,dgam_dta,prop)
            u2 = s2.RES()
            dgam = u2[0]
            dgam_dta = u2[1]
#            res_ssd = u2[2]
            
#---------
            GT_mat = [[[0. for i in range(3)] for j in range(3)] for k in range(n_slip)]

            for k in range(n_slip):
                for i in range(3):
                    for j in range(3):
                        GT_mat[k][i][j] = 0.5*(self.schmid[k][i][j]+self.schmid[k][j][i])*dgam_dta[k]
            self._GT_mat_chache = GT_mat
#-----------
#----------
            RJ_mat = [[[[0. for i in range(3)] for j in range(3)] for k in range(3)] for l in range(3)]
            for i in range(3):
                for j in range(3):
                    for m in range(3):
                        for n in range(3):                  
                            for k in range(n_slip):
                                RJ_mat[i][j][m][n] += self.C_alpha[k][i][j]*GT_mat[k][m][n]
            
#-------------
            delta4 = UMAT.delta_kron4(self)
            for i in range(3):
                for j in range(3):
                    for m in range(3):
                        for n in range(3):
                            RJ_mat[i][j][m][n] += delta4[i][j][m][n]

            RJ_reduced = deflate_matrix_3.reduce_mat(RJ_mat)
            
                        
            RJ_reduced = inverse6.matinv(RJ_reduced,6)
            
            self._RJ_reduced_cache = RJ_reduced

#--------------
#-------------
            GN_vec = calc_GN.calc_GNvec(n_slip,S_star,self.S_trial,dgam,self.C_alpha)
#----------------
#-----------------

            S_vec = deflate_matrix_2.reduce_vec(S_star)

            TMP_RG = [0.0 for i in range(6)]
            TMP_vec = [0.0 for i in range(6)]

            for i in range(6):
                for j in range(6):
                    TMP_RG[i] += RJ_reduced[i][j]*GN_vec[j]

            for i in range(6):
                TMP_vec[i] = S_vec[i]

            for i in range(6):
                S_vec[i] -= TMP_RG[i]
            

            S_star = inflate_vec.inflate_vec(S_vec)
            
         
            dot_product_S_vec = 0.0
            for i in range(6):
                dot_product_S_vec += S_vec[i]*S_vec[i]

            rnorm_s_vec = math.sqrt(dot_product_S_vec)


            dot_product_TMP_vec = 0.0
            for i in range(6):
                dot_product_TMP_vec += TMP_vec[i]*TMP_vec[i]

            rnorm_tmp_vec = math.sqrt(dot_product_TMP_vec)


            if abs(rnorm_tmp_vec) < 1.0:
                ratio_norm = 0.0
            elif abs(rnorm_tmp_vec) >= 1.0:
                ratio_norm = (rnorm_s_vec-rnorm_tmp_vec)/rnorm_tmp_vec
#--------------------

#            print iitr+1,abs(ratio_norm)
#            print

            iitr += 1

#        raw_input()

#------------ Get the resolved shear stress and update state variables with the updated S_star------------
        tau_alpha = [0.0 for k in range(n_slip)]
        tau_alpha = res_shear_stess.resolvedshear(n_slip,S_star,self.schmid)


        s3 = update_statev(tau_alpha,res_ssd,dgam,dgam_dta,prop)
        u3 = s3.RES()
        dgam = u3[0]
        dgam_dta = u3[1]
        res_ssd = u3[2]


        self._dgam_cache = dgam
        self._S_star_cache = S_star
#------------???????????????? transfer dgamma to gauss point (GNDs)----------

#---------- Calculation of Velocity Gradient --------------

        Lp = [[0.0 for i in range(3)] for j in range(3)]
        for i in range(3):
            for j in range(3):
                for k in range(n_slip):
                    Lp[i][j] += dgam[k]*self.schmid[k][i][j]
#----------- Calculation of Fp-tau ---------------------
        Iden = [[0. for i in range(3)] for j in range(3)]
        Fp_tau = [[0. for i in range(3)] for j in range(3)]
        
        Iden[0][0] = Iden[1][1] = Iden[2][2] = 1.0
        for i in range(3):
            for j in range(3):
                Lp[i][j] +=Iden[i][j]

        Fp_tau = matrix_mlutip.matmul(Lp,self.Fp)

        det_Fp_tau=(Fp_tau[0][0]*Fp_tau[1][1]*Fp_tau[2][2]\
                        -Fp_tau[0][0]*Fp_tau[1][2]*Fp_tau[2][1]-Fp_tau[1][0]*Fp_tau[0][1]*Fp_tau[2][2]\
                        +Fp_tau[1][0]*Fp_tau[0][2]*Fp_tau[2][1]+Fp_tau[2][0]*Fp_tau[0][1]*Fp_tau[1][2]\
                        -Fp_tau[2][0]*Fp_tau[0][2]*Fp_tau[1][1])

        oby3=1.0/3.0
        for i in range(3):
            for j in range(3):
                Fp_tau[i][j] = Fp_tau[i][j]/((det_Fp_tau)**oby3)

        
#----------- Calculation of the new Cauchy Stress ------------

        Fe_tau = [[0.0 for i in range(3)] for j in range(3)]
        Fp_tauI = [[0.0 for i in range(3)] for j in range(3)]
        FptauD = inverse3.matinv3(Fp_tau)
        Fp_tauI = FptauD[0] 


        for i in range(3):
            for j in range(3):
                for k in range(3):
                    Fe_tau[i][j] += self.F_tau[i][k]*Fp_tauI[k][j]

    
        TPMc = [[0.0 for i in range(3)] for j in range(3)]                
        TPMc = matrix_mlutip.matmul(S_star,list(zip(*Fe_tau)))
    
        Cauchy = [[0.0 for i in range(3)] for j in range(3)]
        Cauchy = matrix_mlutip.matmul(Fe_tau,TPMc)

        det_F_tau=(self.F_tau[0][0]*self.F_tau[1][1]*self.F_tau[2][2]\
                        -self.F_tau[0][0]*self.F_tau[1][2]*self.F_tau[2][1]-self.F_tau[1][0]*self.F_tau[0][1]*self.F_tau[2][2]\
                        +self.F_tau[1][0]*self.F_tau[0][2]*self.F_tau[2][1]+self.F_tau[2][0]*self.F_tau[0][1]*self.F_tau[1][2]\
                        -self.F_tau[2][0]*self.F_tau[0][2]*self.F_tau[1][1])

        
        for i in range(3):
            for j in range(3):
                Cauchy[i][j] = Cauchy[i][j]/det_F_tau
                

        Iden = [[0. for i in range(3)] for j in range(3)]
        TMPe = [[0. for i in range(3)] for j in range(3)]
        Strain_el = [[0. for i in range(3)] for j in range(3)]
        Iden[0][0] = Iden[1][1] = Iden[2][2] = 1.0

        TMPe = matrix_mlutip.matmul(list(zip(*Fe_tau)),Fe_tau)
        for i in range(3):
            for j in range(3):
                Strain_el[i][j] = 0.5*(TMPe[i][j]-Iden[i][j])
        
        self._Fe_tau_cache = Fe_tau
        return Fe_tau,Fp_tau,Cauchy,res_ssd,tau_alpha,dgam,dgam_dta,S_star


#        return u3
    
#------------------ End of Iteration method for time integration of Fp --------          
#--------------------------------------------------------------

#---------------------- Polar Decomposition of Rlative Gradient Deformation -------------------
    def polardecomp(self):

#-------------- Calculation of relative deformation gradient F_REL = F_tau*(F_t)^-1 -----------
        F_REL = [[0.0 for i in range(3)] for j in range(3)]
        F_tI = [[0. for i in range(3)] for j in range(3)]

        F_tID = inverse3.matinv3(self.F_t)
        F_tI = F_tID[0]


        for i in range(3):
            for j in range(3):
                F_REL[i][j] = 0.0
                for k in range(3):
                    F_REL[i][j] += self.F_tau[i][k]*F_tI[k][j]
           
    
#----------------------- square root of a positive matrix U=sqrt(C) ------------

        R = [[0. for i in range(3)] for j in range(3)]
        C = [[0. for i in range(3)] for j in range(3)]
        Csquare = [[0. for i in range(3)] for j in range(3)]
        Iden = [[0. for i in range(3)] for j in range(3)]
        U = [[0. for i in range(3)] for j in range(3)]
        invU = [[0. for i in range(3)] for j in range(3)]

        Iden[0][0] = Iden[1][1] = Iden[2][2] = 1.0

        F_REL_T = list(zip(*F_REL))
#----------------------- C = FTF --------------------------------------
        for i in range(3):
                for j in range(3):
                    for k in range(3):
                        C[i][j] += F_REL_T[i][k]*F_REL[k][j]

#---------------------- C^2 ------------------------------------------
        for i in range(3):
                for j in range(3):
                    for k in range(3):
                        Csquare[i][j] += C[i][k]*C[k][j]


        I_C = C[0][0]+C[1][1]+C[2][2]
    
#--------------------- Invarients of C -------------------------------
        I_Csquare = Csquare[0][0]+Csquare[1][1]+Csquare[2][2] 
        II_C = 0.5*(I_C**2-I_Csquare)
        III_C = C[0][0]*(C[1][1]*C[2][2]-C[1][2]*C[2][1])-\
            C[0][1]*(C[1][0]*C[2][2]-C[2][0]*C[1][2])+\
            C[0][2]*(C[1][0]*C[2][1]-C[1][1]*C[2][0])

        k = I_C**2-3.0*II_C
        
        l = I_C*(I_C**2-4.5*II_C)+13.5*III_C

######## added for initialization ########## 
        
        if abs(k) < 1.0e-10:
            k = 1.0e-10


###########################################

        t1 = l/(k**1.5)

        if t1<-1.0:
            t1 = -1.0
        elif t1>1.0:
            t1 = 1.0
            
        teta = (1.0/3.0)*math.acos(t1)
            
        alpha = 2.0*math.pi/3.0

    #--------------------- Eigenvalues of U ----------------------------- 
        lamda1 = math.sqrt((1.0/3.0)*(I_C+2.0*math.sqrt(k)*math.cos(teta)))
        lamda2 = math.sqrt((1.0/3.0)*(I_C+2.0*math.sqrt(k)*math.cos(alpha+teta)))
        lamda3 = math.sqrt((1.0/3.0)*(I_C+2.0*math.sqrt(k)*math.cos(2.0*alpha+teta)))

        

    #--------------------- Invarients of U ---------------------------
        I_U = lamda1+lamda2+lamda3
        II_U = lamda1*lamda2+lamda1*lamda3+lamda2*lamda3
        III_U = lamda1*lamda2*lamda3

    #--------------------- U and Inverse of U ------------------------

        for i in range(3):
            for j in range(3):
                U[i][j] = (1.0/(I_U*II_U-III_U))*(I_U*III_U*Iden[i][j]+\
                                                                     (I_U**2-II_U)*C[i][j]-Csquare[i][j])
                invU[i][j] = (1.0/III_U)*(II_U*Iden[i][j]-\
                                                                     I_U*U[i][j]+C[i][j])
#----------------------- R = FU^-1 ----------------------------- 
        for i in range(3):
                for j in range(3):
                    for k in range(3):
                        R[i][j] += F_REL[i][k]*invU[k][j]

#---------------------- End of Polar Decomposition ---------------


        return U,R

#---------------- Calculation of elasto plastic modulus or elastic Jacobian ---------

    def W_mat(self):
#---------Elastic stretch and rotation------------ 
        Poldecomp = UMAT.polardecomp(self)
        Strh_el = Poldecomp[0]
        Rot_el = Poldecomp[1]
#------------------------------------------------
        Fe_tau = self._Fe_tau_cache
        
        
        Fp_t_inv = [[0.0 for i in range(3)] for j in range(3)]
        Fp_tID = inverse3.matinv3(self.Fp)
        Fp_t_inv = Fp_tID[0]
        Fe_t = [[0.0 for i in range(3)] for j in range(3)]
        for i in range(3):
            for j in range(3):
                Fe_t[i][j] = 0.0
                for k in range(3):
                    Fe_t[i][j] += self.F_t[i][k]*Fp_t_inv[k][j]

        
#------ Calculation of L_mat-------
        L_mat = [[[[0. for i in range(3)] for j in range(3)] for k in range(3)] for l in range(3)] 
        
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    for l in range(3):
                        L_mat[i][j][k][l] = 0.0
                        for m in range(3):
                            L_mat[i][j][k][l] += Fe_t[k][i]*Strh_el[l][m]*Fe_t[m][j]+Fe_t[m][i]*Strh_el[m][k]*Fe_t[l][j]

#--------- Calculation of D_mat ---------
        D_mat = [[[[0. for i in range(3)] for j in range(3)] for k in range(3)] for l in range(3)] 
        
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    for l in range(3):
                        D_mat[i][j][k][l] = 0.0
                        for m in range(3):
                            for n in range(3):
                                D_mat[i][j][k][l] += 0.5*self.C_mat[i][j][m][n]*L_mat[m][n][k][l]
 
#-------------- Calculation of G_alpha----------------------
        G_alpha = [[[[[0. for i in range(3)] for j in range(3)] for k in range(3)] for l in range(3)] for alpha in range(len(self.schmid))]
                        
        for alpha in range(len(self.schmid)):
            for m in range(3):
                for n in range(3):
                    for k in range(3):
                        for l in range(3):
                            G_alpha[alpha][m][n][k][l] = 0.0
                            for p in range(3):
                                G_alpha[alpha][m][n][k][l] += L_mat[m][p][k][l]*self.schmid[alpha][p][n]+self.schmid[alpha][m][p]*L_mat[p][n][k][l]
#--------- Calculation of J_alpha ---------
        J_alpha = [[[[[0. for i in range(3)] for j in range(3)] for k in range(3)] for l in range(3)] for alpha in range(len(self.schmid))]
        
        for alpha in range(len(self.schmid)):
            for i in range(3):
                for j in range(3):
                    for k in range(3):
                        for l in range(3):
                            J_alpha[alpha][i][j][k][l] = 0.0
                            for m in range(3):
                                for n in range(3):
                                    J_alpha[alpha][i][j][k][l] += 0.5*self.C_mat[i][j][m][n]*G_alpha[alpha][m][n][k][l]

#------------- Calculation of Q_mat -------------------  
        RJ_reduced = self._RJ_reduced_cache
        K_inv = inflate_matrix_6.inflate_ten(RJ_reduced)
    
        dgam = self._dgam_cache

        TMP_4d = [[[[0. for i in range(3)] for j in range(3)] for k in range(3)] for l in range(3)]
        for m in range(3):
            for n in range(3):
                for k in range(3):
                    for l in range(3):
                        TMP_4d[m][n][k][l] = 0.0
                        for alpha in range(len(self.schmid)):
                            TMP_4d[m][n][k][l] += dgam[alpha]*J_alpha[alpha][m][n][k][l]

        Q_mat = [[[[0. for i in range(3)] for j in range(3)] for k in range(3)] for l in range(3)]
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    for l in range(3):
                        Q_mat[i][j][k][l] = 0.0
                        for m in range(3):
                            for n in range(3):
                                Q_mat[i][j][k][l] += K_inv[i][j][m][n]*(D_mat[m][n][k][l]-TMP_4d[m][n][k][l])


#-------------- Calculation of R_alpha -----------------
        R_alpha = [[[0. for i in range(3)] for j in range(3)] for k in range(len(self.schmid))]
        GT_mat = self._GT_mat_chache
        for alpha in range(len(self.schmid)):
            for i in range(3):
                for j in range(3):
                    R_alpha[alpha][i][j] = 0.0
                    for k in range(3):
                        for l in range(3):
                            R_alpha[alpha][i][j] += GT_mat[alpha][k][l]*Q_mat[k][l][i][j]
#------------- Calculation of S_mat---------------------
        TMP_4d = [[[[0. for i in range(3)] for j in range(3)] for k in range(3)] for l in range(3)]
        for k in range(3):
            for l in range(3):
                for p in range(3):
                    for j in range(3):
                        TMP_4d[k][l][p][j] = 0.0
                        for alpha in range(len(self.schmid)):
                            TMP_4d[k][l][p][j] += R_alpha[alpha][k][l]*self.schmid[alpha][p][j]

        TMP_2d = [[0. for i in range(3)] for j in range(3)]
        for p in range(3):
            for j in range(3):
                TMP_2d[p][j] = 0.0
                for alpha in range(len(self.schmid)):
                    TMP_2d[p][j] += dgam[alpha]*self.schmid[alpha][p][j]

        S_mat = [[[[0. for i in range(3)] for j in range(3)] for k in range(3)] for l in range(3)]
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    for l in range(3):
                        S_mat[i][j][k][l] = Rot_el[i][k]*Fe_t[l][j]

        for i in range(3):
            for j in range(3):
                for k in range(3):
                    for l in range(3):
                        for p in range(3):
                            S_mat[i][j][k][l] -= Rot_el[i][k]*Fe_t[l][p]*TMP_2d[p][j]

        for i in range(3):
            for j in range(3):
                for k in range(3):
                    for l in range(3):
                        for m in range(3):
                            for n in range(3):
                                for p in range(3):
                                    S_mat[i][j][k][l] -= Rot_el[i][m]*Strh_el[m][n]*Fe_t[n][p]*TMP_4d[k][l][p][j]


#------------------ Calculation of inverse of Fe ----------------
        Fe_tau_inv = [[0.0 for i in range(3)] for j in range(3)]
        Fe_tauID = inverse3.matinv3(Fe_tau)
        Fe_tau_inv = Fe_tauID[0] ; det_Fe_tau = Fe_tauID[1]

        S_star = self._S_star_cache
        
 

        W_mat = [[[[0. for i in range(3)] for j in range(3)] for k in range(3)] for l in range(3)]
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    for l in range(3):
                        for m in range(3):
                            for n in range(3):
                                W_mat[i][j][k][l] += S_mat[i][m][k][l]*S_star[m][n]*Fe_tau[j][n]
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    for l in range(3):
                        for m in range(3):
                            for n in range(3):
                                W_mat[i][j][k][l] += Fe_tau[i][m]*Q_mat[m][n][k][l]*Fe_tau[j][n]
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    for l in range(3):
                        for m in range(3):
                            for n in range(3):
                                W_mat[i][j][k][l] += Fe_tau[i][m]*S_star[m][n]*S_mat[j][n][k][l]
        TMP_sf = [[0.0 for i in range(3)] for j in range(3)]
        for k in range(3):
            for l in range(3):
                for p in range(3):
                    for q in range(3):
                        TMP_sf[k][l] += S_mat[p][q][k][l]*Fe_tau_inv[q][p]
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    for l in range(3):
                        for m in range(3):
                            for n in range(3):
                                W_mat[i][j][k][l] -= Fe_tau[i][m]*S_star[m][n]*Fe_tau[j][n]*TMP_sf[k][l]
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    for l in range(3):
                        W_mat[i][j][k][l] = W_mat[i][j][k][l]/det_Fe_tau


        Dep = deflate_matrix_W.reduce_wmat(W_mat)
	
        return Dep



