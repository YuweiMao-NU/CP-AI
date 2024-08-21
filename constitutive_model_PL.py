# Crystal Plasticity Code.

import math,sys





####################### Flow rule #########################
class update_statev:
    def __init__(self,tau_alpha,res_ssd,dgam,dgam_dta,prop):
        self.tau_alpha = tau_alpha
        self.res_ssd = res_ssd
        self.dgam = dgam
        self.dgam_dta = dgam_dta
        self.prop = prop

    def RES(self):
        

        const_w1 = self.prop[4] ; const_w2 = self.prop[5]
        const_ss = self.prop[6]
        const_a = self.prop[7]
        const_h0 = self.prop[8] ; const_m = self.prop[9]
        g0dot = self.prop[10]
        dt = self.prop[11]

        n_slip = len(self.tau_alpha)
        res = [0.0 for i in range(n_slip)]
        res_t = [0.0 for i in range(n_slip)]
        dgamma = [0.0 for i in range(n_slip)]
        dgam_dtau = [0.0 for i in range(n_slip)]

 #       for i in range(n_slip):
        res0 = self.res_ssd
        dgamma = self.dgam
        dgam_dtau = self.dgam_dta
    

        for k in range(n_slip):
            res[k] = 0.0
            for i in range(n_slip):
                if i == k:
                    const_qab = const_w1
                elif i != k:
                    const_qab = const_w2
                
                ratio_res=1.0-(res0[i]/const_ss)

                res[k] += const_qab*const_h0*abs(dgamma[i])*((ratio_res)**const_a)

        for k in range(n_slip):
            res[k] += res0[k]


        for k in range(n_slip):
            if res[k] >= 1.0:
                
                ratio_alpha = self.tau_alpha[k]/res[k]
                
                if self.tau_alpha[k] >= 0.0:
                    const_sign = 1.0
                elif self.tau_alpha[k] < 0.0:
                    const_sign = -1.0
    
                m_inv = 1.0/const_m
                dgamma[k] = dt*const_sign*g0dot*((abs(ratio_alpha))**m_inv)
                

######    Calculation of dgamma_dtau (The sgn(Tau_alpha) is multiplied with another sgn(Tau_alpha) from the derivative)
                res_inv = 1.0/res[k]
                
                dgam_dtau[k] = dt*res_inv*g0dot*m_inv*((abs(ratio_alpha))**(m_inv-1.0))

            elif res[k] < 1.0:
               
                dgamma[k] = 0.0
                dgam_dtau[k] = 0.0
    
        return dgamma,dgam_dtau,res

