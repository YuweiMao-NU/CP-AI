# Crystal Plasticity Code.

import math,sys



################  Resolved Shear Stress ##########
def resolvedshear(n_slip,stress,schmidfactor):
    tau_alpha = [0.0 for k in range(n_slip)]
    for k in range(n_slip):
        for i in range(3):
            for j in range(3):
                tau_alpha[k] += stress[i][j]*schmidfactor[k][i][j]
 
    return tau_alpha

#################### End of Resolved Shear Stress ##############

