# Crystal Plasticity Code.

import math,sys

from FE_functions import res_shear_stess

from Mat_functions import deflate_matrix_3
from Mat_functions import deflate_matrix_2
from Mat_functions import inflate_vec
from Mat_functions import inflate_matrix_6
from Mat_functions import inverse6
from Mat_functions import matrix_mlutip
from Mat_functions import inverse3


from constitutive_model_PL import *

def update_constitutive(schmid,res,props,S_star):

    n_slip = len(schmid)

    tau_alpha = [0.0 for k in range(n_slip)]
    tau_alpha = res_shear_stess.resolvedshear(n_slip,S_star,schmid)

    res_ssd = res
    dgam = [0.0 for k in range(n_slip)]
    dgam_dta = [0.0 for k in range(n_slip)]
   
        
#---------------------- Update state variables ----------------
    prop = props
    s1 = update_statev(tau_alpha,res_ssd,dgam,dgam_dta,prop)
    u1 = s1.RES()
    dgam = u1[0]
    dgam_dta = u1[1]

    return res,dgam,dgam_dta


