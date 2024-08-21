# Crystal Plasticity Code.

import math,sys




############ making global internal forces  #############
def asembl_vec(gf_nod,InF,lnode):
    node_indx = [0 for i in range(24)]
    
    for i in range(8):
        node_indx[3*i] = 3*lnode[i] ; node_indx[3*i+1] = 3*lnode[i]+1 ; node_indx[3*i+2] = 3*lnode[i]+2

    for i in range(24):
        gf_nod[node_indx[i]] += InF[i]

    return gf_nod
############ end of making global internal forces  #############
