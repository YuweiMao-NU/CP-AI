# Crystal Plasticity Code.

import math,sys



######### Inflate 1st 6 vector to 2nd 3*3 tensor ################
def inflate_vec(as_vec):
    as_star = [[0.0 for i in range(3)] for j in range(3)]
    for i in range(3):
        as_star[i][i] = as_vec[i]

    as_star[0][1] = as_vec[3]
    as_star[1][2] = as_vec[4]    
    as_star[0][2] = as_vec[5]

    as_star[1][0] = as_vec[3]
    as_star[2][1] = as_vec[4]
    as_star[2][0] = as_vec[5]

    return as_star

######### end of Inflate 1st 6 vector to 2nd 3*3 tensor ################

