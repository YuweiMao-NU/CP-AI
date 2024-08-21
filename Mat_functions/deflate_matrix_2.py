# Crystal Plasticity Code.

import math,sys




#################  Reduce 2nd order 3*3 tensor to 1st 6 vector ########

def reduce_vec(as_star):
    as_vec = [0.0 for i in range(6)]
    for i in range(3):
        as_vec[i] = as_star[i][i]
    
    as_vec[3] = 0.5*(as_star[0][1]+as_star[1][0])
    as_vec[4] = 0.5*(as_star[1][2]+as_star[2][1])
    as_vec[5] = 0.5*(as_star[0][2]+as_star[2][0])

    return as_vec

#################  end of Reduce 2nd order 3*3 tensor to 1st 6 vector ########


