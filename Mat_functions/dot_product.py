# Crystal Plasticity Code.

import math,sys



############ dot product of a vector  #############
def dot_product(gu,neq):
    sf = 0.0
    for i in range(neq):
        sf += gu[i]*gu[i]

    return sf

############ end of dot product of a vector  #############


