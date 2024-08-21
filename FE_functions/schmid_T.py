# Crystal Plasticity Code.

import math,sys




###########################################################################
##################### Schmid Tensor for FCC materials #####################
###########################################################################

def calc_schmid(crystal_type,n_slip,qrot):
    schmid = [[[0. for i in range(3)] for j in range(3)] for k in range(n_slip)]

    r1=1.0
    r2=1.0/math.sqrt(2.0)
    r3=1.0/math.sqrt(3.0)
    r6=1.0/math.sqrt(6.0)
 
    if crystal_type == 'fcc':

        cns = [[0. for ii in range(3)] for jj in range(n_slip)]
        cms = [[0. for ii in range(3)] for jj in range(n_slip)]
        vec1 = [[0. for ii in range(3)] for jj in range(1)]
        vec2 = [[0. for ii in range(1)] for jj in range(3)]

#---------------fcc slip system 1----------
    ind_slip = 0

    cns[ind_slip][0] = r3
    cns[ind_slip][1] = r3
    cns[ind_slip][2] = r3

    cms[ind_slip][0] = 0.0
    cms[ind_slip][1] = r2
    cms[ind_slip][2] = -r2

    for ij in range(3):
       vec2[ij][0] = cns[ind_slip][ij]


    for ik in range(3):
        vec1[0][ik] = cms[ind_slip][ik]

    tmp_v1 = [[0. for ii in range(3)] for jj in range(1)]
    tmp_v2 = [[0. for ii in range(1)] for jj in range(3)]
    tmp1 = [[0. for ii in range(3)] for jj in range(3)]

    for i in range(3):
        for j in range(3):
            tmp_v1[0][i] += qrot[i][j]*vec1[0][j]
            tmp_v2[i][0] += qrot[i][j]*vec2[j][0]

    for i in range(3):
        for j in range(3):
            tmp1[i][j] = tmp_v1[0][i]*tmp_v2[j][0]
        
    for i in range(3):
        for j in range(3):
            schmid[ind_slip][i][j] = tmp1[i][j]

    #---------------fcc slip system 2----------
    ind_slip = 1

    cns[ind_slip][0] = r3
    cns[ind_slip][1] = r3
    cns[ind_slip][2] = r3

    cms[ind_slip][0] = -r2
    cms[ind_slip][1] = 0.0
    cms[ind_slip][2] = r2

    for ij in range(3):
       vec2[ij][0] = cns[ind_slip][ij]


    for ik in range(3):
        vec1[0][ik] = cms[ind_slip][ik]

    tmp_v1 = [[0. for ii in range(3)] for jj in range(1)]
    tmp_v2 = [[0. for ii in range(1)] for jj in range(3)]
    tmp1 = [[0. for ii in range(3)] for jj in range(3)]

    for i in range(3):
        for j in range(3):
            tmp_v1[0][i] += qrot[i][j]*vec1[0][j]
            tmp_v2[i][0] += qrot[i][j]*vec2[j][0]


    for i in range(3):
        for j in range(3):
            tmp1[i][j] = tmp_v1[0][i]*tmp_v2[j][0]
        
    for i in range(3):
        for j in range(3):
            schmid[ind_slip][i][j] = tmp1[i][j]
                 
    #---------------fcc slip system 3----------
    ind_slip = 2

    cns[ind_slip][0] = r3
    cns[ind_slip][1] = r3
    cns[ind_slip][2] = r3

    cms[ind_slip][0] = r2
    cms[ind_slip][1] = -r2
    cms[ind_slip][2] = 0.0

    for ij in range(3):
       vec2[ij][0] = cns[ind_slip][ij]


    for ik in range(3):
        vec1[0][ik] = cms[ind_slip][ik]

    tmp_v1 = [[0. for ii in range(3)] for jj in range(1)]
    tmp_v2 = [[0. for ii in range(1)] for jj in range(3)]
    tmp1 = [[0. for ii in range(3)] for jj in range(3)]

    for i in range(3):
        for j in range(3):
            tmp_v1[0][i] += qrot[i][j]*vec1[0][j]
            tmp_v2[i][0] += qrot[i][j]*vec2[j][0]

    for i in range(3):
        for j in range(3):
            tmp1[i][j] = tmp_v1[0][i]*tmp_v2[j][0]
        
    for i in range(3):
        for j in range(3):
            schmid[ind_slip][i][j] = tmp1[i][j]
                 
    #---------------fcc slip system 4----------
    ind_slip = 3

    cns[ind_slip][0] = r3
    cns[ind_slip][1] = -r3
    cns[ind_slip][2] = -r3

    cms[ind_slip][0] = 0.0
    cms[ind_slip][1] = -r2
    cms[ind_slip][2] = r2

    for ij in range(3):
       vec2[ij][0] = cns[ind_slip][ij]


    for ik in range(3):
        vec1[0][ik] = cms[ind_slip][ik]

    tmp_v1 = [[0. for ii in range(3)] for jj in range(1)]
    tmp_v2 = [[0. for ii in range(1)] for jj in range(3)]
    tmp1 = [[0. for ii in range(3)] for jj in range(3)]

    for i in range(3):
        for j in range(3):
            tmp_v1[0][i] += qrot[i][j]*vec1[0][j]
            tmp_v2[i][0] += qrot[i][j]*vec2[j][0]

    for i in range(3):
        for j in range(3):
            tmp1[i][j] = tmp_v1[0][i]*tmp_v2[j][0]
        
    for i in range(3):
        for j in range(3):
            schmid[ind_slip][i][j] = tmp1[i][j]
                 
    #---------------fcc slip system 5----------
    ind_slip = 4

    cns[ind_slip][0] = r3
    cns[ind_slip][1] = -r3
    cns[ind_slip][2] = -r3

    cms[ind_slip][0] = -r2
    cms[ind_slip][1] = 0.0
    cms[ind_slip][2] = -r2

    for ij in range(3):
       vec2[ij][0] = cns[ind_slip][ij]


    for ik in range(3):
        vec1[0][ik] = cms[ind_slip][ik]

    tmp_v1 = [[0. for ii in range(3)] for jj in range(1)]
    tmp_v2 = [[0. for ii in range(1)] for jj in range(3)]
    tmp1 = [[0. for ii in range(3)] for jj in range(3)]

    for i in range(3):
        for j in range(3):
            tmp_v1[0][i] += qrot[i][j]*vec1[0][j]
            tmp_v2[i][0] += qrot[i][j]*vec2[j][0]

    for i in range(3):
        for j in range(3):
            tmp1[i][j] = tmp_v1[0][i]*tmp_v2[j][0]
        
    for i in range(3):
        for j in range(3):
            schmid[ind_slip][i][j] = tmp1[i][j]

    #---------------fcc slip system 6----------
    ind_slip = 5

    cns[ind_slip][0] = r3
    cns[ind_slip][1] = -r3
    cns[ind_slip][2] = -r3

    cms[ind_slip][0] = r2
    cms[ind_slip][1] = r2
    cms[ind_slip][2] = 0.0

    for ij in range(3):
       vec2[ij][0] = cns[ind_slip][ij]


    for ik in range(3):
        vec1[0][ik] = cms[ind_slip][ik]


    tmp_v1 = [[0. for ii in range(3)] for jj in range(1)]
    tmp_v2 = [[0. for ii in range(1)] for jj in range(3)]
    tmp1 = [[0. for ii in range(3)] for jj in range(3)]

    for i in range(3):
        for j in range(3):
            tmp_v1[0][i] += qrot[i][j]*vec1[0][j]
            tmp_v2[i][0] += qrot[i][j]*vec2[j][0]

    for i in range(3):
        for j in range(3):
            tmp1[i][j] = tmp_v1[0][i]*tmp_v2[j][0]
        
    for i in range(3):
        for j in range(3):
            schmid[ind_slip][i][j] = tmp1[i][j]

    #---------------fcc slip system 7----------
    ind_slip = 6

    cns[ind_slip][0] = -r3
    cns[ind_slip][1] = r3
    cns[ind_slip][2] = -r3

    cms[ind_slip][0] = 0.0
    cms[ind_slip][1] = r2
    cms[ind_slip][2] = r2

    for ij in range(3):
       vec2[ij][0] = cns[ind_slip][ij]

       
    for ik in range(3):
        vec1[0][ik] = cms[ind_slip][ik]

    tmp_v1 = [[0. for ii in range(3)] for jj in range(1)]
    tmp_v2 = [[0. for ii in range(1)] for jj in range(3)]
    tmp1 = [[0. for ii in range(3)] for jj in range(3)]

    for i in range(3):
        for j in range(3):
            tmp_v1[0][i] += qrot[i][j]*vec1[0][j]
            tmp_v2[i][0] += qrot[i][j]*vec2[j][0]

    for i in range(3):
        for j in range(3):
            tmp1[i][j] = tmp_v1[0][i]*tmp_v2[j][0]
        
    for i in range(3):
        for j in range(3):
            schmid[ind_slip][i][j] = tmp1[i][j]
                 
    #---------------fcc slip system 8----------
    ind_slip = 7

    cns[ind_slip][0] = -r3
    cns[ind_slip][1] = r3
    cns[ind_slip][2] = -r3

    cms[ind_slip][0] = r2
    cms[ind_slip][1] = 0.0
    cms[ind_slip][2] = -r2

    for ij in range(3):
       vec2[ij][0] = cns[ind_slip][ij]


    for ik in range(3):
        vec1[0][ik] = cms[ind_slip][ik]

    tmp_v1 = [[0. for ii in range(3)] for jj in range(1)]
    tmp_v2 = [[0. for ii in range(1)] for jj in range(3)]
    tmp1 = [[0. for ii in range(3)] for jj in range(3)]

    for i in range(3):
        for j in range(3):
            tmp_v1[0][i] += qrot[i][j]*vec1[0][j]
            tmp_v2[i][0] += qrot[i][j]*vec2[j][0]

    for i in range(3):
        for j in range(3):
            tmp1[i][j] = tmp_v1[0][i]*tmp_v2[j][0]
        
    for i in range(3):
        for j in range(3):
            schmid[ind_slip][i][j] = tmp1[i][j]

    #---------------fcc slip system 9----------
    ind_slip = 8

    cns[ind_slip][0] = -r3
    cns[ind_slip][1] = r3
    cns[ind_slip][2] = -r3

    cms[ind_slip][0] = -r2
    cms[ind_slip][1] = -r2
    cms[ind_slip][2] = 0.0

    for ij in range(3):
       vec2[ij][0] = cns[ind_slip][ij]


    for ik in range(3):
        vec1[0][ik] = cms[ind_slip][ik]

    tmp_v1 = [[0. for ii in range(3)] for jj in range(1)]
    tmp_v2 = [[0. for ii in range(1)] for jj in range(3)]
    tmp1 = [[0. for ii in range(3)] for jj in range(3)]

    for i in range(3):
        for j in range(3):
            tmp_v1[0][i] += qrot[i][j]*vec1[0][j]
            tmp_v2[i][0] += qrot[i][j]*vec2[j][0]

    for i in range(3):
        for j in range(3):
            tmp1[i][j] = tmp_v1[0][i]*tmp_v2[j][0]
        
    for i in range(3):
        for j in range(3):
            schmid[ind_slip][i][j] = tmp1[i][j]

    #---------------fcc slip system 10----------
    ind_slip = 9

    cns[ind_slip][0] = -r3
    cns[ind_slip][1] = -r3
    cns[ind_slip][2] = r3

    cms[ind_slip][0] = 0.0
    cms[ind_slip][1] = -r2
    cms[ind_slip][2] = -r2

    for ij in range(3):
       vec2[ij][0] = cns[ind_slip][ij]


    for ik in range(3):
        vec1[0][ik] = cms[ind_slip][ik]

    tmp_v1 = [[0. for ii in range(3)] for jj in range(1)]
    tmp_v2 = [[0. for ii in range(1)] for jj in range(3)]
    tmp1 = [[0. for ii in range(3)] for jj in range(3)]

    for i in range(3):
        for j in range(3):
            tmp_v1[0][i] += qrot[i][j]*vec1[0][j]
            tmp_v2[i][0] += qrot[i][j]*vec2[j][0]

    for i in range(3):
        for j in range(3):
            tmp1[i][j] = tmp_v1[0][i]*tmp_v2[j][0]
        
    for i in range(3):
        for j in range(3):
            schmid[ind_slip][i][j] = tmp1[i][j]

    #---------------fcc slip system 11----------
    ind_slip = 10

    cns[ind_slip][0] = -r3
    cns[ind_slip][1] = -r3
    cns[ind_slip][2] = r3

    cms[ind_slip][0] = r2
    cms[ind_slip][1] = 0.0
    cms[ind_slip][2] = r2

    for ij in range(3):
       vec2[ij][0] = cns[ind_slip][ij]


    for ik in range(3):
        vec1[0][ik] = cms[ind_slip][ik]

    tmp_v1 = [[0. for ii in range(3)] for jj in range(1)]
    tmp_v2 = [[0. for ii in range(1)] for jj in range(3)]
    tmp1 = [[0. for ii in range(3)] for jj in range(3)]

    for i in range(3):
        for j in range(3):
            tmp_v1[0][i] += qrot[i][j]*vec1[0][j]
            tmp_v2[i][0] += qrot[i][j]*vec2[j][0]

    for i in range(3):
        for j in range(3):
            tmp1[i][j] = tmp_v1[0][i]*tmp_v2[j][0]
        
    for i in range(3):
        for j in range(3):
            schmid[ind_slip][i][j] = tmp1[i][j]

    #---------------fcc slip system 12----------
    ind_slip = 11

    cns[ind_slip][0] = -r3
    cns[ind_slip][1] = -r3
    cns[ind_slip][2] = r3

    cms[ind_slip][0] = -r2
    cms[ind_slip][1] = r2
    cms[ind_slip][2] = 0.0

    for ij in range(3):
       vec2[ij][0] = cns[ind_slip][ij]


    for ik in range(3):
        vec1[0][ik] = cms[ind_slip][ik]

    tmp_v1 = [[0. for ii in range(3)] for jj in range(1)]
    tmp_v2 = [[0. for ii in range(1)] for jj in range(3)]
    tmp1 = [[0. for ii in range(3)] for jj in range(3)]

    for i in range(3):
        for j in range(3):
            tmp_v1[0][i] += qrot[i][j]*vec1[0][j]
            tmp_v2[i][0] += qrot[i][j]*vec2[j][0]

    for i in range(3):
        for j in range(3):
            tmp1[i][j] = tmp_v1[0][i]*tmp_v2[j][0]
        
    for i in range(3):
        for j in range(3):
            schmid[ind_slip][i][j] = tmp1[i][j]
    
    return schmid

############## End of Schmid Tensor for FCC materials #####################


