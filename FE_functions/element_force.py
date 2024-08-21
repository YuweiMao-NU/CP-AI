# Crystal Plasticity Code.

import math,sys
from FE_functions import gauss
from FE_functions import shp3d



################### element reaction or internal forces ##########
def element_re(ngauss,nnode,xyz,strs):

    InF = [0.0 for i in range(nnode*3)]

    for ig in range(ngauss):
        gaussw = gauss.pgauss()
        SHP = shp3d.SHP3D(gaussw[0][ig],gaussw[1][ig],gaussw[2][ig],xyz)
        xsJ = SHP[1]
        J1 = 0
        for ind in range(nnode):

            InF[J1] += (SHP[0][0][ind]*strs[ig][0]+SHP[0][1][ind]*strs[ig][3]+SHP[0][2][ind]*strs[ig][5])*xsJ

            InF[J1+1] += (SHP[0][1][ind]*strs[ig][1]+SHP[0][0][ind]*strs[ig][3]+SHP[0][2][ind]*strs[ig][4])*xsJ

            InF[J1+2] += (SHP[0][2][ind]*strs[ig][2]+SHP[0][1][ind]*strs[ig][4]+SHP[0][0][ind]*strs[ig][5])*xsJ

            J1 += 3

    return InF
######################### end of element reaction forces #######################

