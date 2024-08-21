# Crystal Plasticity Code.

import math,sys
from pylab import *

from solver_init_CP import *
from solver_ML import *
import pickle
import numpy as np
import sklearn
import csv

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++ reading the base input file for init_ML_step ++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

f = open('data_py_init.15','r')
a = [int(x) for x in f.readline().split()]

elements = []
nnode = a[0] ; nelem = a[1] ; tnode = a[2] ; nbc = a[3]
nconv = a[4]-1  ; n_slip = a[5]
for i in range(nelem):
    elements.append([int(x)-1 for x in f.readline().split()])

coors = []
for i in range(tnode):
    b = [float(x) for x in f.readline().split()]
    b[0] = int(b[0])-1
    coors.append(b)

bcs = []
for i in range(nbc):
    bcs.append([int(x) for x in f.readline().split()])

props = [float(x) for x in f.readline().split()]

tstep = [float(x) for x in f.readline().split()]
tstep[0] = int(tstep[0])
init_ML_step = tstep[0]
props.append(tstep[1])
#----------------------------------------------
f1 = open('texture_init.15','r')
first_list = range(0, 90, 10)
second_list = range(0, 90, 10)
for first in first_list:
    for second in second_list:
        try:
            angle = [[first, second, 0] for _ in range(8)]
            # for i in range(nelem):
            #     angle.append([float(x) for x in f1.readline().split()])
            print(angle)

    #---------------- end of reading input file --------------------

            ssy = open('SS.txt','w')

            ssy.write("%f   %f" % (0.0,0.0))
            ssy.write('\n')




            m_init_CP = Mesh_init(nnode,tnode,coors,nelem,elements)
            ga = m_init_CP.solve_init_CP(nbc,bcs,angle,props,tstep,nconv,n_slip)




            #--------------- volume stress-strain--------------
            x = [] ; y = []
            for istep in range(init_ML_step):
                x.append(ga[0][istep][1]) ; y.append(ga[1][istep][1]/1.0e6)
                ssy.write("%f   %f" % (ga[0][istep][1],ga[1][istep][1]/1.0e6))
                ssy.write('\n')

            #ssy_base.close()
            #--------------- end of volume stress-strain--------------


            #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
            #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
            #++++++++++++++ End of simulations for the base microstructure ++++++++++++++++++++#
            #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
            #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

            #### At this point the full history of Ft and Fp will be avialble
            # Fp_list and F_list is a 3D list witt these diementions Fp_list[init_ML_step][element=1][ngauss=8]

            gdisp = ga[2]
            Fp_list = ga[4]
            F_list = ga[5]
            K_i = ga[6]
            garray_list = ga[7]
            Cauchy_i = ga[8]
            Load_i = init_ML_step*tstep[1]*tstep[2]

            # print(garray_list)
            Fp_list = np.array(Fp_list).reshape(init_ML_step, 8, 8, 3, 3)
            F_list = np.array(F_list).reshape(init_ML_step, 8, 8, 3, 3)

            array1 = np.array(x)
            array2 = np.array(y)

            combined_array = np.concatenate((array1, array2), axis=0)

            datapath = 'data/'
            file_path = datapath + str(first) + '_' + str(second)+ '_' + 'stress-strain.csv'

            with open(file_path, 'w', newline='') as csvfile:
                writer = csv.writer(csvfile)
                # writer.writerow(['Combined Array'])
                writer.writerow(combined_array)


            np.savetxt(datapath + str(first) + '_' + str(second)+ '_' + 'F2000.npy', F_list.reshape(-1, 9))
            np.savetxt(datapath + str(first) + '_' + str(second)+ '_' + 'Fp2000.npy', Fp_list.reshape(-1, 9))
        except:
            print('error ', first, second)
