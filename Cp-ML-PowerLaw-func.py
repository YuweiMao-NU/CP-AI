# Crystal Plasticity Code.

import math,sys
from pylab import *
import matplotlib.pyplot as plt

from solver_init_CP import *
from solver_ML import *
import pickle
import numpy as np
import sklearn

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++ reading the base input file for init_ML_step ++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

f = open('py_init.15','r')
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
angle = []
for i in range(nelem):
    angle.append([float(x) for x in f1.readline().split()])

#---------------- end of reading input file --------------------

ssy = open('SS.txt','w')

ssy.write("%f   %f" % (0.0,0.0))
ssy.write('\n')


m_init_CP = Mesh_init(nnode,tnode,coors,nelem,elements)
ga = m_init_CP.solve_init_CP(nbc,bcs,angle,props,tstep,nconv,n_slip)

#--------------- volume stress-strain--------------
x = [] ; y = []
for istep in range(init_ML_step):
    x.append(ga[0][istep][1])
    y.append(ga[1][istep][1]/1.0e6)
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

#***********************************************************************************#
#***********************************************************************************#
#************ Reading the input file for the actual microstructure *****************#
#***********************************************************************************#
#***********************************************************************************#

f = open('py.15','r')
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
nstep = tstep[0]
props.append(tstep[1])
#----------------------------------------------
f1 = open('texture.15','r')
angle = []
for i in range(nelem):
    angle.append([float(x) for x in f1.readline().split()])

#---------------- end of reading input file --------------------

#ssy = open('SS.txt','w')

# load ML model and parameters
input_size = 50
model = pickle.load(open('Ridge.model', 'rb'))

# print(garray_list)
Fp_list = np.array(Fp_list).reshape(init_ML_step, 8, 8, 3, 3)
F_list = np.array(F_list).reshape(init_ML_step, 8, 8, 3, 3)
# print(Fp_list[-3:, 0, 0, :, :].reshape(-1, 9))

# datapath = ''
# np.savetxt(datapath + 'F.npy', F_list.reshape(-1, 9))
# np.savetxt(datapath + 'Fp.npy', Fp_list.reshape(-1, 9))


# print(Fp_list.shape, F_list.shape)
Fp_list = Fp_list.tolist()
F_list = F_list.tolist()
m = Mesh(nnode,tnode,coors,nelem,elements)
ga = m.solve(input_size, model, nbc, bcs, angle, props, tstep, nconv,n_slip,init_ML_step,Fp_list,F_list,K_i,Load_i,garray_list,Cauchy_i,gdisp)

#--------------- volume stress-strain--------------
#x = [] ; y = []

#for istep in range(init_ML_step,nstep,1):

for istep in range(nstep-init_ML_step-1):
    x.append(ga[0][istep+1][1])
    y.append(ga[1][istep+1][1]/1.0e6)
    ssy.write("%f   %f" % (ga[0][istep+1][1],ga[1][istep+1][1]/1.0e6))
    ssy.write('\n')

ssy.close()
#--------------- end of volumeC stress-strain--------------

# plot(x,y)
# ylabel('Stress(MPa)')
# xlabel('Strain')
#
# show()


# save the figure instead of showing it
# plt.savefig('stress_strain_plot.png')  # save the plot to a file

#***********************************************************************************#
#***********************************************************************************#
#**************** End of simulations for the actual microstructure *****************#
#***********************************************************************************#
#***********************************************************************************#

