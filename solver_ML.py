# Crystal Plasticity Code.

from FE_functions import euler_rotation
from FE_functions import schmid_T
from FE_functions import ciklj
from FE_functions import gauss
from FE_functions import shp3d
from FE_functions import def_grad_F
from FE_functions import calc_ABC
from FE_functions import element_force
from FE_functions import stiffness_matrix_24
from FE_functions import user_print

from Mat_functions import asembl_force
from Mat_functions import dot_product

import numpy as np

import solver

from node import *
from element import *
from umat_CP import *
from array import *
from solver_init_CP import *
from umat_ML import *
from Cauchy_ML import *
from W_mat import *
from update_constitutive import *



class Mesh:
    def __init__(self,nnode,tnode,coors,nelem,elements):
        self.nodelist = []
        self.ellist = []
        
        for inode in range(tnode):
            node = Node(coors[inode][0],coors[inode][1],coors[inode][2],coors[inode][3])
            self.nodelist.append(node)

        element_index = 0
        for iele in range(nelem):
            nodes = [self.nodelist[elements[iele][j+2]] for j in range(nnode)]
            self.ellist.append(Element(element_index,nodes))
            element_index += 1

############ starting element loop #####################
    def solve(self,input_size, model, nbc,bcs,angle,props,tstep,nconv,n_slip,init_ML_step,Fp_list,F_list,K_i,Load_i,garray_list,Cauchy_i,gdisp):

##### Calculate the global load vector(here it is displacement control only) #####
### nbc=total bcs, kbc[node number][direction x=1,y=2,z=3
        nstep = tstep[0] ; dt = tstep[1] ; strain_rate0 = tstep[2]

        neq = 3*len(self.nodelist)

        kbc_n = [0 for i in range(nbc)]
        kbc_d = [0 for i in range(nbc)]
        kbc_v = [0 for i in range(nbc)]

        for ibc in range(nbc):
            kbc_n[ibc] = bcs[ibc][0]-1
            kbc_d[ibc] = bcs[ibc][1]
            kbc_v[ibc] = bcs[ibc][2]
        
        f_bc = [0.0 for i in range(neq)]


        for ibc in range(nbc):
            f_bc[3*kbc_n[ibc]+kbc_d[ibc]-1] = kbc_v[ibc]


        c11 = props[0] ; c12 = props[1] ; c44 = props[2]
        res0_ssd = props[3]


        nbw = 0


##### end of load vector #######
####### initialisation ######
        crystal_type = 'fcc'
        
        XYZ = [[0.0 for i in range(3)] for j in range(8)]
        xyz = [[0.0 for i in range(3)] for j in range(8)]
    

        # garray_tau = []
        vol_g = []
        
        



        neq = 3*len(self.nodelist)
        
        gstiff = [[0.0 for i in range(neq)] for j in range(neq)]
        b_rhs = [0.0 for i in range(neq)]
        gf_nod = [0.0 for i in range(neq)]
        unb_nod = [0.0 for i in range(neq)]
        
        
#--------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------#
#------------------------------- Building Initial Stiffness Matrix --------------------------#
#--------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------#


        for ielem in self.ellist:

            nnode = 8

            lnode = [ielem.nodes[i].inode for i in range(nnode)]

            nb1 = max(lnode)
            nb2 = min(lnode)
            nb = 3*(nb1-nb2+1)
            if nb > nbw:
                nbw = nb

            array_tau = []
            

            qrot = [[0. for ii in range(3)] for jj in range(3)]
            qrot[0][0] = 1.0 ; qrot[1][1] = 1.0 ; qrot[2][2] = 1.0


            angle_el = angle[ielem.index]
            qrot = euler_rotation.euler_rot(angle_el,qrot)

            res_ssd = [res0_ssd for ir in range(n_slip)]
            dgam = [0.0 for ir in range(n_slip)]
            dgam_dta = [0.0 for ir in range(n_slip)]
            S_star0 = [[0.0 for ir in range(3)] for jr in range(3)]

            

            for i in range(nnode):
                XYZ[i][0] = ielem.nodes[i].coorx
                XYZ[i][1] = ielem.nodes[i].coory
                XYZ[i][2] = ielem.nodes[i].coorz

#
                for j in range(3):
                    xyz[i][j] = XYZ[i][j]


#
        ngauss = 8
#
#
##---------------------- Start Calculation for Gauss Points
##
            
######## asembel total stiffness matrix ##############
        gstiff = K_i
    
######## end of asembel total stiffness matrix ##################
#--------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------#
#-------------------------End of Building Initial Stiffness Matrix --------------------------#
#--------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------#
        tnode = len(self.nodelist)
        nelem = len(self.ellist)
        
        strsvol = []
        strnvol = []
        
#        gdisp = [[0.0 for i in range(3)]for j in range(tnode)]

#????????????? it should be modified
 #       time = 0.0
 
        time = init_ML_step*dt

        facload = [dt for i in range(nstep)]
        strain_rate = [strain_rate0 for i in range(nstep)]
        
        
 #===========================================================================================#
 #===========================================================================================#
 #===========================================================================================#
 #===========================================================================================#
 #===========================================================================================#
 #===========================================================================================#
 #===========================================================================================#
 #============ Entering to the time step loop for the actual microstructure =================#
 #===========================================================================================#
 #===========================================================================================#
 #===========================================================================================#
 #===========================================================================================#
 #===========================================================================================#
 #===========================================================================================#
 #===========================================================================================#
 
        for istep in range(init_ML_step,nstep,1):

            for ibc in range(nbc):
                f_bc[3*kbc_n[ibc]+kbc_d[ibc]-1] = kbc_v[ibc]
      

#        factor1 = facload[istep]*strain_rate[istep]
#        dt = facload[istep]; time = time+dt

            factor1 = facload[init_ML_step]*strain_rate[init_ML_step]
            dt = facload[init_ML_step]; time = time+dt
            

#            if istep == init_ML_step:
#                factor1 = Load_i
                
            

#----------------- apply boundary conditions and solve system equations -----------
            b_rhs = solver.solvesyseq(nbc,gstiff,factor1,f_bc,neq,kbc_n,kbc_d,b_rhs,nbw)
#----------------------------------------------------------------------------------
            
            for i in range(tnode):
                for j in range(3):
                    gdisp[i][j] += b_rhs[3*i+j]


#            for i in range(tnode):
#                for j in range(3):
#                    print (i,j,gdisp[i][j])
                

            iconv = 0
            ratio_norm = 1.0e10

            while (iconv <= nconv and abs(ratio_norm) >= 1.0e9):
                if iconv > 0:

                    for i in range(neq):
                        b_rhs[i] = unb_nod[i]

                    f_bc = [0.0 for i in range(neq)]
    #                   b_rhs = solvesmallmatrix(nbc,gstiff,factor1,f_bc,neq,kbc_n,kbc_d,b_rhs)
                    b_rhs = solver.solvesyseq(nbc,gstiff,factor1,f_bc,neq,kbc_n,kbc_d,b_rhs,nbw)


                    for i in range(tnode):
                        for j in range(3):
                            gdisp[i][j] += b_rhs[3*i+j]


                gf_nod = [0.0 for i in range(neq)]
                unb_nod = [0.0 for i in range(neq)]
                gstiff = [[0.0 for i in range(neq)] for j in range(neq)]

                garray_t = garray_list[-1]
                garray_tau = []
                vol_g = []

                ngauss = 8

                Fp_tau_ellist = []
                F_tau_ellist = []
    #----------------------- start of loop element --------------------------
                for ielem in self.ellist:


                    array_tau = []
                    disp = [0.0 for i in range(24)]
                    XYZ = [[0.0 for i in range(3)] for j in range(8)]
                    xyz = [[0.0 for i in range(3)] for j in range(8)]


                    qrot = [[0. for ii in range(3)] for jj in range(3)]
                    qrot[0][0] = 1.0 ; qrot[1][1] = 1.0 ; qrot[2][2] = 1.0


                    angle_el = angle[ielem.index]
                    qrot = euler_rotation.euler_rot(angle_el,qrot)


    #------------ Coordinates in the refrence configuration
    #------------ Assign element-node relationship
                    for i in range(nnode):
                        XYZ[i][0] = ielem.nodes[i].coorx
                        XYZ[i][1] = ielem.nodes[i].coory
                        XYZ[i][2] = ielem.nodes[i].coorz


                        disp[i*3+0] = gdisp[ielem.nodes[i].inode][0]
                        disp[i*3+1] = gdisp[ielem.nodes[i].inode][1]
                        disp[i*3+2] = gdisp[ielem.nodes[i].inode][2]

                        for j in range(3):
                            xyz[i][j] = XYZ[i][j]+disp[i*3+j]

                    schmid = schmid_T.calc_schmid(crystal_type,n_slip,qrot)

                    C_mat = [[[[0. for i in range(3)] for j in range(3)] for k in range(3)] for l in range(3)]

                    C_mat = ciklj.cmat(c11,c12,c44,qrot)

                    strs = [[0.0 for i in range(6)] for j in range(ngauss)]

                    gaussw = gauss.pgauss()


    #---------------------- Start Calculation for Gauss Points

    #------- Get the Deformation Gradient at time = tau for gauss points for the given element

                    Fp_tau_ngauss = []
                    F_tau_ngauss = []
                    for ig in range(ngauss):

    #                        print ig+1

    #???????????????? Transfer deformation gradient to the previous value ?????????

    #---------------------- F_t = F_tau previous -------------


                        # if istep == init_ML_step:
                        #     Ft = F_list[init_ML_step-2][0]
                        # else:
                        Ft = garray_t[ielem.index][ig].defg
                        # print('Ft', Ft)
    #
                        SHP = shp3d.SHP3D(gaussw[0][ig],gaussw[1][ig],gaussw[2][ig],XYZ)

                        F_tau = def_grad_F.calc_F(xyz,SHP[0])
                        
    #
                        # if istep == init_ML_step:
                        #     Fp = Fp_list[init_ML_step-2][0]
                        # else:
                        Fp = garray_t[ielem.index][ig].pdefg

                        # print('Fp', Fp)
                        # print('F_tau', F_tau)
    # #
    # #------------------- Calculate matrix [A]
                        A_mat = calc_ABC.calc_A(F_tau,Fp)
    #
    # #------------------  Calculate matrix [B_alpha]
                        B_alpha = calc_ABC.calc_B(n_slip,A_mat,schmid)
    #
    # #----------------- Calculate the trial elastic stress S_trial
    #                         S_trial = calc_ABC.stress_trial(C_mat,A_mat)
    #
    # #---------------- Calculate the matrix C_alpha for each slip system
                        C_alpha = calc_ABC.calc_C_alpha(n_slip,C_mat,B_alpha)
    #
    # ######################## Enter the User Material ###################
    #                     input_size = 2
                        Fp_list_array = np.array(Fp_list)
                        # print(Fp_list_array.shape)
                        Fp_list_array = Fp_list_array[-input_size:, ielem.index, ig, :, :]
                        # print(Fp_list_array.shape)
                        # print('Fp', Fp_list_array[-1])
                        Fp_list_array = Fp_list_array.reshape(input_size, 9)

                        F_list_array = np.array(F_list)
                        F_list_array = F_list_array[-input_size:, ielem.index, ig, :, :]

                        F_list_array = F_list_array.reshape(input_size, 9)

                        # print('Ft', F_list_array[-1])
                        F_list_array = np.concatenate((F_list_array, np.array(F_tau).reshape(1, 9)))
                        # print(F_list_array.shape)

                        Fp_tau = UMAT_ML(model, Fp_list_array, F_list_array)
                        # print(Fp_tau[1][1])
                        # input()

                        U1, S_star = Cauchy_ML(F_tau, Fp_tau, C_mat)
                        # print(Fp_tau[1][1], U1[1][1], F_tau[1][1])
                        # print(Fp_tau)

                        res_ssd = garray_t[ielem.index][ig].res_ssd
                        res,dgam,dgam_dta = update_constitutive(schmid,res_ssd,props,S_star)

                        U2 = W_mat(Ft,F_tau,Fp_tau,Fp,C_mat,schmid,S_star,dgam,dgam_dta,C_alpha)

                        # U = UMAT(S_trial,C_mat,Fp_t,C_alpha,schmid,F_t,F_tau,res_ssd,dgam,dgam_dta,props,S_star0)
                        #
                        #
                        # U1 = U.itr()
                        # U2 = U.W_mat()

                        strs[ig][0] = U1[0][0] ; strs[ig][1] = U1[1][1] ; strs[ig][2] = U1[2][2]
                        strs[ig][3] = U1[0][1] ; strs[ig][4] = U1[1][2] ; strs[ig][5] = U1[0][2]

    #####################     U1[0]=Fe_tau, U1[1]=Fp_tau, U1[2]=Cauchy, U1[3]=res_ssd, U1[5]=dgam , U1[6]=dgam_dta
                        mg = array(ielem.index,ig,U1,F_tau,Fp_tau,U2, res_ssd,dgam,dgam_dta,S_star)

                        array_tau.append(mg)

                        Fp_tau_ngauss.append(Fp_tau)
                        F_tau_ngauss.append(F_tau)
                        # print(ielem.index, ig, F_tau[1][1])
    ########### End of gauss point loop ##############

    ############### Element reaction forces or internal forces ################

                    InF = element_force.element_re(ngauss,nnode,xyz,strs)

                    lnode = [ielem.nodes[i].inode for i in range(nnode)]

                    gf_nod = asembl_force.asembl_vec(gf_nod,InF,lnode)


    ############## End of Element reaction forces or internal forces ################

                    garray_tau.append(array_tau)

    ############# Calculate the stiffness matrix ###############

                    dep = [garray_tau[ielem.index][ig].dep for ig in range(ngauss)]

                    ev = stiffness_matrix_24.STIF3D(ngauss,nnode,xyz,strs,dep)
                    estiff = ev[0]
                    vol_g.append(ev[1])

    ######## asembel total stiffness matrix ##############
                    for i in range(nnode):
                        for j in range(nnode):
                            for k in range(3):
                                for l in range(3):
                                    gstiff[3*ielem.nodes[i].inode+k][3*ielem.nodes[j].inode+l] += estiff[3*i+k][3*j+l]

                    Fp_tau_ellist.append(Fp_tau_ngauss)
                    F_tau_ellist.append(F_tau_ngauss)
    #----------------- end of element calculations ----------------------
                for i in range(neq):
                    unb_nod[i] = -gf_nod[i]
    ########## Check for L2 norm of unb_nod for convergence ##############

                for ibc in range(nbc):
                    jbc = 3*kbc_n[ibc]+kbc_d[ibc]-1
                    unb_nod[jbc] = 0.0

                norm_gs = math.sqrt(dot_product.dot_product(gf_nod,neq))
                norm_us = math.sqrt(dot_product.dot_product(unb_nod,neq))

                ratio_norm = norm_us/norm_gs

                for i in range(neq):
                    unb_nod[i] = -gf_nod[i]

                print (istep+1,iconv+1,norm_us,ratio_norm)

                iconv += 1

            # print(np.array(Fp_tau_ellist).shape)
            Fp_list.append(Fp_tau_ellist)
            # print(np.array(Fp_list).shape)
            F_list.append(F_tau_ellist)
            # print(np.array(F_list).shape)
            # print(np.array(F_list)[-1, 0, 0, 0, 1])

            garray_list.append(garray_tau)


            strsstrnv = user_print.user_print(tnode,nelem,ngauss,gdisp,garray_tau,vol_g)
            strsv = strsstrnv[0] ; strnv = strsstrnv[1]
            strsvol.append(strsv)
            strnvol.append(strnv)
            
            if istep == nstep-1:
                Cauchy_stress = []
                for iel in range(nelem):
                    Cauchy_stress_el = []
                    for ig in range(ngauss):
                        Cauchy_stress_el.append(garray_tau[iel][ig].cauchy)
                    Cauchy_stress.append(Cauchy_stress_el)

        return strnvol,strsvol,gdisp,Cauchy_stress

    
 #===========================================================================================#
 #===========================================================================================#
 #===========================================================================================#
 #===========================================================================================#
 #===========================================================================================#
 #===========================================================================================#
 #===========================================================================================#
 #=========+++++=== End of the time step loop for the base microstructure ===================#
 #===========================================================================================#
 #===========================================================================================#
 #===========================================================================================#
 #===========================================================================================#
 #===========================================================================================#
 #===========================================================================================#
 #===========================================================================================#
