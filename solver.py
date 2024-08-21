# Crystal Plasticity Code.

import math,sys

      

################### apply boundary conditions and solve system equations ##################
def solvesyseq(nbc,gstiff,factor1,f_bc,neq,kbc_n,kbc_d,b_rhs,nbw):

    gstiffm = [[0.0 for i in range(nbw)] for j in range(neq)]

    for ibc in range(nbc):
        jbc = 3*kbc_n[ibc]+kbc_d[ibc]-1
        b_rhs[jbc] = gstiff[jbc][jbc]*1.0e20*factor1*f_bc[jbc]
        gstiff[jbc][jbc] = gstiff[jbc][jbc]*1.0e20

#---------- solve with bandwidth-----------------    
#------ storing total stiffness matrix into neq by nbw matrix ------
#--------- storing first nbw rows---------
    for i in range(neq-nbw+1):
        k = -1
        for j in range(i,i+nbw,1):
            k = k+1
            gstiffm[i][k] = gstiff[i][j]
#--------- storing nbw+1 to neq rows---------
    for i in range(neq-nbw+1,neq):
        for j in range(neq-i):
            gstiffm[i][j] = gstiff[i][j+i]	
#------ End of storing total stiffness matrix into neq by nbw matrix ------
#------- Solve system equations bt Gauss elimination method-------
    for i in range(neq-1):
        nbk = min(neq-i,nbw)
        for j in range(i+1,nbk+i,1):
            j1=j-i
            c = gstiffm[i][j1]/gstiffm[i][0]
            for l in range(j,nbk+i,1):
                l1=l-j
                l2=l-i
                gstiffm[j][l1] += -c*gstiffm[i][l2]
        
            b_rhs[j] += -c*b_rhs[i]

    b_rhs[neq-1] = b_rhs[neq-1]/gstiffm[neq-1][0]

    for ii in range(neq-1):
        i = neq-ii-2
        nbi = min(neq-i,nbw)
        up = 0.0
        for j in range(1,nbi,1):
            up += gstiffm[i][j]*b_rhs[i+j-1+1]
        b_rhs[i] = (b_rhs[i]-up)/gstiffm[i][0]

#------------ End of solveing with bandwidth ---------                    
                

            
    return b_rhs


################### end of apply boundary conditions and solve system equations ##################

