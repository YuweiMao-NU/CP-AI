# Crystal Plasticity Code.

import math,sys



      


################ Calculation of logarithmic strain #############
def true_strain(a):
    v = [[0.0 for i in range(3)] for j in range(3)]
    beta = [0.0 for i in range(3)]
    
    p1 = a[0][1]**2 + a[0][2]**2 + a[1][2]**2

    if p1 < 1.0e-10:
        for ivec in range(3):
            beta[ivec] = a[ivec][ivec]
            v[ivec][ivec] = 1.0
    elif p1 >= 1.0e-10:

        q = 0.0
        for i in range(3):
            q += a[i][i]
        q = q/3.0
        Iden = [[0.0 for i in range(3)] for j in range(3)]
        Iden[0][0] = Iden[1][1] = Iden[2][2] = 1.0
        tmp1 = [[0.0 for i in range(3)] for j in range(3)]
        tmp2 = [[0.0 for i in range(3)] for j in range(3)]

        for i in range(3):
            for j in range(3):
                tmp1[i][j] = a[i][j]-q*Iden[i][j]

        for i in range(3):
            for j in range(3):
                for k in range(3):
                    tmp2[i][j] += tmp1[i][k]*tmp1[k][j]
                tmp2[i][j] = tmp2[i][j]/6.0
        p = 0.0
        for i in range(3):
            p += tmp2[i][i]
        p = math.sqrt(p)

        b = [[0.0 for i in range(3)] for j in range(3)]
        for i in range(3):
            for j in range(3):
                b[i][j] = tmp1[i][j]/p
       
        det_b = b[0][0]*b[1][1]*b[2][2]-b[0][0]*b[1][2]*b[2][1]-b[1][0]*b[0][1]\
                  *b[2][2]+b[1][0]*b[0][2]*b[2][1]+b[2][0]*b[0][1]*b[1][2]-b[2][0]\
                  *b[0][2]*b[1][1]
        r = det_b/2.0
        if r <= -1.0:
            phi = math.pi/3.0
        elif r >= 1.0:
            phi = 0.0
        elif r > -1.0 and r < 1.0:
           phi = math.acos(r)/3.0
       

        beta[0] = q + 2.0 * p * math.cos(phi)
        beta[1] = q + 2.0 * p * math.cos(phi + (2.0*math.pi/3.0))
        beta[2] = 3.0 * q - beta[0] - beta[1]

        tmp3 = [[0.0 for i in range(3)] for j in range(3)]

        for ivec in range(3):
            for i in range(3):
                for j in range(3):
                    tmp3[i][j] = a[i][j]-beta[ivec]*Iden[i][j]

            if tmp3[2][2] != 0.0:
                alpha1 = 1.0
                alpha21 = -tmp3[1][0]+(tmp3[1][2]*tmp3[2][0]/tmp3[2][2])
                alpha22 = tmp3[1][1]-(tmp3[1][2]*tmp3[2][1]/tmp3[2][2])
                alpha2 = alpha21/alpha22
                alpha3 = (-tmp3[2][0]/tmp3[2][2])-(tmp3[2][1]/tmp3[2][2])*alpha2
            if tmp3[2][2] == 0.0 and tmp3[2][1] != 0.0:
                alpha1 = 1.0
                alpha2 = -tmp3[2][0]/tmp3[2][1]
                if tmp3[1][2] != 0.0:
                    alpha3 = (tmp3[1][0]-tmp3[1][1]*alpha2)/tmp3[1][2]
                elif tmp3[1][2] == 0.0 and tmp3[0][2] != 0.0:
                    alpha3 = (tmp3[0][0]-tmp3[0][1]*alpha2)/tmp3[1][2]
                elif tmp3[1][2] == 0.0 and tmp3[0][2] == 0.0:
                    alpha3 = 0.0

            if tmp3[2][2] == 0.0 and tmp3[2][1] == 0.0 and tmp3[2][1] != 0.0:
                alpha1 = 0.0
                if tmp3[1][2] != 0.0:
                    alpha2 = 1.0
                    alpha3 = -tmp3[1][1]/tmp3[1][2]
                elif tmp3[1][2] == 0.0 and tmp3[0][2] != 0.0:
                    alpha2 = 0.0
                    alpha3 = 1.0

            alpha = math.sqrt(alpha1**2+alpha2**2+alpha3**2)
            alpha1 = alpha1/alpha ; alpha2 = alpha2/alpha ; alpha3 = alpha3/alpha
            v[ivec][0] = alpha1 ; v[ivec][1] = alpha2 ; v[ivec][2] = alpha3

    strain = [[0.0 for i in range(3)] for j in range(3)]

    tmpv1 = [[0.0 for i in range(3)] for j in range(3)]
    tmpv2 = [[0.0 for i in range(3)] for j in range(3)]
    tmpv3 = [[0.0 for i in range(3)] for j in range(3)]

    for i in range(3):
        for j in range(3):
            tmpv1[i][j] = v[0][i]*v[0][j]
            
    for i in range(3):
        for j in range(3):
            tmpv2[i][j] = v[1][i]*v[1][j]

    for i in range(3):
        for j in range(3):
            tmpv3[i][j] = v[2][i]*v[2][j]
    for i in range(3):
        for j in range(3):
            strain[i][j] = 0.0
            strain[i][j] = 0.5*math.log(beta[0])*tmpv1[i][j]+\
                           0.5*math.log(beta[1])*tmpv2[i][j]+\
                           0.5*math.log(beta[2])*tmpv3[i][j]

    return strain

################ end of calculation of logarithmic strain ######


