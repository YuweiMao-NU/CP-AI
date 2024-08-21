# Crystal Plasticity Code.

import math,sys


############## Rotation Matrix based on Euler angles #####################

def euler_rot(angle_el,qrot):
    const_pi = math.acos(-1.0)
  
    phi = angle_el[0]*const_pi/180.0
    theta = angle_el[1]*const_pi/180.0
    omega = angle_el[2]*const_pi/180.0
    
    sp = math.sin(phi)
    cp = math.cos(phi)
    st = math.sin(theta)
    ct = math.cos(theta)
    so = math.sin(omega)
    co = math.cos(omega)

    qrot[0][0] = co*cp-so*sp*ct
    qrot[1][0] = co*sp+so*ct*cp
    qrot[2][0] = so*st
    qrot[0][1] = -so*cp-sp*co*ct
    qrot[1][1] = -so*sp+ct*co*cp
    qrot[2][1] = co*st
    qrot[0][2] = sp*st
    qrot[1][2] = -st*cp
    qrot[2][2] = ct


    return qrot

###################### End of Rotation Matrix #####################

