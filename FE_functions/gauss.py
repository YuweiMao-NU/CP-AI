
import math,sys

########## Loading gauss points for 8 integration points per element ############
def pgauss():
    g = 1.0/math.sqrt(3.0)
    ls = [-1,1,1,-1,-1,1,1,-1]
    lt = [-1,-1,1,1,-1,-1,1,1]
    lz = [-1,-1,-1,-1,1,1,1,1]
    sg = [] ; tg = [] ; zg = [] ; wg = [] 
    for i in range(8):
        s = g*ls[i] ; sg.append(s)
        t = g*lt[i] ; tg.append(t)
        z = g*lz[i] ; zg.append(z)
        w = 1.0     ; wg.append(w)
    return sg,tg,zg,wg
##################### End of loading gauss points 
