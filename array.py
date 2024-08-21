# Crystal Plasticity Code.

import math,sys

class array:
    def __init__(self,ielem,igauss,cauchy,defg,pdefg,dep,res_ssd,dgam,dgam_dta,S_star):
        self.ielem = ielem
        self.igauss = igauss
        self.cauchy = cauchy

        self.defg = defg
        self.pdefg = pdefg
        self.dep = dep
        self.res_ssd = res_ssd
        self.dgam = dgam
        self.dgam_dta = dgam_dta
        self.S_star = S_star

    def __repr__(self):
        return "(%d,%d,%s,%s,%s,%s,%s)" % (self.ielem,self.igauss,self.cauchy,self.defg,self.pdefg,self.dep,self.res_ssd,self.dgam,self.dgam_dta,self.S_star)

   

