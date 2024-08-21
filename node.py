# Crystal Plasticity Code.

import math,sys



class Node:
    def __init__(self,inode,coorx,coory,coorz):
        self.inode = inode
        self.coorx = coorx
        self.coory = coory
        self.coorz = coorz
    def __repr__(self):
        return "(%d,%f,%f,%f)" % (self.inode,self.coorx,self.coory,self.coorz)

