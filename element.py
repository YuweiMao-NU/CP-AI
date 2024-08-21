# Crystal Plasticity Code.

import math,sys






class Element:
    def __init__(self,idx,nodelist=[]):
        self.index = idx
        self.nodes = nodelist[:]

    def __repr__(self):
        return "Element(%d,%s)" % (self.index, self.nodes)


