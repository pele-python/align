'''
Created on 12 Sep 2013

@author: vr274
'''

import numpy as np
from align import Transformation

class TransformCluster3D(Transformation):
    ''' base interface for transformations in 3D space '''
    
    def __init__(self, nsites=None):
        self.A = np.eye(4)
        self.permlist = None
        if nsites is not None:
            self.nsites = nsites
        
    def identity(self):
        self.A = np.eye(4)
        return self            
    
    def apply(self, coords):
        c = coords.reshape([-1,3])        
        c[:] = np.dot(self.A[0:3,0:3], c.transpose()).transpose() + self.A[0:3,-1]
        if self.permlist is not None:
            c[:] = c[self.permlist]
        return coords
    
    def transform(self, transformation):
        self.A = np.dot(transformation.A, self.A)
        if transformation.permlist is not None:
            if self.permlist is not None:
                self.permlist = self.permlist[transformation.permlist]
            else:
                self.permlist = transformation.permlist.copy()            
        return self
    
    def translate(self, u):
        ''' apply a translation 
        
            Parameters
            ----------
            u : np.array 3
                translation vector
            
            Returns
            -------
            self
        '''
        self.A[0:3,-1]+=u
        return self
        
    def rotate(self, R):
        ''' apply a rotation 
        
            Parameters
            ----------
            R : np.array 3x3
                rotation matrix
            
            Returns
            -------
            self
        '''
        self.A[0:3,0:3] = np.dot(R, self.A[0:3,0:3])
        self.A[0:3,-1]=np.dot(R, self.A[0:3,-1])
        return self
    
    def permute(self, permutations):
        assert(self.permlist is not None)
        assert(len(self.permlist) == len(permutations))
        self.permlist[:] = self.permlist[permutations]
        return self
        
    @property
    def nsites(self):
        return len(self.permlist)

    @nsites.setter
    def nsites(self, nsites):
        if not self.permlist is None:
            raise RuntimeError("cannot change number of sites in transformation class")
        self.permlist = np.arange(nsites)