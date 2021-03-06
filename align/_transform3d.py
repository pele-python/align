'''
Created on 12 Sep 2013

@author: vr274
'''

import numpy as np
from align import Transformation

class TransformCluster3D(Transformation):
    ''' base interface for transformations in 3D space
    
    Parameters
    ----------
    nsites : int
        the number of sites in the system.  This must be set if permutations
        are allowed
    '''
    
    def __init__(self, nsites=None):
        self.A = np.eye(4)
        self.permutation = None
        if nsites is not None:
            self.nsites = nsites
        
    def identity(self):
        self.A = np.eye(4)
        return self            
    
    def apply(self, coords):
        c = coords.reshape([-1,3])        
        c[:] = np.dot(self.A[0:3,0:3], c.transpose()).transpose() + self.A[0:3,-1]
        if self.permutation is not None:
            c[:] = c[self.permutation]
        return coords
    
    def transform(self, transformation):
        self.A = np.dot(transformation.A, self.A)
        if transformation.permutation is not None:
            if self.permutation is not None:
                self.permutation = self.permutation[transformation.permutation]
            else:
                self.permutation = transformation.permutation.copy()            
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
    
    def permute(self, permutation):
        if self.permutation is None:
            self.nsites = len(permutation)
        assert(len(self.permutation) == len(permutation))
        self.permutation[:] = self.permutation[permutation]
        return self
    
    def invert(self):
        self.A[:3,] *= -1
        return self
        
    @property
    def nsites(self):
        return len(self.permutation)

    @nsites.setter
    def nsites(self, nsites):
        if not self.permutation is None:
            raise RuntimeError("cannot change number of sites in transformation class")
        self.permutation = np.arange(nsites)
        
    def apply_rotation(self, x, rotation_matrix):
        """apply a rotation matrix to a set of coordinates without modifying self"""
        tform = self.__class__()
        tform.rotate(rotation_matrix).apply(x)
        
    def apply_permutation(self, x, permutation):
        """apply a permutation to a set of coordinates without modifying self"""
        tform = self.__class__()
        tform.permute(permutation).apply(x)
        
    def apply_inversion(self, x):
        """invert a set of coordinates without modifying self"""
        tform = self.__class__()
        tform.invert().apply(x)
    
    def apply_translation(self, x, translation):
        """apply a translation to a set of coordinates without modifying self"""
        tform = self.__class__()
        tform.translate(translation).apply(x)
        
        
        
        
        