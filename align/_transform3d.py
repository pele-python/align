'''
Created on 12 Sep 2013

@author: vr274
'''

import numpy as np
from align import Transformation

class TransformCluster3D(Transformation):
    ''' base interface for transformations in 3D space '''
    
    def __init__(self):
        self.A = np.eye(4)            
    
    def transform(self, transformation):
        self.A = np.dot(transformation.A, self.A)
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
        raise NotImplementedError
        
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
        return self