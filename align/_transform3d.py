'''
Created on 12 Sep 2013

@author: vr274
'''

class Transform3D(object):
    ''' base interface for transformations in 3D space '''
        
    
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
