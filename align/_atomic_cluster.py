"""
This file will hold jake's implementation of the alignment routines for atomic clusters.
I.e. this is a place for jake to experiment.
Once we're happy with the interface we can dump this one or merge them. 
"""
import numpy as np

from _interfaces import ApproximateTransformation

class TransformPolicy(object):
    def __init__(self, can_rotate=False, can_translate=False, 
                 can_invert=False, can_permute=False):
        self._can_rotate = can_rotate
        self._can_translate = can_translate
        self._can_invert = can_invert
        self._can_permute = can_permute
    
    def can_rotate(self):
        return self._can_rotate
    
    def can_translate(self):
        return self._can_translate

    def can_invert(self):
        return self._can_invert

    def can_permute(self):
        return self._can_permute
    
    def rotate(self, coords, rotation):
        raise NotImplementedError

    def translate(self, coords, translation):
        raise NotImplementedError

    def permute(self, coords, permutation):
        raise NotImplementedError

    def invert(self, coords):
        raise NotImplementedError


class TransformPolicyCartesian(TransformPolicy):
    @staticmethod
    def translate(X, d):
        Xtmp = X.reshape([-1,3])
        Xtmp += d #js850> will this work if 
    
    @staticmethod
    def rotate(X, mx,):
        Xtmp = X.reshape([-1,3])
        Xtmp = np.dot(mx, Xtmp.transpose()).transpose()
        X[:] = Xtmp.reshape(X.shape)
    
    @staticmethod        
    def permute(X, perm):
        a = X.reshape(-1,3)[perm].flatten()
        # now modify the passed object, X
        X[:] = a[:]
    
    @staticmethod
    def invert(X):
        X[:] = -X


class MeasureCartesian(object):
    ''' measure rules for atomic clusters '''
    
    def __init__(self, permlist=None):
        self.permlist = permlist
    
    def get_com(self, X):
        X = np.reshape(X, [-1,3])
        natoms = len(X[:,0])
        com = X.sum(0) / natoms
        return com

    def get_dist(self, X1, X2):
        return np.linalg.norm(X1.flatten()-X2.flatten())
    
    def find_permutation(self, X1, X2):
        return find_best_permutation(X1, X2, self.permlist)
    
    def find_rotation(self, X1, X2):
        dist, mx = findrotation(X1, X2)
        return dist, mx
    
class TransformPolicyAtomicCluster(TransformPolicyCartesian):
    def __init__(self):
        super(TransformPolicyAtomicCluster, self).__init__(can_translate=True,
                can_rotate=True, can_invert=True, can_permute=True)




class ApproximateTransformationAtomicCluster(object):



