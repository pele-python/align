import numpy as np
from align._transform3d import TransformCluster3D
from align._interfaces import StructuralAlignment
from align._optimize_permutations import optimize_permutations

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

    def new_transformation(self):
        raise NotImplementedError

class TransformPolicyCluster(TransformPolicy):
    def __init__(self):
        super(TransformPolicyCluster, self).__init__(can_translate=True,
                can_rotate=True, can_invert=True, can_permute=True)

class TransformPolicyAtomicCluster(TransformPolicyCluster):
    def __init__(self):
        super(TransformPolicyAtomicCluster, self).__init__(can_translate=True,
                can_rotate=True, can_invert=True, can_permute=True)
        
    def new_transformation(self):
        return TransformCluster3D()

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
        return optimize_permutations(X1, X2, self.permlist)
    
    def find_rotation(self, X1, X2):
        return find_rotation_kabsch(X1, X2)

class ClusterAlignment(StructuralAlignment):
    def ClusterAlignment(self):
        self.measuer = None
        self.transform = TransformationPolicy()
    
    def get_alignment(self, coords1, coords2):
        assert(self.measure is not None)
        assert(self.transform is not None)
        

        # setup new transformation object to store results
        align = self.transform.new_transformation()
                
        coords1 = coords1.copy()
        coords2 = coords2.copy()
        
        x1 = np.copy(coords1)
        x2 = np.copy(coords2)
    
        com1 = self.measure.get_com(x1)
        self.transform.new_transformation.translate(-com1).apply(x1)
        com2 = self.measure.get_com(x2)
        self.transform.new_transformation.translate(-com1).apply(x2)

        # add com removal to transformation
        align.translate(-com1)
        
        align_best = self.transform.new_transformation()
        self.distbest = self.measure.get_dist(x1, x2)
        
                
        # if we didn't find a perfect match here, try random rotations to optimize the match
        for i in range(self.niter):
            rot = rotations.aa2mx(rotations.random_aa())
            self.check_match(x1, x2, rot, False)
            if(self.transform.can_invert()):
                self.check_match(x1, x2, rot, True)
       
        # add the best transformaton and translate back to com
        return align.transform(align_best).translate(com1)

    