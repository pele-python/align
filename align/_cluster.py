import numpy as np
from align._transform3d import TransformCluster3D
from align._interfaces import StructuralAlignment
from align._optimize_permutations import optimize_permutations
from align._optimize_rotation import findrotation_kabsch

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
    def __init__(self, can_translate=True, can_rotate=True, can_invert=True, 
                  can_permute=True):
        super(TransformPolicyCluster, self).__init__(can_translate=can_translate,
                can_rotate=can_rotate, can_invert=can_invert, can_permute=can_permute)

class TransformPolicyAtomicCluster(TransformPolicyCluster):
#    def __init__(self):
#        super(TransformPolicyAtomicCluster, self).__init__(can_translate=True,
#                can_rotate=True, can_invert=True, can_permute=True)
        
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
        if self.permlist is None:
            raise RuntimeError("can't optimize permutations if permlist is None")
        return optimize_permutations(X1, X2, self.permlist)
    
    def find_rotation(self, X1, X2):
        return findrotation_kabsch(X1, X2)

class ClusterAlignment(StructuralAlignment):
    def __init__(self):
        self.measure = None
        self.policy = TransformationPolicy()
    
    def get_alignment(self, coords1, coords2):
        assert(self.measure is not None)
        assert(self.policy is not None)
        

        # setup new transformation object to store results
        align = self.policy.new_transformation()
                
        coords1 = coords1.copy()
        coords2 = coords2.copy()
        
        x1 = np.copy(coords1)
        x2 = np.copy(coords2)
    
        com1 = self.measure.get_com(x1)
        self.policy.new_transformation().translate(-com1).apply(x1)
        com2 = self.measure.get_com(x2)
        self.policy.new_transformation().translate(-com2).apply(x2)

        # add com removal to transformation
        align.translate(-com2)
        
        align_best = self.policy.new_transformation()
        self.distbest = self.measure.get_dist(x1, x2)
        
                
        # if we didn't find a perfect match here, try random rotations to optimize the match
        for i in range(self.niter):
            rot = rotations.aa2mx(rotations.random_aa())
            self.check_match(x1, x2, rot, False)
            if(self.policy.can_invert()):
                self.check_match(x1, x2, rot, True)
       
        # add the best transformaton and translate back to com
        return align.transform(align_best).translate(com1)

class AlignClusterSimple(StructuralAlignment):
    """this class will be used in the trivial case where either permutations or rotations are not possible"""
    def __init__(self, transform_policy, measure):
#        self.measure = MeasureCartesian(permlist)
#        self.policy = TransformPolicyAtomicCluster()
        self.measure = measure
        self.policy = transform_policy
        if self.policy.can_permute() and self.policy.can_rotate():
            raise ValueError("ClusterAlignmentSimple cannot be used if both rotations and permutations are possible")
    
    
    def _optimize(self, x1, x2, invert=False):
        if self.policy.can_permute() and self.policy.can_rotate():
            raise RuntimeError("ClusterAlignmentSimple cannot be used if both rotations and permutations are possible")
        x2 = x2.copy()
        tform = self.policy.new_transformation()
        
        if invert:
            tform.invert()
            tform.apply(x2)
        
        if self.policy.can_rotate():
            rot = self.measure.find_rotation(x1, x2)
            tform.rotate(rot)
        
        if self.policy.can_permute():
            perm = self.measure.find_permutation(x1, x2)
            tform.permute(perm)
        
        return tform
            
            
    def get_transformation(self, coords1, coords2):
        assert(self.measure is not None)
        assert(self.policy is not None)
        if self.policy.can_permute() and self.policy.can_rotate():
            raise RuntimeError("ClusterAlignmentSimple cannot be used if both rotations and permutations are possible")
        

        # setup new transformation object to store results
        align = self.policy.new_transformation()
                
        coords1 = coords1.copy()
        coords2 = coords2.copy()
        
        x1 = np.copy(coords1)
        x2 = np.copy(coords2)
    
        trans = lambda : self.transformation.new_transformation()
        # move the center of mass to the origin
        if self.policy.can_translate():
            com1 = self.measure.get_com(x1)
            self.policy.new_transformation().translate(-com1).apply(x1)
            com2 = self.measure.get_com(x2)
            self.policy.new_transformation().translate(-com2).apply(x2)
    
            # add com removal to transformation
            align.translate(-com2)
        
        
        # optimize for rotations OR permutations
        align_noinv = self._optimize(x1, x2, invert=False)
        align_best = align_noinv
        
        if self.policy.can_invert():
            # invert the structure and reoptimize for rotations OR permutations
            align_inv = self._optimize(x1, x2, invert=True)
            
            # compute the distances
            x2_noinv = x2.copy()
            align_noinv.apply(x2_noinv)
            dist_noinv = self.measure.get_dist(x1, x2_noinv)
            
            x2_inv = x2.copy()
            align_inv.apply(x2_inv)
            dist_inv = self.measure.get_dist(x1, x2_inv)
            
            # if the inverted structure gives the better result then save that transformation
            if dist_inv < dist_noinv:
                align_best = align_inv
        
        align.transform(align_best)
            
            
        # add the best transformation and translate back to com
        if self.policy.can_translate():
            align.translate(com1)
        
        return align

    