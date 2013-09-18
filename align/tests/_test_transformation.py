import unittest
import itertools
import numpy as np

import _utils
from align._atomic_cluster import TransformPolicyAtomicCluster

class TestTransformation(unittest.TestCase):
    def setUP(self):
        self.natoms = 11
        self.policy = TransformPolicyAtomicCluster()
        
        self.transform_list = [self.invert, self.rrot, self.rtrans, self.rperm]
        
    def rrot(self, tform, coords):
        """apply random rotation to both coords and transformation object"""
        mx = _utils.random_mx()
        self.policy.rotate(coords, mx)
        tform.apply_rotation(mx)
    
    def rtrans(self, tform, coords):
        """apply random translation to both coords and transformation object"""
        xyz = _utils.random_translation()
        self.policy.translate(coords, xyz)
        tform.apply_translation(xyz)
    
    def invert(self, tform, coords):
        """apply inversion to both coords and transformation object"""
        tform.apply_inversion()
        self.policy.invert(coords)
    
    def rperm(self, tform, coords):
        """apply random permutation to both coords and transformation object"""
        perm = _utils.random_permutation(len(coords)/3)
        tform.apply_permutation(perm)
    
    def do_test(self, f):
        """test that that a transofrmation object can reproduce the transformation in f"""
        x = _utils.random_configuration(self.natoms * 3)
        tform = Transformation()
        
        xbkup = x.copy()
        x2 = x.copy()
        
        # do transormation
        f(tform, x)
        
        tform.apply(x2)
        
        #x2 and x should be the same
        self.assertLess(np.abs(x-x2).max(), 1e-3)
        
        
    
    def fchain(self, flist):
        """return a function which applies all of the functions in flist to the input"""
        def function_chain(tform, x):
            for f in reversed(flist):
                f(tform, x)
        return function_chain

    def test_1(self):
        """run do_tests() on transformations
        """
        for f in self.transform_list:
            self.apply_tests(f)


    def test_2(self):
        """run do_tests() on all combinations of length 2 of the transformations
        """
        for flist in itertools.product(self.transform_list, repeat=2):
            tform = self.fchain(flist)
            self.apply_tests(tform)


    def test_3(self):
        """run do_tests() on all combinations of length 2 of the transformations
        """
        for flist in itertools.product(self.transform_list, repeat=3):
            tform = self.fchain(flist)
            self.apply_tests(tform)

        
        