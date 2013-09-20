import unittest
import numpy as np
import itertools

import align._utils as utils
from align._cluster import TransformPolicyAtomicCluster
from align._cluster import AlignClusterSimple, MeasureCartesian

class TestAlignClusterSimple(unittest.TestCase):
    def setUp(self):
        policy = TransformPolicyAtomicCluster(can_permute=True, can_rotate=False, can_invert=True, can_translate=True)
        self.setUp1(policy)

    def setUp1(self, policy):
        self.policy = policy
        self.natoms = 11
        
        if self.policy.can_permute():
            self.permlist = [range(self.natoms)]
        else:
            self.permlist = None
        
        self.measure = MeasureCartesian(self.permlist)
        self.align = AlignClusterSimple(self.policy, self.measure)
        
        self.transform_list = []
        if self.policy.can_invert():
            self.transform_list.append(utils.invert)
        if self.policy.can_translate():
            self.transform_list.append(utils.translate_randomly)
        if self.policy.can_rotate():
            self.transform_list.append(utils.rotate_randomly)
        if self.policy.can_permute():
            self.transform_list.append(self.permute_randomly)
            
    
    def permute_randomly(self, x):
        utils.permute_randomly(x, self.permlist, self.natoms)
    
    def basic_test(self, x0, x1):
        """find the best alignment and do some basic tests"""
        x0bk = x0.copy()
        x1bk = x1.copy()
        
        tform = self.align.get_transformation(x0, x1)
        
        # test the passed structures are not altered
        self.assertTrue((x0 == x0bk).all())
        self.assertTrue((x1 == x1bk).all())
        
        # apply the transformation to x1 and test that distance is not worse
        tform.apply(x1)
        self.assertLessEqual(self.measure.get_dist(x0, x1), self.measure.get_dist(x0bk, x1bk) + 1e-10)

    def test1(self):
        """align two random configurations"""
        x0 = utils.random_configuration(3 * self.natoms)
        x1 = utils.random_configuration(3 * self.natoms)
        
        self.basic_test(x0, x1)
    

    def fchain(self, flist):
        """return a function which applies all of the functions in flist to the input"""
        def function_chain(x):
            for f in reversed(flist):
#                print f.__name__
                f(x)
        return function_chain

    def exact_match_test(self, f):
        x0 = utils.random_configuration(3 * self.natoms)
        x1 = x0.copy()

        # apply transformations to x1
        f(x1) 

        self.basic_test(x0, x1)
        self.assertLess(self.measure.get_dist(x0, x1), 1e-3)
    
    def test_1(self):
        """run exact_match_test() for all transformations"""
        for f in self.transform_list:
            self.exact_match_test(f)
            
    def test_2(self):
        """run exact_match_test() on all combinations of length 2 of the transformations
        """
        for flist in itertools.product(self.transform_list, repeat=2):
            tform = self.fchain(flist)
            self.exact_match_test(tform)

    def test_3(self):
        """run exact_match_test() on all combinations of length 2 of the transformations
        """
        for flist in itertools.product(self.transform_list, repeat=2):
            tform = self.fchain(flist)
            self.exact_match_test(tform)

    def test_4(self):
        """run exact_match_test() on all combinations of length 2 of the transformations
        """
        for flist in itertools.product(self.transform_list, repeat=4):
            tform = self.fchain(flist)
            self.exact_match_test(tform)

    def test_10(self):
        """run do_tests() on combinations of length 10 of the transformations
        """
        maxiter = 1000
        i = 0
        for flist in itertools.product(self.transform_list, repeat=10):
            tform = self.fchain(flist)
            self.exact_match_test(tform)
            i += 1
            if i > maxiter:
                break


class TestAlignClusterSimple1(TestAlignClusterSimple):
    def setUp(self):
        policy = TransformPolicyAtomicCluster(can_permute=False, can_rotate=True, can_invert=True, can_translate=True)
        self.setUp1(policy)

class TestAlignClusterSimple2(TestAlignClusterSimple):
    def setUp(self):
        policy = TransformPolicyAtomicCluster(can_permute=False, can_rotate=True, can_invert=False, can_translate=True)
        self.setUp1(policy)

class TestAlignClusterSimple3(TestAlignClusterSimple):
    def setUp(self):
        policy = TransformPolicyAtomicCluster(can_permute=False, can_rotate=True, can_invert=True, can_translate=False)
        self.setUp1(policy)

class TestAlignClusterSimple4(TestAlignClusterSimple):
    def setUp(self):
        policy = TransformPolicyAtomicCluster(can_permute=True, can_rotate=False, can_invert=False, can_translate=False)
        self.setUp1(policy)


if __name__ == "__main__":
    unittest.main()     
        