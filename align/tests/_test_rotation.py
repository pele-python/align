import unittest

import numpy as np

import _utils
from align._atomic_cluster import TransformPolicyAtomicCluster, MeasureCartesian
from align._optimize_rotation import findrotation_kabsch, findrotation_kearsley


class TestFindRotation(unittest.TestCase):
    def setUp(self):
        self.natoms = 11
#        self.x0 = _utils.random_configuration(3*11)
#        self.x1 = _utils.random_configuration(3*11)
#        
        self.policy = TransformPolicyAtomicCluster()
        self.measure = MeasureCartesian()
#    
#        self.rot = findrotation_kabsch(self.x0, self.x1, align_com=True)
    
    def basic_test(self, x0, x1):
        x0bk = x0.copy()
        x1bk = x1.copy()
        
#        rot = findrotation_kabsch(x0, x1)
        rot = findrotation_kearsley(x0, x1)
        
        self.assertTrue((x0 == x0bk).all())
        self.assertTrue((x1 == x1bk).all())
        
        self.policy.rotate(x1, rot)
        self.assertLessEqual(self.measure.get_dist(x0, x1), self.measure.get_dist(x0bk, x1bk))
        
    
    def test1(self):
        """align two random configurations"""
        x0 = _utils.random_configuration(3 * self.natoms)
        x1 = _utils.random_configuration(3 * self.natoms)
        
        # subtract center of mass
        com = self.measure.get_com(x0)
        self.policy.translate(x0, -com)
        com = self.measure.get_com(x1)
        self.policy.translate(x1, -com)

        self.basic_test(x0, x1)
    
    def test2(self):
        """align two configurations that should match exactly
        """
        x0 = _utils.random_configuration(3 * self.natoms)
        com = self.measure.get_com(x0)
        self.policy.translate(x0, -com)
        
        x1 = x0.copy()
        self.policy.rotate(x1, _utils.random_rotation())
        
        self.basic_test(x0, x1)
        
        self.assertLess(np.abs(x0-x1).max(), 1e-3)
        
        
        
        
        
        

if __name__ == "__main__":
    unittest.main()