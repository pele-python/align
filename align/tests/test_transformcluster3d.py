'''
Created on 12 Sep 2013

@author: vr274
'''
import unittest
import numpy as np
from align import TransformCluster3D


class TestTransformCluster3D(unittest.TestCase):
    
    def setUp(self):
        pass 

    def tearDown(self):
        pass


    def testReturnTypes(self):
        trans = TransformCluster3D()
        self.assertTrue(isinstance(trans.identity(), TransformCluster3D))
        self.assertTrue(isinstance(trans.translate(np.ones(3)), TransformCluster3D))
        self.assertTrue(isinstance(trans.rotate(np.eye(3)), TransformCluster3D))

    def testInitializationToIdentity(self):
        coords = np.random.random(3*12)
        coords2 = coords.copy()
        TransformCluster3D().apply(coords2)
        self.assertAlmostEqual(coords, coords2)
        
    def testTranslation(self):
        coords = np.array([1., 2., 3., 4., 5., 6.])
        # X
        coords2 = coords.copy() 
        TransformCluster3D().translate(np.array([1., 0, 0])).apply(coords2)
        self.assertAlmostEqual(coords, coords + np.array([2., 2., 3., 5., 5., 6.]))
        # Y
        coords2[:] = coords
        TransformCluster3D().translate(np.array([0., 1, 0])).apply(coords2)
        self.assertAlmostEqual(coords, coords + np.array([1., 3., 3., 4., 6., 6.]))
        # Z
        coords2[:] = coords
        TransformCluster3D().translate(np.array([0., 0, 1])).apply(coords2)
        self.assertAlmostEqual(coords, coords + np.array([1., 2., 4., 4., 5., 7.]))
        # -x-y-z
        coords2[:] = coords
        TransformCluster3D().translate(np.array([-1., -2, -3])).apply(coords2)
        self.assertAlmostEqual(coords, coords + np.array([0., 0., 0., 1., 1., 1.]))
        
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()