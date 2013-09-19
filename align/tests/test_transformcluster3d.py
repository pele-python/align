'''
Created on 12 Sep 2013

@author: vr274
'''
import unittest
import numpy as np
from align import TransformCluster3D
from _utils import random_rotation


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
        self.assertLess(np.linalg.norm(coords - coords2), 1e-10)
        
    def testTranslation(self):
        coords = np.array([1., 2., 3., 4., 5., 6.])
        # X
        coords2 = coords.copy() 
        TransformCluster3D().translate(np.array([1., 0, 0])).apply(coords2)
        self.assertLess(np.linalg.norm(coords2 - np.array([2., 2., 3., 5., 5., 6.])), 1e-10)
        # Y
        coords2[:] = coords
        TransformCluster3D().translate(np.array([0., 1, 0])).apply(coords2)
        self.assertLess(np.linalg.norm(coords2 - np.array([1., 3., 3., 4., 6., 6.])), 1e-10)
        # Z
        coords2[:] = coords
        TransformCluster3D().translate(np.array([0., 0, 1])).apply(coords2)
        self.assertLess(np.linalg.norm(coords2 - np.array([1., 2., 4., 4., 5., 7.])), 1e-10)
        # -x-y-z
        coords2[:] = coords
        TransformCluster3D().translate(np.array([-1., -2, -3])).apply(coords2)
        self.assertLess(np.linalg.norm(coords2 - np.array([0., 0., 0., 3., 3., 3.])), 1e-10)
        
    def testPermute(self):
        from _utils import random_permutation
        for i in xrange(100):
            coords2 = np.random.random(3*20)
            coords = coords2.copy().reshape([-1,3])
            coords3 = coords2.copy()
            transform = TransformCluster3D(nsites=20)
            transform2 = TransformCluster3D(nsites=20)
            
            perm = random_permutation(20)
            transform.permute(perm)
            coords = coords[perm]
            transform2.transform(TransformCluster3D(nsites=20).permute(perm))
            
            perm = random_permutation(20)
            transform.permute(perm)
            coords = coords[perm]
            transform2.transform(TransformCluster3D(nsites=20).permute(perm))

            transform.apply(coords2)
            transform2.apply(coords3)
            
            self.assertLess(np.linalg.norm(coords2 - coords.flatten()), 1e-10)
            self.assertLess(np.linalg.norm(coords3 - coords.flatten()), 1e-10)
        
        
    def testCombine(self):
        from _utils import random_rotation
        for i in xrange(100):
            dx1 = np.random.random(3)
            R1 = random_rotation()
            dx2 = np.random.random(3)
            R2 = random_rotation()
            
            transform = TransformCluster3D().translate(dx1) 
            transform.rotate(R1)
            transform.translate(dx2)
            transform.rotate(R2)
            
            x = np.random.random(3*20)
            x2 = x.reshape([-1,3]).copy()            
            x2 = np.dot(R1, (x2+dx1).transpose()).transpose()
            x2 = np.dot(R2, (x2+dx2).transpose()).transpose()
                    
            transform.apply(x)            
            self.assertLess(np.linalg.norm(x - x2.flatten()), 1e-10)            
    
    def testTransform(self):
        # first combine 2 translations
        coords = np.array([1., 2., 3., 4., 5., 6.])
        coords2 = coords.copy()
        
        transform = TransformCluster3D().translate(np.array([1., 2, 3]))
        transform.transform(TransformCluster3D().translate(np.array([2., 3, 4])))
        transform.apply(coords2)
        self.assertLess(np.linalg.norm(coords2 - np.array([4., 7., 10., 7., 10., 13.])), 1e-10)
        
        for i in xrange(100):
            dx1 = np.random.random(3)
            R1 = random_rotation()
            dx2 = np.random.random(3)
            R2 = random_rotation()
                
            transform = TransformCluster3D().translate(dx1) 
            transform.transform(TransformCluster3D().rotate(R1))
            transform.transform(TransformCluster3D().translate(dx2))
            transform.transform(TransformCluster3D().rotate(R2))
            
            x = np.random.random(3*20)
            x2 = x.reshape([-1,3]).copy()            
            x2 = np.dot(R1, (x2+dx1).transpose()).transpose()
            x2 = np.dot(R2, (x2+dx2).transpose()).transpose()
                    
            transform.apply(x)            
            self.assertLess(np.linalg.norm(x - x2.flatten()), 1e-10)            
        
        
    def testRotate(self):        
        for i in xrange(100):
            R = random_rotation()
            transform = TransformCluster3D().rotate(R) 
            x = np.random.random(3*20)
            x2 = x.reshape([-1,3]).copy()
            x2 = np.dot(R, x2.transpose()).transpose()
            transform.apply(x)            
            self.assertLess(np.linalg.norm(x - x2.flatten()), 1e-10)            
        
    
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()