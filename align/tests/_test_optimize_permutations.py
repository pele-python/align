import unittest

import numpy as np

import _utils
from align._src._cost_matrix import _cost_matrix_cartesian
from align._optimize_permutations import optimize_permutations, _optimize_permutations_hungarian, _optimize_permutations_optim
from align._atomic_cluster import TransformPolicyAtomicCluster, MeasureCartesian

def _make_cost_matrix(X1, X2):
    """
    return the cost matrix for use in the hungarian algorithm.
    
    the cost matrix is the distance matrix (squared) for all atoms in atomlist
    """
    cost = (((X1[np.newaxis,:] - X2[:,np.newaxis,:])**2).sum(2))
    return cost

class TestCostMatrix(unittest.TestCase):
    def setUp(self):
        self.natoms = 4
        self.x1 = _utils.random_configuration(self.natoms * 3).reshape(-1, 3)
        self.x2 = _utils.random_configuration(self.natoms * 3).reshape(-1, 3)
    
    def test1(self):
        m1 = _make_cost_matrix(self.x1, self.x2)
        m2 = _cost_matrix_cartesian(self.x1, self.x2)
        print m1
        print m2
        self.assertLess(np.abs(m1 - m2).max(), 1e-3)


class TestOptimizePermutations(unittest.TestCase):
    def setUp(self):
        self.setUp1()

    def setUp1(self):
        self.natoms = 11
        self.policy = TransformPolicyAtomicCluster()
        self.permlist = [range(self.natoms)]
        self.measure = MeasureCartesian(permlist=self.permlist)
        self.wrapped_optimize_permutations = lambda x1, x2: optimize_permutations(x1, x2, self.permlist) 

    def basic_test(self, x0, x1):
        x0bk = x0.copy()
        x1bk = x1.copy()
        
        perm = self.wrapped_optimize_permutations(x0, x1)
        
        self.assertTrue((x0 == x0bk).all())
        self.assertTrue((x1 == x1bk).all())
        
        self.policy.permute(x1, perm)
        self.assertLessEqual(self.measure.get_dist(x0, x1), self.measure.get_dist(x0bk, x1bk))


    def test1(self):
        x0 = _utils.random_configuration(self.natoms * 3)
        x1 = _utils.random_configuration(self.natoms * 3)
        
        self.basic_test(x0, x1)
        
        
    def test2(self):
        """align two configurations that should match exactly
        """
        x0 = _utils.random_configuration(3 * self.natoms)
        x1 = x0.copy()
        self.policy.permute(x1, _utils.random_permutation_permlist(self.permlist, self.natoms))
        
        self.basic_test(x0, x1)
        
        self.assertLess(np.abs(x0-x1).max(), 1e-3)
        
        
class TestOptimizePermutationsCostMatrix(TestOptimizePermutations):
    def setUp(self):
        self.setUp1()
        self.wrapped_optimize_permutations = lambda x1, x2: optimize_permutations(x1, x2, self.permlist, user_cost_matrix=_cost_matrix_cartesian)

class TestOptimizePermutationsHungarian(TestOptimizePermutations):
    def setUp(self):
        self.setUp1()
        self.wrapped_optimize_permutations = lambda x1, x2: _optimize_permutations_hungarian(x1, x2, make_cost_matrix=_cost_matrix_cartesian)

class TestOptimizePermutationsOPTIM(TestOptimizePermutations):
    def setUp(self):
        self.setUp1()
        self.wrapped_optimize_permutations = lambda x1, x2: _optimize_permutations_optim(x1, x2)

class TestOptimizePermutationsBinary(TestOptimizePermutations):
    def setUp(self):
        self.natoms = 31
        self.ntypeA = int(self.natoms * 0.8)
        self.policy = TransformPolicyAtomicCluster()
        self.permlist = [range(self.ntypeA), range(self.ntypeA, self.natoms)]
        self.measure = MeasureCartesian(permlist=self.permlist)
        self.wrapped_optimize_permutations = lambda x1, x2: optimize_permutations(x1, x2, self.permlist)


  
if __name__ == "__main__":
    unittest.main()