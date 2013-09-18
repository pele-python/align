import numpy as np


def random_q():
    """
    uniform random rotation in angle axis formulation
    input: 3 uniformly distributed random numbers
    uses the algorithm given in
    K. Shoemake, Uniform random rotations, Graphics Gems III, pages 124-132. Academic, New York, 1992.
    This first generates a random rotation in quaternion representation. We should substitute this by
    a direct angle axis generation, but be careful: the angle of rotation in angle axis representation
    is NOT uniformly distributed
    """
    from numpy import sqrt, sin, cos, pi
    u = np.random.uniform(0,1,[3])
    q = np.zeros(4, np.float64)
    q[0] = sqrt(1.-u[0]) * sin(2.*pi*u[1])
    q[1] = sqrt(1.-u[0]) * cos(2.*pi*u[1])
    q[2] = sqrt(u[0]) * sin(2.*pi*u[2])
    q[3] = sqrt(u[0]) * cos(2.*pi*u[2])
    return q

def q2mx( qin ):
    """quaternion to rotation matrix"""
    Q = qin / np.linalg.norm(qin)
    RMX = np.zeros([3,3], np.float64)
    Q2Q3 = Q[1]*Q[2];
    Q1Q4 = Q[0]*Q[3];
    Q2Q4 = Q[1]*Q[3];
    Q1Q3 = Q[0]*Q[2];
    Q3Q4 = Q[2]*Q[3];
    Q1Q2 = Q[0]*Q[1];

    RMX[0,0] = 2.*(0.5 - Q[2]*Q[2] - Q[3]*Q[3]);
    RMX[1,1] = 2.*(0.5 - Q[1]*Q[1] - Q[3]*Q[3]);
    RMX[2,2] = 2.*(0.5 - Q[1]*Q[1] - Q[2]*Q[2]);
    RMX[0,1] = 2.*(Q2Q3 - Q1Q4);
    RMX[1,0] = 2.*(Q2Q3 + Q1Q4);
    RMX[0,2] = 2.*(Q2Q4 + Q1Q3);
    RMX[2,0] = 2.*(Q2Q4 - Q1Q3);
    RMX[1,2] = 2.*(Q3Q4 - Q1Q2);
    RMX[2,1] = 2.*(Q3Q4 + Q1Q2);
    return RMX

def random_mx():
    return q2mx(random_q())

def sum_internal_distances(x):
    """return the sqrt of the sum of the squared internal distances
    
    this can be used to check that the relative positions of the atoms are unchanged 
    """
    s = ((x[:,np.newaxis] - x[np.newaxis,:])**2).sum()
    return np.sqrt(s)

def random_configuration(n):
    return np.random.uniform(-1, 1, n)

def random_translation():
    return np.random.uniform(-1,1,3)

def random_rotation():
    return q2mx(random_q())

def random_permutation(natoms):
    perm = range(natoms)
    np.random.shuffle(perm)
    return perm
