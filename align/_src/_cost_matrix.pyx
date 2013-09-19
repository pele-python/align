import numpy as np
cimport numpy as np
cimport cython


@cython.boundscheck(False)
@cython.wraparound(False)
def _cost_matrix_cartesian(np.ndarray[double, ndim=2, mode="c"] x1, 
                           np.ndarray[double, ndim=2, mode="c"] x2
                           ):
    """return a matrix of the squared cartesian distances between the atoms"""
    cdef unsigned int natoms, d, i, j, k
    cdef double r2, r

    assert x1.shape[0] == x2.shape[0]
    assert x1.shape[1] == x2.shape[1]
    natoms, d = x1.shape[0], x1.shape[1]

#    np.ndarray[double, ndim=1, mode="c"] dr = np.zeros(d)
#    np.ndarray[double, ndim=1, mode="c"] iboxvec = 1. / boxvec
#    cdef int periodic = 0
    
    cdef np.ndarray[double, ndim=2, mode="c"] cost_matrix = np.zeros([natoms, natoms])

    
    
    for i in xrange(natoms):
        for j in xrange(natoms):
            r2 = 0.
            for k in xrange(d):
                r = x1[i,k] - x2[j,k]
                r2 += r * r
            cost_matrix[j,i] = r2
            
    return cost_matrix
