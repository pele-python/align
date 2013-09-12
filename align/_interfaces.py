'''
Created on 12 Sep 2013

@author: vr274
'''

class StructuralAlignment(object):
    '''
    base interface for structural alignment algorithms
    '''
        
    def get_transformation(self, coords1, coords2):
        ''' determine the transformation to match coordinates
        
            The function calculates the transformation which transforms
            coords2 to the best match of coords1.
        '''
        raise NotImplementedError
    
    def apply_transformation(self, coords, transformation):
        ''' Apply a transformation to a set of coordinates.
        
            should this act in place or create and return a copy of 
            the coordinates?
        '''
        raise NotImplementedError
    
class ExactStructuralAlignment(StructuralAlignment):
    ''' Interface for structural alignment, where the structures should
        match within numerical errors
    '''
    
    def enumerate(self, coords):
        ''' enumerate symmetry operations
        
            This function enumerates all symmetry operations which do not
            contain an inversion. The inversion just doubles these enumerated
            entries. The point group order of a structure is the 2*length or 
            length of this list for systems with or without inversion symmetry,
            respectively.
            
            TODO: does this function fit here?
        '''  
        raise NotImplementedError
    
    def has_inversion_symmetry(self, coords):
        ''' determine if a specific structure has inversion symmetry 
        
            TODO: do we want this function here?
        '''
        raise NotImplementedError 
   
class ApproximateTransformation(StructuralAlignment):
    ''' Interface for structural alignment, where the structures can be
        different

        ApproximateAlignment tries to find the best transformation which e.g.
        minimizes the rms of 2 structures
    '''
