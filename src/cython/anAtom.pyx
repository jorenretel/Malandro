

cdef class anAtom :

  cdef object ccpnAtom

  cdef str atomName

  cdef aResidue residue

  cdef list assignmentPossibilityDimensions
  
  cdef dict labelInfoTemp 

  def __init__(self, residue, ccpnAtom):
    
    self.residue = residue
    
    self.ccpnAtom = ccpnAtom
    
    self.atomName = ccpnAtom.chemAtom.name

    self.assignmentPossibilityDimensions = []
    
    self.labelInfoTemp = {}
