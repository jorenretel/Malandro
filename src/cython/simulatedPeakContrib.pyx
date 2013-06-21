
cdef class simulatedPeakContrib:
  
  
  cdef str ccpCode
  
  cdef str atomName
  
  cdef str isotopeCode
  
  cdef int dimNumber
  
  cdef int firstOrSecondResidue
  
  cdef object pySimulatedPeakContrib
    
  def __init__(self):

    self.ccpCode = None
    
    self.atomName = None
    
    self.isotopeCode = None
    
  cdef void createPythonStyleObject(self) :
    
    self.pySimulatedPeakContrib = pySimulatedPeakContrib()
    
    self.pySimulatedPeakContrib.ccpCode = self.ccpCode
    
    self.pySimulatedPeakContrib.atomName = self.atomName
    
    self.pySimulatedPeakContrib.dimNumber =  self.dimNumber
    
    self.pySimulatedPeakContrib.firstOrSecondResidue = self.firstOrSecondResidue
