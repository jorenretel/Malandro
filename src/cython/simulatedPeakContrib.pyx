
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
    
    
  def __reduce__(self) :

    return (generalFactory, (type(self),), self.__getstate__())
    
  def __getstate__(self) :
    
    return (self.ccpCode, self.atomName, self.dimNumber, self.firstOrSecondResidue)
  
  def __setstate__(self, state) :

    self.ccpCode, self.atomName, self.dimNumber, self.firstOrSecondResidue = state
