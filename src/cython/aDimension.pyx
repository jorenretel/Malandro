
cdef class aDimension :
  
  cdef object ccpnDim, pyDimension
  
  cdef double ppmValue, tolerance
  
  cdef int dimNumber
  
  cdef list possibleContributions, nonLabelledResonances
  
  cdef aPeak peak
  
  
  def __init__(self, peak, ccpnDim):

    self.peak = peak
    self.ccpnDim = ccpnDim
    self.ppmValue = ccpnDim.value
    self.dimNumber = ccpnDim.dataDim.expDim.refExpDim.dim
    self.tolerance = getAnalysisDataDim(ccpnDim.dataDim).assignTolerance
      
    self.possibleContributions = []                         # All resonances in the resonanceList that could potentially contribute to this dimension of the peak 
    self.nonLabelledResonances = []                     # Here all resonances are gathered that can not contribute to the peak because of the labelling scheme. They are collected anywya to search for peaks that explicitely should NOT be there.

    self.peak = None
    
    
  cdef void createPythonStyleObject(self) :
    
    self.pyDimension = pyDimension()
          
    self.pyDimension.ppmValue = self.ppmValue
    
    self.pyDimension.dimNumber = self.dimNumber
    
    
  def __reduce__(self) :

    return (generalFactory, (type(self),), self.__getstate__())
    
  def __getstate__(self) :
    
    return (self.dimNumber, self.ppmValue, self.peak)
  
  def __setstate__(self, state) :

    self.dimNumber, self.ppmValue, self.peak = state
    
  def getChemicalShift(self) :
    
    return self.ppmValue
    
  def getDimNumber(self) :
    
    return self.dimNumber
  
