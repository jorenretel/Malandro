
cdef class simulatedPeak :
  
  cdef double colabelling
  
  cdef list simulatedPeakContribs
  
  cdef aSpectrum spectrum
  
  cdef object pySimulatedPeak
  
  
    
  def __init__(self):

    self.colabelling = 0.0
    
    self.simulatedPeakContribs = []
    
    self.spectrum = None
    
  cdef void createPythonStyleObject(self) :
    
    cdef simulatedPeakContrib contrib
    
    self.pySimulatedPeak = pySimulatedPeak()
    
    self.pySimulatedPeak.colabelling = self.colabelling
    
    for contrib in self.simulatedPeakContribs :
      
      self.pySimulatedPeak.simulatedPeakContribs.append(contrib.pySimulatedPeakContrib)
   
  def __reduce__(self) :

    return (generalFactory, (type(self),), self.__getstate__())
    
  def __getstate__(self) :
    
    return (self.spectrum, self.simulatedPeakContribs, self.colabelling)
  
  def __setstate__(self, state) :

    self.spectrum, self.simulatedPeakContribs, self.colabelling = state