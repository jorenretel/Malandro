
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
   