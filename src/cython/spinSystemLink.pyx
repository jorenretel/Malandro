
cdef class spinSystemLink :

  cdef int score
  
  cdef list realPeaks, simulatedPeaks
  
  cdef mySpinSystem spinSystem1, spinSystem2
  
  cdef list notFoundSimulatedPeaks
  
  cdef object pySpinSystemLink
  
  cdef list peakLinks
  
  cdef list notFoundPeakLinks
  
  def __init__(self):
    
    self.score = 0
    self.spinSystem1 = None
    self.spinSystem2 = None
    self.simulatedPeaks = []
    self.realPeaks = []
    self.notFoundSimulatedPeaks = []
    self.peakLinks = []
    self.notFoundPeakLinks = []
    
  cdef void createPythonStyleObject(self) :
    
    cdef aPeak peak
    
    cdef simulatedPeak simulatedPeak
    
    self.pySpinSystemLink = pySpinSystemLink()
    
    self.pySpinSystemLink.spinSystem1 = self.spinSystem1.pySpinSystem
    
    self.pySpinSystemLink.spinSystem2 = self.spinSystem2.pySpinSystem
    
    for peak in self.realPeaks :
      
      self.pySpinSystemLink.realPeaks.append(peak.pyPeak)
      
    for simulatedPeak in self.simulatedPeaks :
      
      self.pySpinSystemLink.simulatedPeaks.append(simulatedPeak.pySimulatedPeak)
      
      
  def __reduce__(self) :

    return (generalFactory, (type(self),), self.__getstate__())
    
  def __getstate__(self) :
    
    return (self.peakLinks, self.notFoundPeakLinks)
  
  def __setstate__(self, state) :

    self.peakLinks, self.notFoundPeakLinks = state
    