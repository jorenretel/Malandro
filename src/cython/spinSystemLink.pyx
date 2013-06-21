
cdef class spinSystemLink :

  cdef int score
  
  cdef list realPeaks
  
  cdef mySpinSystem spinSystem1
  
  cdef mySpinSystem spinSystem2
  
  cdef list simulatedPeaks
  
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
    