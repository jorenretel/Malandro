
cdef class spinSystemLink :

  cdef int score
  
  cdef list realPeaks, simulatedPeaks
  
  cdef mySpinSystem spinSystem1, spinSystem2
  
  cdef list notFoundSimulatedPeaks
  
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
 
  def __reduce__(self) :

    return (generalFactory, (type(self),), self.__getstate__())
    
  def __getstate__(self) :
    
    return (self.peakLinks)#, self.notFoundPeakLinks)
  
  def __setstate__(self, state) :

    self.peakLinks = state #self.notFoundPeakLinks
    
  def getContributingResonances(self) :
    
    resonances = []
    
    for pl in self.peakLinks :
      
      resonances.extend(pl.getResonances())
      
    return set(resonances)  
  
  def getAllResonances(self) :
  
    resonances = []
    
    for pl in self.peakLinks :
      
      resonances.extend(pl.getResonances())
      
    for pl in self.notFoundPeakLinks :
      
      resonances.extend(pl.getResonances())
      
    return set(resonances)  
    
  def getSimulatedPeaks(self):        #TODO: update, simulatedPeaks is not used anymore
        
    return self.simulatedPeaks
  
  def getRealPeaks(self):
    
    return self.realPeaks
  
  def getPeakLinks(self) :
    
    return self.peakLinks
  
  def getNotFoundPeakLinks(self) :
    
    return self.notFoundPeakLinks
  
  def getAllPeakLinks(self) :
    
    return self.peakLinks + self.notFoundPeakLinks