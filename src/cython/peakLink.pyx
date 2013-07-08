
cdef class peakLink :
  
  cdef aPeak peak
  
  cdef simulatedPeak simulatedPeak
  
  cdef double score
  
  cdef list resonances
  
  def __init__(self, aPeak peak, simulatedPeak simPeak, list resonances, double score):
    
    self.peak = peak
    self.simulatedPeak = simPeak
    self.score = score
    self.resonances = resonances

  def __reduce__(self) :

    return (generalFactory, (type(self),), self.__getstate__())
    
  def __getstate__(self) :
    
    return (self.peak, self.simulatedPeak, self.resonances, self.score)
  
  def __setstate__(self, state) :

    self.peak, self.simulatedPeak, self.resonances, self.score = state
    
  def getPeak(self):
    
    return self.peak
    
  def getSimulatedPeak(self) :
    
    return self.simulatedPeak
  
  def getResonances(self) :
    
    return self.resonances
  