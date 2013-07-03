
cdef class peakLink :
  
  cdef aPeak peak
  
  cdef simulatedPeak simPeak
  
  cdef double score
  
  cdef list resonances
  
  def __init__(self, aPeak peak, simulatedPeak simPeak, double score):
    
    self.peak = peak
    self.simPeak = simPeak
    self.score = score
    self.resonances = []

  def __reduce__(self) :

    return (generalFactory, (type(self),), self.__getstate__())
    
  def __getstate__(self) :
    
    return (self.peak, self.simPeak, self.score, self.resonances)
  
  def __setstate__(self, state) :

    self.peak, self.simPeak, self.score, self.resonances = state