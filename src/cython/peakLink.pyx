
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
