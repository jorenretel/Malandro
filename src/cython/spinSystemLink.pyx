
cdef class SpinSystemLink :

  cdef int score
  
  cdef list realPeaks, simulatedPeaks
  
  cdef Residue residue1, residue2
  
  cdef SpinSystem spinSystem1, spinSystem2
  
  cdef list notFoundSimulatedPeaks
  
  cdef list peakLinks
  
  cdef list notFoundPeakLinks
  
  def __init__(self, residue1=None,residue2=None,spinSystem1=None,spinSystem2=None):
    
    self.score = 0
    self.residue1 = residue1
    self.residue2 = residue2
    self.spinSystem1 = spinSystem1
    self.spinSystem2 = spinSystem2
    self.simulatedPeaks = []
    self.realPeaks = []
    self.notFoundSimulatedPeaks = []
    self.peakLinks = []
    self.notFoundPeakLinks = []

  cdef void addPeak(self,Peak realPeak, SimulatedPeak simPeak, list resonances, double score) :
    
    if realPeak is not None :
      
      newPeakLink = PeakLink(realPeak,simPeak,resonances, score) 
      self.peakLinks.append(newPeakLink)
      
    else :
      
      self.notFoundSimulatedPeaks.append(simPeak)
 
  def __reduce__(self) :

    return (generalFactory, (type(self),), self.__getstate__())
    
  def __getstate__(self) :
    
    return (self.peakLinks, self.notFoundSimulatedPeaks, self.spinSystem1, self.spinSystem2) #, self.notFoundPeakLinks)#, self.notFoundPeakLinks)
  
  def __setstate__(self, state) :

    self.peakLinks, self.notFoundSimulatedPeaks, self.spinSystem1, self.spinSystem2 = state # , self.notFoundPeakLinks = state #self.notFoundPeakLinks
    
  def getContributingResonances(self) :
    
    resonances = []
    
    for pl in self.peakLinks :
      
      resonances.extend(pl.getResonances())
      
    return set(resonances)  
  
  def getAllResonances(self) :                        #TODO: change this, first of all it is not used and second of all self.notFoundPeakLinks is empty now.
  
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
    
    cdef SimulatedPeak simPeak
    cdef SimulatedPeakContrib contrib
    cdef SpinSystem spinSystem

    notFoundPeakLinks = []
    
    spinSystems = {1:self.spinSystem1,2:self.spinSystem2}
    
    for simPeak in self.notFoundSimulatedPeaks :
      
      resonances = []
      
      for contrib in simPeak.simulatedPeakContribs :
        
        spinSystem = spinSystems[contrib.firstOrSecondResidue]
        
        resonances.append(spinSystem.getResonanceForAtomName(contrib.atomName))
      
      notFoundPeakLinks.append(PeakLink(None, simPeak,resonances,0.0))
    
    return notFoundPeakLinks
  
  def getAllPeakLinks(self) :
    
    #return self.peakLinks + self.notFoundPeakLinks
    
    return self.getPeakLinks() + self.getNotFoundPeakLinks()
  
  def getResidues(self) :
    
    return (self.residue1, self.residue2)
  
  def getSpinSystems(self) :
    
    return (self.spinSystem1, self.spinSystem2)
  
  