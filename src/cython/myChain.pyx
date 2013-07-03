
cdef class myChain :

  cdef list residues 
    
  cdef dict residuesByCcpCode
  
  cdef dict residueTypeFrequencyDict

  cdef object ccpnChain 
  
  cdef object pyChain
  

  
  def __init__(self, ccpnChain):
    
    self.ccpnChain = ccpnChain

    self.residues = []
    
    self.residueTypeFrequencyDict = {}
    
    self.residuesByCcpCode = {}

    self.setupResidues()
    
    self.countResidueTypeFrequency()
    
    self.linkResiduesTogether()
    
    self.addDummyResiduesAtEnds()
    
  #def __getstate__(self):
  #  state = dict(self.__dict__)
  #  if 'ccpnChain' in state :
  #    del state['ccpnChain']

  cdef void setupResidues(self):
  
    for res in self.ccpnChain.sortedResidues() :
      
      newresidue = aResidue(self, res)

      self.residues.append(newresidue)
      
      if res.ccpCode in self.residuesByCcpCode :
          
        self.residuesByCcpCode[res.ccpCode].append(newresidue)
        
      else :
        
        self.residuesByCcpCode[res.ccpCode] = [newresidue]
        
  cdef void addDummyResiduesAtEnds(self) :
    '''
    Put a residue before the beginning and after the end. Why do this:
    this saves a lot of checks during the annealing procedure for not getting
    a None-type object for .previousResidue or .nextResidue. It is not added to
    self.residues, so normally you don't even notice those dummy residues are
    there at all.
    '''
    cdef aResidue res
    cdef mySpinSystem spinSystem
    cdef aResidue firstResidue
    cdef aResidue lastResidue
    
    firstResidue = self.residues[0]
    lastResidue = self.residues[-1]
    
    res = aResidue(self,None)
    spinSystem = mySpinSystem()
    #spinSystem.isJoker = True
    
    res.currentSpinSystemAssigned = spinSystem
    
    firstResidue.previousResidue = res
    lastResidue.nextResidue = res
    
  cdef void countResidueTypeFrequency(self):
    
    cdef aResidue res
    
    for res in self.residues :
      
      if res.ccpCode in self.residueTypeFrequencyDict:
        
        self.residueTypeFrequencyDict[res.ccpCode] = (self.residueTypeFrequencyDict[res.ccpCode] + 1)
        
      else :
        
        self.residueTypeFrequencyDict[res.ccpCode] = 1
        
  cdef void linkResiduesTogether(self):
    
    cdef list residues
    
    cdef int i
    
    cdef aResidue res
    
    cdef aResidue nextResidue
    
    residues = self.residues
    
    for i,  res in enumerate(residues[:-1])  :
      
      nextResidue = residues[i+1]
      
      res.nextResidue = nextResidue
      
      nextResidue.previousResidue = res

  cdef void createPythonStyleObject(self):
    
    cdef aResidue res
    
    self.pyChain = pyChain(self.ccpnChain)
    
    for res in self.residues :
      
      self.pyChain.residues.append(res.pyResidue)
      
  def __reduce__(self) :

    return (generalFactory, (type(self),), self.__getstate__())
    
  def __getstate__(self) :
    
    return (self.residues)
  
  def __setstate__(self, state) :

    self.residues = state
