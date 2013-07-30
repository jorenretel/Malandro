
cdef class Chain :

  cdef public list residues 
    
  cdef dict residuesByCcpCode
  
  cdef dict residueTypeFrequencyDict

  cdef object ccpnChain 
  
  cdef object pyChain
  
  cdef str ccpnChainCode, molSystemCode
  
  def __init__(self, ccpnChain):
    
    self.ccpnChain = ccpnChain
    
    self.ccpnChainCode = ccpnChain.code
    
    self.molSystemCode = ccpnChain.molSystem.code

    self.residues = []
    
    self.residueTypeFrequencyDict = {}
    
    self.residuesByCcpCode = {}

    self.setupResidues()
    
    self.countResidueTypeFrequency()
    
    self.linkResiduesTogether()
    
    self.addDummyResiduesAtEnds()
    
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
    cdef SpinSystem spinSystem
    cdef aResidue firstResidue
    cdef aResidue lastResidue
    
    firstResidue = self.residues[0]
    lastResidue = self.residues[-1]
    
    res = aResidue(self,None)
    spinSystem = SpinSystem()
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
      
  def __reduce__(self) :

    return (generalFactory, (type(self),), self.__getstate__())
    
  def __getstate__(self) :
    
    return (self.residues, self.molSystemCode, self.ccpnChainCode)
  
  def __setstate__(self, state) :

    self.residues, self.molSystemCode, self.ccpnChainCode = state
   
  def connectToProject(self, project) :
    '''(Re)set the connection to the chain in the analysis project,
        i.e. set self.ccpnChain atribute.
    '''
    cdef aResidue residue
    
    molSystem = project.findFirstMolSystem(code=self.molSystemCode)
    
    if molSystem :
    
      ccpnChain = molSystem.findFirstChain(code=self.ccpnChainCode)
      
      if ccpnChain :
        
        self.ccpnChain = ccpnChain
        
    if not self.ccpnChain :
      
      print 'It seems like the ' + self.molSystemCode + ' ' + self.ccpnChainCode + ' chain was removed.'
      return
    
    for residue in self.residues :
      
      residue.connectToProject()
      
      
    
  def getResidues(self) :
    
    return self.residues
