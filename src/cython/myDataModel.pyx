
cdef class myDataModel :

  cdef list spectra

  cdef myChain myChain

  cdef dict mySpinSystems, previouslyAssignedSpinSystems, justTypedSpinSystems, tentativeSpinSystems, untypedSpinSystems, typeProbSpinSystems, allSpinSystemsWithoutAssigned, jokerSpinSystems
  
  cdef Malandro auto
  
  cdef object pyDataModel
  
  cdef double minIsoFrac


  def __init__(self, Malandro auto):
    
    self.auto = auto

    self.spectra = []

    self.myChain = None
        
    self.mySpinSystems = {}
      
    self.previouslyAssignedSpinSystems = {}
        
    self.justTypedSpinSystems = {}
    
    self.tentativeSpinSystems = {}
    
    self.typeProbSpinSystems = {}
    
    self.untypedSpinSystems = {}
    
    self.allSpinSystemsWithoutAssigned = {}
    
    self.jokerSpinSystems = {}
    
    self.minIsoFrac = auto.minIsoFrac

  def setupChain(self):
    
    self.myChain = myChain(self.auto.chain)
        
  def setupSpectra(self):

    for temporary_spectrum_object in self.auto.selectedSpectra :

      newspectrum = aSpectrum(temporary_spectrum_object)
      
      newspectrum.DataModel = self

      self.spectra.append(newspectrum)

  def setupSpinSystems(self) :
    
    cdef mySpinSystem newSpinSystem
    
    for resonanceGroup in self.auto.nmrProject.resonanceGroups :   # taking all spinsystems in the project, except for the ones that have no resonances

      if not resonanceGroup.resonances :
        
        continue
      
      if resonanceGroup.residue and resonanceGroup.residue.chain is self.myChain.ccpnChain :                                  # SpinSystem is assigned to a residue in the selected chain
        
        seqCode = int(resonanceGroup.residue.seqCode)
        ccpCode = getResidueCode(resonanceGroup.residue.molResidue)
        
        newSpinSystem = mySpinSystem(DataModel=self, ccpnResonanceGroup=resonanceGroup, ccpnSeqCode = seqCode, ccpCode=ccpCode)
        
        self.addToDictWithLists(self.previouslyAssignedSpinSystems, ccpCode, newSpinSystem)
        self.addToDictWithLists(self.mySpinSystems, ccpCode, newSpinSystem)
        
      elif resonanceGroup.residueProbs :                                                                                                      # SpinSystem has one or more tentative assignments. Got this piece out of EditSpinSystem.py in popups.
        
        ccpCodes = []
        seqCodes = []
        
        for residueProb in resonanceGroup.residueProbs:
          if not residueProb.weight:
            continue
            
          residue = residueProb.possibility
          
          seq = residue.seqCode
          resCode = residue.ccpCode
          
          if residue.chain is self.myChain.ccpnChain :
  
            ccpCodes.append(resCode)                                                                   # Only consider the tentative assignments to residues in the selected chain.                                                         
            seqCodes.append(seq)
            
        if seqCodes :

          newSpinSystem = mySpinSystem(DataModel=self, ccpnResonanceGroup=resonanceGroup, tentativeSeqCodes = seqCodes,tentativeCcpCodes=ccpCodes)
          
          for ccpCode in set(ccpCodes) :
          
            self.addToDictWithLists(self.tentativeSpinSystems, ccpCode, newSpinSystem)
            self.addToDictWithLists(self.mySpinSystems, ccpCode, newSpinSystem)
            self.addToDictWithLists(self.allSpinSystemsWithoutAssigned, ccpCode, newSpinSystem)
      
      elif resonanceGroup.ccpCode :                                                                                                              # Residue is just Typed
  
        ccpCode = resonanceGroup.ccpCode
        
        newSpinSystem = mySpinSystem(DataModel=self, ccpnResonanceGroup=resonanceGroup, ccpCode=ccpCode)
        
        self.addToDictWithLists(self.justTypedSpinSystems, ccpCode, newSpinSystem)
        self.addToDictWithLists(self.mySpinSystems, ccpCode, newSpinSystem)
        self.addToDictWithLists(self.allSpinSystemsWithoutAssigned, ccpCode, newSpinSystem)

      elif resonanceGroup.residueTypeProbs :
        
        typeProbCcpCodes = [residueTypeProb.possibility.ccpCode for residueTypeProb in resonanceGroup.residueTypeProbs]
        
        newSpinSystem = mySpinSystem(DataModel=self, ccpnResonanceGroup=resonanceGroup,typeProbCcpCodes=typeProbCcpCodes)
        
        for ccpCode in typeProbCcpCodes :
        
          self.addToDictWithLists(self.typeProbSpinSystems, ccpCode, newSpinSystem)
          self.addToDictWithLists(self.mySpinSystems, ccpCode, newSpinSystem)
          self.addToDictWithLists(self.allSpinSystemsWithoutAssigned, ccpCode, newSpinSystem)
        
      elif self.auto.typeSpinSystems :
        
        newSpinSystem = mySpinSystem(DataModel=self, ccpnResonanceGroup=resonanceGroup, typeSpinSystem=True)
        
        for ccpCode in newSpinSystem.aminoAcidProbs.keys() :
          
          self.addToDictWithLists(self.untypedSpinSystems, ccpCode, newSpinSystem)
          self.addToDictWithLists(self.mySpinSystems, ccpCode, newSpinSystem)
          self.addToDictWithLists(self.allSpinSystemsWithoutAssigned, ccpCode, newSpinSystem)
        
  cdef void addToDictWithLists(self, dict dictToAddTo, key,value) :
    
    if key in dictToAddTo :       
      
      dictToAddTo[key].append(value)
      
    else :                                                                                                                          
      
      dictToAddTo[key] = [value]

  cdef list getResonancesForAtomSiteName(self, str atomSiteName) :
    
    cdef mySpinSystem spinSystem
    
    resonances = []
    
    spinSystems = []
        
    for spinSystemList in self.mySpinSystems.values() :
    
      spinSystems.extend(spinSystemList)
      
      
    spinSystems = set(spinSystems)
    
    for spinSystem in spinSystems :
      
      resonances.extend(spinSystem.getResonancesForAtomSiteName(atomSiteName))
      
    return resonances  
    
  cdef void createPythonStyleObject(self):
    
    cdef aSpectrum spec
    
    cdef mySpinSystem spinSystem 
    
    cdef list listWithSpinSystems
    
    cdef list listWithPySpinSystems
    
    cdef str key
    
    self.pyDataModel = pyDataModel()
    
    self.pyDataModel.chain = self.myChain.pyChain
    
    for spec in self.spectra :
    
      self.pyDataModel.spectra.append(spec.pySpectrum)
      
    for key, listWithSpinSystems in self.mySpinSystems.items():
    
      listWithPySpinSystems = []
    
      for spinSystem in listWithSpinSystems :
        
        listWithPySpinSystems.append(spinSystem.pySpinSystem)
        
      self.pyDataModel.spinSystems[key] = listWithPySpinSystems