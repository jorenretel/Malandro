
cdef class myDataModel :

  cdef list spectra

  cdef public myChain myChain

  cdef dict mySpinSystems, previouslyAssignedSpinSystems, justTypedSpinSystems, tentativeSpinSystems, untypedSpinSystems, typeProbSpinSystems, allSpinSystemsWithoutAssigned, jokerSpinSystems
  
  cdef Malandro auto
  
  cdef double minIsoFrac
  
  cdef object nmrProject, project


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
    
    aminoAcids = ['Ala', 'Arg', 'Asn', 'Asp', 'Cys', 'Gln', 'Glu',
                  'Gly', 'His', 'Ile', 'Leu', 'Lys', 'Met', 'Phe',
                  'Pro', 'Ser', 'Thr', 'Trp', 'Tyr', 'Val']
    
    for aa in aminoAcids :
      
      self.mySpinSystems[aa] = []
      self.previouslyAssignedSpinSystems[aa] = []
      self.justTypedSpinSystems[aa] = []
      self.tentativeSpinSystems[aa] = []
      self.typeProbSpinSystems[aa] = []
      self.untypedSpinSystems[aa] = []
      self.allSpinSystemsWithoutAssigned[aa] = []
      self.jokerSpinSystems[aa] = []

    
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
  
  def setupLinks(self) :
  
    cdef aResidue resA, resB
    cdef dict linkDict, intraDict
    cdef str ccpCodeA, ccpCodeB
    cdef mySpinSystem spinSystemA, spinSystemB
    cdef list residues
    
    residues = self.myChain.residues
    
    for resA, resB in zip(residues,residues[1:]):

      ccpCodeA = resA.ccpCode
      ccpCodeB = resB.ccpCode
      
      linkDict = resA.linkDict
      intraDict = resA.intraDict
      
      for spinSystemA in self.mySpinSystems[ccpCodeA] :
        
        intraDict[spinSystemA.spinSystemNumber] = spinSystemLink(residue1=resA,residue2=resA,spinSystem1=spinSystemA,spinSystem2=spinSystemA)
        
        for spinSystemB in self.mySpinSystems[ccpCodeB] :
          
          linkDict[spinSystemA.spinSystemNumber*10000+spinSystemB.spinSystemNumber] = spinSystemLink(residue1=resA,residue2=resB,spinSystem1=spinSystemA,spinSystem2=spinSystemB)
          
          
          
          
      
      
        
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
          
  def __reduce__(self) :

    return (generalFactory, (myDataModel,), self.__getstate__())
    
  def __getstate__(self) :

    return (self.spectra, self.myChain, self.mySpinSystems)
  
  def __setstate__(self, state) :

    self.spectra, self.myChain, self.mySpinSystems = state
    
  def connectToProject(self, project, nmrProject) :
    
    '''This method can be called after a unpickling from file. Because not the whole ccpn analysis project
       should be pickled (it would also not be possible), the link between the malandro objects and the analysis
       objects is lost. By running this method, the attributes that link to objects in the ccpn project are
       set again. Of course when the project has changed a lot in the mean time and objects have been removed,
       they can not be set again.
    '''
    
    self.project = project
    self.nmrProject = nmrProject
    
    # Connect the chain, residues and atoms to the corresponding chain, residues and atoms in ccpn.
    self.myChain.connectToProject(project)
    
    # Connect spinsystems and resonances
    spinSystems = set([spinSystem for sublist in self.mySpinSystems.values() for spinSystem in sublist])

    for spinSystem in spinSystems :
      
      spinSystem.connectToProject(nmrProject)
      
    # Connect spectra, peaks and peak dimensions  
    for spectrum in self.spectra :
      
      spectrum.connectToProject(nmrProject)
      
  def getChain(self) :
  
    return self.myChain
  
  def getSpinSystems(self) :
    
    return self.mySpinSystems
  
  def getSpectra(self) :
    
    return self.spectra