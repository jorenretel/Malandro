
cdef class myDataModel :
  '''Root of the  model for Malandro. Most objects in the model are bascically
     wrappers around objects from CCPN. This is done like this, so features can
     be added to the objects without changing the original object.
     This class does not contain much more than references to spectra, spin systems,
     and the chain. Most methods are involved in setting up these objects.
     
  '''

  cdef list spectra
  cdef public Chain chain
  cdef dict spinSystems, previouslyAssignedSpinSystems, justTypedSpinSystems, tentativeSpinSystems, untypedSpinSystems, typeProbSpinSystems, allSpinSystemsWithoutAssigned, jokerSpinSystems
  cdef Malandro auto  
  cdef double minIsoFrac
  cdef object nmrProject, project

  def __init__(self, Malandro auto):
    
    self.auto = auto
    self.spectra = []
    self.chain = None
    self.spinSystems = {}
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
      
      self.spinSystems[aa] = []
      self.previouslyAssignedSpinSystems[aa] = []
      self.justTypedSpinSystems[aa] = []
      self.tentativeSpinSystems[aa] = []
      self.typeProbSpinSystems[aa] = []
      self.untypedSpinSystems[aa] = []
      self.allSpinSystemsWithoutAssigned[aa] = []
      self.jokerSpinSystems[aa] = []

    
    self.minIsoFrac = auto.minIsoFrac

  def setupChain(self):
    '''Generates new Chain by based on a chain object from CCPN.'''
    
    self.chain = Chain(self.auto.chain)
        
  def setupSpectra(self):
    '''Sets up spectra objects for all selected spectra.'''
    
    for temporary_spectrum_object in self.auto.selectedSpectra :

      newspectrum = Spectrum(temporary_spectrum_object)
      
      newspectrum.DataModel = self

      self.spectra.append(newspectrum)

  def setupSpinSystems(self, minTypeScore=1.0) :
    '''Sets up spin system objects based on CCPN resonanceGroup.
       kwarg: minTypeScore can be passed. This is used when the possible
              residue types of the spin system has to be determined. The residue
              typing algorithm of CCPN analysis is used for this. The minTypeScore
              is a percentage (i.e. 0-100). All residue types scoring higher than
              this cut-off will be consider possible.
       Spin systems can have different levels of assignment. Depending on their
       assignment state they are sorted into different dicts, to be easily accessible
       later.
       
    '''
    
    cdef SpinSystem newSpinSystem
    cdef bint reTypeSpinSystems
    
    reTypeSpinSystems = self.auto.reTypeSpinSystems
    
    for resonanceGroup in self.auto.nmrProject.resonanceGroups :   # taking all spinsystems in the project, except for the ones that have no resonances

      if not resonanceGroup.resonances :
        
        continue
      
      if resonanceGroup.residue and resonanceGroup.residue.chain is self.chain.ccpnChain and not reTypeSpinSystems:                                  # SpinSystem is assigned to a residue in the selected chain
        
        seqCode = int(resonanceGroup.residue.seqCode)
        ccpCode = getResidueCode(resonanceGroup.residue.molResidue)
        
        newSpinSystem = SpinSystem(DataModel=self, ccpnResonanceGroup=resonanceGroup, ccpnSeqCode = seqCode, ccpCode=ccpCode)
        
        self.addToDictWithLists(self.previouslyAssignedSpinSystems, ccpCode, newSpinSystem)
        self.addToDictWithLists(self.spinSystems, ccpCode, newSpinSystem)
        
      elif resonanceGroup.residueProbs and not reTypeSpinSystems:                                                                                                      # SpinSystem has one or more tentative assignments. Got this piece out of EditSpinSystem.py in popups.
        
        ccpCodes = []
        seqCodes = []
        
        for residueProb in resonanceGroup.residueProbs:
          if not residueProb.weight:
            continue
            
          residue = residueProb.possibility
          
          seq = residue.seqCode
          resCode = residue.ccpCode
          
          if residue.chain is self.chain.ccpnChain :
  
            ccpCodes.append(resCode)                                                                   # Only consider the tentative assignments to residues in the selected chain.                                                         
            seqCodes.append(seq)
            
        if seqCodes :

          newSpinSystem = SpinSystem(DataModel=self, ccpnResonanceGroup=resonanceGroup, tentativeSeqCodes = seqCodes,tentativeCcpCodes=ccpCodes)
          
          for ccpCode in set(ccpCodes) :
          
            self.addToDictWithLists(self.tentativeSpinSystems, ccpCode, newSpinSystem)
            self.addToDictWithLists(self.spinSystems, ccpCode, newSpinSystem)
            self.addToDictWithLists(self.allSpinSystemsWithoutAssigned, ccpCode, newSpinSystem)
      
      elif resonanceGroup.ccpCode and not reTypeSpinSystems:                                                                                                              # Residue is just Typed
  
        ccpCode = resonanceGroup.ccpCode
        
        newSpinSystem = SpinSystem(DataModel=self, ccpnResonanceGroup=resonanceGroup, ccpCode=ccpCode)
        
        self.addToDictWithLists(self.justTypedSpinSystems, ccpCode, newSpinSystem)
        self.addToDictWithLists(self.spinSystems, ccpCode, newSpinSystem)
        self.addToDictWithLists(self.allSpinSystemsWithoutAssigned, ccpCode, newSpinSystem)

      elif resonanceGroup.residueTypeProbs and not reTypeSpinSystems:
        
        typeProbCcpCodes = [residueTypeProb.possibility.ccpCode for residueTypeProb in resonanceGroup.residueTypeProbs]
        
        newSpinSystem = SpinSystem(DataModel=self, ccpnResonanceGroup=resonanceGroup,typeProbCcpCodes=typeProbCcpCodes)
        
        for ccpCode in typeProbCcpCodes :
        
          self.addToDictWithLists(self.typeProbSpinSystems, ccpCode, newSpinSystem)
          self.addToDictWithLists(self.spinSystems, ccpCode, newSpinSystem)
          self.addToDictWithLists(self.allSpinSystemsWithoutAssigned, ccpCode, newSpinSystem)
        
      elif self.auto.typeSpinSystems or reTypeSpinSystems:
        
        newSpinSystem = SpinSystem(DataModel=self, ccpnResonanceGroup=resonanceGroup, typeSpinSystem=True, minTypeScore=minTypeScore)
        
        for ccpCode in newSpinSystem.aminoAcidProbs.keys() :
          
          self.addToDictWithLists(self.untypedSpinSystems, ccpCode, newSpinSystem)
          self.addToDictWithLists(self.spinSystems, ccpCode, newSpinSystem)
          self.addToDictWithLists(self.allSpinSystemsWithoutAssigned, ccpCode, newSpinSystem)
  
  def setupLinks(self) :
    '''Setup all links between spin systems. Two dicts are created.
       One intraDict for intra-residual links, which contains intra-residual
       peaks (or peak links, to be exact). In this case the key in the
       dict is just the serial of the spin system. The linkDict contains
       sequential links between spin system. Very simple keys are created
       from both spin system serials.
       
    '''
    cdef Residue resA, resB
    cdef dict linkDict, intraDict
    cdef str ccpCodeA, ccpCodeB
    cdef SpinSystem spinSystemA, spinSystemB
    cdef list residues
    
    residues = self.chain.residues
    
    for resA, resB in zip(residues,residues[1:]):

      ccpCodeA = resA.ccpCode
      ccpCodeB = resB.ccpCode
      
      linkDict = resA.linkDict
      intraDict = resA.intraDict
      
      for spinSystemA in self.spinSystems[ccpCodeA] :
        
        intraDict[spinSystemA.spinSystemNumber] = SpinSystemLink(residue1=resA,residue2=resA,spinSystem1=spinSystemA,spinSystem2=spinSystemA)
        
        for spinSystemB in self.spinSystems[ccpCodeB] :
          
          linkDict[spinSystemA.spinSystemNumber*10000+spinSystemB.spinSystemNumber] = SpinSystemLink(residue1=resA,residue2=resB,spinSystem1=spinSystemA,spinSystem2=spinSystemB)
    
    # Also doing the last residue here, otherwise not reached in loop.
    intraDict = resB.intraDict
    
    for spinSystemB in self.spinSystems[ccpCodeB] :
    
      intraDict[spinSystemB.spinSystemNumber] = SpinSystemLink(residue1=resB,residue2=resB,spinSystem1=spinSystemB,spinSystem2=spinSystemB)
          
  cdef void addToDictWithLists(self, dict dictToAddTo, key,value) :
    
    if key in dictToAddTo :       
      
      dictToAddTo[key].append(value)
      
    else :                                                                                                                          
      
      dictToAddTo[key] = [value]

  cdef list getResonancesForAtomSiteName(self, str atomSiteName) :
    '''Return all resonances of all spin systems in the model
       that fit the atomSite described by the atomSiteName.
       
    '''
    
    cdef SpinSystem spinSystem
    
    resonances = []
    
    spinSystems = []
        
    for spinSystemList in self.spinSystems.values() :
    
      spinSystems.extend(spinSystemList)
      
      
    spinSystems = set(spinSystems)
    
    for spinSystem in spinSystems :
      
      resonances.extend(spinSystem.getResonancesForAtomSiteName(atomSiteName))
      
    return resonances  
          
  def __reduce__(self) :

    return (generalFactory, (myDataModel,), self.__getstate__())
    
  def __getstate__(self) :

    return (self.spectra, self.chain, self.spinSystems)
  
  def __setstate__(self, state) :

    self.spectra, self.chain, self.spinSystems = state
    
  def connectToProject(self, project, nmrProject) :
    
    '''This method can be called after an unpickling from file. Because not the whole ccpn analysis project
       should be pickled (it would also not be possible), the link between the malandro objects and the analysis
       objects is lost. By running this method, the attributes that link to objects in the ccpn project are
       set again. Of course when the project has changed and objects have been removed,
       they can not be set again.
       
    '''
    
    self.project = project
    self.nmrProject = nmrProject
    
    # Connect the chain, residues and atoms to the corresponding chain, residues and atoms in ccpn.
    self.chain.connectToProject(project)
    
    # Connect spinsystems and resonances
    spinSystems = set([spinSystem for sublist in self.spinSystems.values() for spinSystem in sublist])

    for spinSystem in spinSystems :
      
      spinSystem.connectToProject(nmrProject)
      
    # Connect spectra, peaks and peak dimensions  
    for spectrum in self.spectra :
      
      spectrum.connectToProject(nmrProject)
      
  def getChain(self) :
    '''Resturn chain.'''
  
    return self.chain
  
  def getSpinSystems(self) :
    '''Return spin systems in dict, keys are three-letter amino acid codes'''
    
    return self.spinSystems
  
  def getSpectra(self) :
    '''Return list with spectra.'''
    return self.spectra