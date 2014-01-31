
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
  #cdef double minIsoFrac
  cdef object nmrProject, project

  def __init__(self, Malandro auto):
    
    self.auto = auto
    self.spectra = []
    self.chain = None
    self.spinSystems = {}

    
    aminoAcids = ['Ala', 'Arg', 'Asn', 'Asp', 'Cys', 'Gln', 'Glu',
                  'Gly', 'His', 'Ile', 'Leu', 'Lys', 'Met', 'Phe',
                  'Pro', 'Ser', 'Thr', 'Trp', 'Tyr', 'Val']
    
    for aa in aminoAcids :
      
      self.spinSystems[aa] = []

    #self.minIsoFrac = auto.minIsoFrac

  def setupChain(self, ccpnChain):
    '''Generates new Chain by based on a chain object from CCPN.'''
    
    self.chain = Chain(ccpnChain)
        
  def setupSpectra(self, selectedSpectra):
    '''Sets up spectra objects for all selected spectra.'''
    
    for temporary_spectrum_object in selectedSpectra :

      newspectrum = Spectrum(temporary_spectrum_object)
      
      newspectrum.DataModel = self

      self.spectra.append(newspectrum)

  def setupSpinSystems(self, resonanceGroups, useAssignments=True, useTentative=True,
                       useType=True,includeUntypedSpinSystems=True, minTypeScore=0.0, makeJokers=False) :
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
    
    cdef SpinSystem spinSystem
        
    assignedResidues = set()
    
    assignedResonanceGroups = []
    tentativeResonanceGroups = []
    typedResonanceGroups = []
    multipleTypeResonanceGroups = []
    untypedResonanceGroups = []
    
    ignoreAssignedResidues = useAssignments
    
    for resonanceGroup in resonanceGroups :
      # Sorting the resonanceGroups by the type of
      # information they contain.

      if not resonanceGroup.resonances :
        # resonanceGroups that contain no resonances are basically placeholders
        # in the analysis project that we are not interested in.
        continue
      
      if resonanceGroup.residue:
        # resonanceGroup is assigned to a residue.
        if resonanceGroup.residue.chain is self.chain.ccpnChain:
          # we should not use resonanceGroups that
          # are assigned to another chain.
          assignedResonanceGroups.append(resonanceGroup)
        
      elif resonanceGroup.residueProbs:
        # resonanceGroup has tentative assignment to residues
        # that can be used to constraint the possible residue
        # assignments.
        tentativeResonanceGroups.append(resonanceGroup)
      
      elif resonanceGroup.ccpCode:
        # resonanceGroup is just amino acid typed. 
        typedResonanceGroups.append(resonanceGroup)
  
      elif resonanceGroup.residueTypeProbs:
        # resonanceGroup has a set of options for the amino acid type.
        multipleTypeResonanceGroups.append(resonanceGroup)

      elif includeUntypedSpinSystems:
        # resonanceGroup does not hold information
        # about residue assignment or type.
        untypedResonanceGroups.append(resonanceGroup)

    for resonanceGroup in assignedResonanceGroups:
      
      seqCode = int(resonanceGroup.residue.seqCode)
      residue = self.chain.residues[seqCode-1]
      residues = [residue]
      ccpCodes = [residue.ccpCode] #getResidueCode(resonanceGroup.residue.molResidue)
      
      if useAssignments:
        self._makeSpinsystem(resonanceGroup=resonanceGroup, residues=residues)
        assignedResidues.add(residue)
      elif useType:
        self._makeSpinSystem(resonanceGroup, ccpCodes)
      else:
        self._makeSpinSystem(resonanceGroup, minTypeScore=minTypeScore)
 
    for resonanceGroup in tentativeResonanceGroups:
      
      residues = []
      
      for residueProb in resonanceGroup.residueProbs:
        if residueProb.weight:
          ccpnResidue = residueProb.possibility
          if ccpnResidue.chain is self.chain.ccpnChain :
            seqCode = ccpnResidue.seqCode
            residue = self.chain.residues[seqCode-1]
            residues.append(residue)

      if residues :
        
        if useTentative:
          self._makeSpinSystem(resonanceGroup, residues)
        elif useType:
          ccpCodes = [residue.ccpCode for residue in residues]
          self._makeSpinSystem(resonanceGroup, ccpCodes, unAllowedResidues=assignedResidues)
        else :
          self._makeSpinSystem(resonanceGroup, unAllowedResidues=assignedResidues, minTypeScore=minTypeScore)

    for resonanceGroup in typedResonanceGroups:
      
      if useType:
        ccpCodes = [resonanceGroup.ccpCode]
        self._makeSpinSystem(resonanceGroup, ccpCodes, unAllowedResidues=assignedResidues)
      else:
        self._makeSpinSystem(resonanceGroup, unAllowedResidues=assignedResidues, minTypeScore=minTypeScore)
    
    for resonanceGroup in multipleTypeResonanceGroups:
      
      if useType :
        ccpCodes = [residueTypeProb.possibility.ccpCode for residueTypeProb in resonanceGroup.residueTypeProbs]
        self._makeSpinSystem(resonanceGroup, ccpCodes, unAllowedResidues=assignedResidues)
      else:
        self._makeSpinSystem(resonanceGroup, unAllowedResidues=assignedResidues, minTypeScore=minTypeScore)
        
    for resonanceGroup in untypedResonanceGroups:
      
      self._makeSpinSystem(resonanceGroup, unAllowedResidues=assignedResidues, minTypeScore=minTypeScore)
      
    if makeJokers:
      
      for ccpCode, residues in self.chain.residuesByCcpCode:
      
        for i in range(len(residues)):
          
          self._makeSpinSystem(None,ccpCodes=[ccpCode], unAllowedResidues=assignedResidues)
      
    
  def _makeSpinsystem(self, resonanceGroup, residues=None, ccpCodes=None, unAllowedResidues=None, minTypeScore=0.0):
    
    allowedResidues = set()
    if residues:
      ccpCodes = set([residue.ccpCode for residue in residues])
      allowedResidues = set(residues)
    elif ccpCodes:
      ccpCodes = set(ccpCodes)
      allowedResidues = set([self.chain.residuesByCcpCode[ccpCode] for ccpCode in ccpCodes])
    else:
      aminoAcidProbs = runAminoAcidTyping(resonanceGroup, self.auto.shiftsList, self.chain.ccpnChain, minTypeScore)
      ccpCodes = set(aminoAcidProbs.keys())
      allowedResidues = set([self.chain.residuesByCcpCode[ccpCode] for ccpCode in ccpCodes])
 
    if unAllowedResidues:
      allowedResidues -= unAllowedResidues
    
    newSpinSystem = SpinSystem(DataModel=self, ccpnResonanceGroup=resonanceGroup, allowedResidues=allowedResidues, ccpCodes=ccpCodes)
    newSpinSystem.allowedResidues = allowedResidues
    
    for residue in allowedResidues:
      residue.allowedSpinSystems.add(newSpinSystem)
    for ccpCode in ccpCodes:
      self.addToDictWithLists(self.spinSystems, ccpCode, newSpinSystem)
  
 
  def setupLinks(self) :
    '''Setup all links between spin systems. Two dicts are created.
       One intraDict for intra-residual links, which contains intra-residual
       peaks (or peak links, to be exact). In this case the key in the
       dict is just the serial of the spin system. The linkDict contains
       sequential links between spin system. Very simple keys are created
       from both spin system serials. Jokers spin systems are excluded.
       
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
        
        if spinSystemA.isJoker:
          continue
        
        intraDict[spinSystemA.spinSystemNumber] = SpinSystemLink(residue1=resA,residue2=resA,spinSystem1=spinSystemA,spinSystem2=spinSystemA)
        
        for spinSystemB in self.spinSystems[ccpCodeB] :
          
          if spinSystemB.isJoker:
            continue
          
          linkDict[spinSystemA.spinSystemNumber*10000+spinSystemB.spinSystemNumber] = SpinSystemLink(residue1=resA,residue2=resB,spinSystem1=spinSystemA,spinSystem2=spinSystemB)
    
    # Also doing the last residue here, otherwise not reached in loop.
    intraDict = resB.intraDict
    
    for spinSystemB in self.spinSystems[ccpCodeB] :
      
      if spinSystemB.isJoker:
        continue
      
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
  
  def getSpinSystemSet(assigned=True, typed=True, ccpCode=None):
    pass
  
  def getSpinSystems(self) :
    '''Return spin systems in dict, keys are three-letter amino acid codes'''
    
    return self.spinSystems
  
  def getSpectra(self) :
    '''Return list with spectra.'''
    return self.spectra