
cdef class mySpinSystem :
  
  cdef myDataModel DataModel

  cdef public aResidue currentResidueAssignment
  
  cdef public int spinSystemNumber
  
  cdef public list exchangeSpinSystems
  
  cdef str ccpCode
  
  cdef dict resonanceDict
  
  cdef object ccpnResonanceGroup
  
  cdef int ccpnSeqCode
  
  cdef bint isJoker
  
  cdef list solutions
  
  #cdef list userDefinedSolutions
  
  cdef list tentativeCcpCodes
  
  cdef list tentativeSeqCodes
  
  cdef list typeProbCcpCodes
  
  cdef set allowedResidues
  
  cdef dict aminoAcidProbs
  
  cdef object pySpinSystem
  
  cdef dict resonancesByAtomSiteName
  
  cdef bint [:] allowedResidueView
  

  def __init__(self, DataModel=None, ccpnResonanceGroup=None, ccpnSeqCode = 0, ccpCode=None,
               tentativeSeqCodes = [],tentativeCcpCodes=[] , typeProbCcpCodes=[], typeSpinSystem=False) :
    
    self.DataModel = DataModel
    self.ccpnResonanceGroup = ccpnResonanceGroup
    self.ccpCode = ccpCode
    self.ccpnSeqCode = ccpnSeqCode
    self.typeProbCcpCodes = typeProbCcpCodes
    self.tentativeCcpCodes = tentativeCcpCodes
    self.tentativeSeqCodes = tentativeSeqCodes
    self.resonanceDict = {}
    
    self.currentResidueAssignment = None   
    
    self.isJoker = False
    
    self.solutions = []
    
    self.aminoAcidProbs = {}
    
    self.exchangeSpinSystems = []
    
    self.allowedResidues = set()                                      # This is later on going to be a Frozen Set for fast membership testing during the annealing run. If the set is empty that means everything residue is allowed, if has members, only these residues are allowed.
    
    self.resonancesByAtomSiteName = {}

    if not ccpnResonanceGroup :
      
      self.isJoker = True
      
      return
    
    self.spinSystemNumber = ccpnResonanceGroup.serial
    

      
    if typeSpinSystem : 

      self.runAminoAcidTyping()
 
    self.setupResonances()

    self.groupResonancesByAtomSite()


  cdef set getCcpCodes(self) :
    
    if self.ccpCode :
      
      ccpCodes = [self.ccpCode]
    
    else :
      
      ccpCodes = self.tentativeCcpCodes or self.typeProbCcpCodes or self.aminoAcidProbs.keys() or []
      
    return set(ccpCodes)  
      
  cdef void setupResonances(self) :

    for resonance in self.ccpnResonanceGroup.resonances :
          
      if resonance.assignNames :
    
        newResonance = myResonance(self,resonance)
        
        self.resonanceDict[newResonance.atomName] = newResonance

  cdef void runAminoAcidTyping(self) :
    
    shiftList = self.DataModel.auto.shiftList
    ccpnChain = self.DataModel.myChain.ccpnChain
    
    shifts = []
    for resonance in self.ccpnResonanceGroup.resonances:
      if resonance.isotopeCode in ('1H',  '13C',  '15N') :
        shift = resonance.findFirstShift(parentList=shiftList)
        if shift:
          shifts.append(shift)
  
    scores = getShiftsChainProbabilities(shifts, ccpnChain)
    total = sum(scores.values())

    scoreList = []
    
    scoreDict = {}
    
    if total:
        
      for ccpCode, score in scores.items() :
        
        if score > 2*(float(self.DataModel.myChain.residueTypeFrequencyDict[ccpCode])/len(self.DataModel.myChain.residues)) :
          
          scoreDict[ccpCode] = score
          
    self.aminoAcidProbs = scoreDict
    
  cdef void setupAllowedResidueView(self) :
    
    cdef myChain chain
    
    cdef int resNumber
    
    chain = self.DataModel.myChain
    
    cythonArray = cvarray(shape=(len(chain.residues) + 1,), itemsize=sizeof(bint), format="i")
    
    self.allowedResidueView = cythonArray
    
    self.allowedResidueView[:] = False
    
    for resNumber in self.allowedResidues :
      
      self.allowedResidueView[resNumber] = True
      
  cdef myResonance getResonanceForAtomName(self,str atomName) :       # Used in matchSpectrum()
    
    if atomName in self.resonanceDict :
                
      return self.resonanceDict[atomName]
    
    else :
      
      return None
    
  cdef list getResonancesForAtomSiteName(self, str atomSiteName) :
    
    return self.resonancesByAtomSiteName.get(atomSiteName, [])

  cdef void setupAllowedResidues(self, bint useAssignments, bint useTentative) :
    
    cdef dict residuesByCcpCode
    cdef str ccpCode
    cdef aResidue res
    residuesByCcpCode = self.DataModel.myChain.residuesByCcpCode

    if useAssignments and self.ccpnSeqCode :

      self.allowedResidues = set([self.ccpnSeqCode])
    
    elif useTentative and self.tentativeSeqCodes :

      self.allowedResidues |= set(self.tentativeSeqCodes)
        
    elif self.ccpCode :

      self.allowedResidues |= set([res.seqCode for res in residuesByCcpCode[self.ccpCode]])
      
    elif self.tentativeCcpCodes :                        # Only end up in here when not useTentative

      for ccpCode in self.tentativeCcpCodes :

        self.allowedResidues |= set([res.seqCode for res in residuesByCcpCode[ccpCode]])
        
    elif self.aminoAcidProbs :

      for ccpCode in self.aminoAcidProbs.keys() :

        if ccpCode in residuesByCcpCode :

          self.allowedResidues |= set([res.seqCode for res in residuesByCcpCode[ccpCode]])
          
  cdef void setupExchangeSpinSystems(self, bint useAssignments) :
    
    cdef list spinSysList
    cdef mySpinSystem spinSys
    cdef dict spinSystemDict
    
    spinSystemDict = (useAssignments and self.DataModel.allSpinSystemsWithoutAssigned) or self.DataModel.mySpinSystems
    
    for spinSysList in spinSystemDict.values() :
      
      for spinSys in spinSysList :
        
        if not spinSys is self and not ( self.isJoker and spinSys.isJoker ) and  self.allowedResidues & spinSys.allowedResidues :
          
          self.exchangeSpinSystems.append(spinSys)

  cdef void groupResonancesByAtomSite(self) :                                 #TODO: finish
    
    cdef myResonance resonance
    
    resonancesByAtomSiteName = self.resonancesByAtomSiteName
    resonancesByName = self.resonanceDict
    
    ccpCodes = self.getCcpCodes()
    
    HAs = []
    HBs = []
    CXs = []
    Calis = []
    Caros = []
    COs = []
    
    resonancesByAtomSiteName['CA'] = []
    resonancesByAtomSiteName['CO'] = []
    
    
    # C,F,Br,H,P,N, CX
    for name, resonance in resonancesByName.items() :
      
      elementSymbol = resonance.ccpnResonance.isotope.chemElement.symbol
      
      resonancesByAtomSiteName[elementSymbol] = resonancesByAtomSiteName.get(elementSymbol, []) + [resonance]
      
      if elementSymbol == 'C' :
        
        CXs.append(resonance)
    
    # CA
    if 'CA' in resonancesByName :
  
      resonancesByAtomSiteName['CA'] = [resonancesByName['CA']]
      
    # CB
    if 'CB' in resonancesByName :
  
      resonancesByAtomSiteName['CB'] = [resonancesByName['CB']]
      
    # CO  
    if 'C' in resonancesByName :
  
      COs.append(resonancesByName['C'])
      
    if ('Gln' in ccpCodes or 'Glu' in ccpCodes) and 'CD' in resonancesByName :
      
      COs.append(resonancesByName['CD'])
      
    if ('Asp' in ccpCodes or 'Asn' in ccpCodes) and 'CG' in resonancesByName :
      
      COs.append(resonancesByName['CG'])
      
    # HA  
    for name in ['HA','HA1','HA2','HA3'] :
      
      if name in resonancesByName :
      
        HAs.append(resonancesByName[name])
    
    # HB    
    for name in ['HB','HB1','HB2','HB3', 'HB*'] :
      
      if name in resonancesByName :
      
        HBs.append(resonancesByName[name])
        
        
    # Caros
    if 'Phe' in ccpCodes or 'Tyr' in ccpCodes or 'Trp' in ccpCodes or 'His' in ccpCodes :
          
      resonancesByAtomSiteName['Caro'] = []
      
      for name in ['CG','CD1','CD2','CE1','CE2','CE3','CZ','CZ2','CZ3','CH2','CE*','CD*'] :
        
        if name in resonancesByName :
          
          resonance = resonancesByName[name]
          
          if 100 < resonance.CS < 170 :
          
            Caros.append(resonancesByName[name])
          
    #Calis
    else :
      
      for resonance in resonancesByAtomSiteName['C'] :
        
        if 0 < resonance.CS < 100 :
        #if resonance.atomName != 'C' :
          Calis.append(resonance)
          
    resonancesByAtomSiteName['HA'] = HAs
    resonancesByAtomSiteName['HB'] = HBs
    resonancesByAtomSiteName['CX'] = CXs
    resonancesByAtomSiteName['Cali'] = Calis
    resonancesByAtomSiteName['Caro'] = Caros
    resonancesByAtomSiteName['CO'] = COs
    
  cdef void createPythonStyleObject(self) :
    
    cdef aResidue res
    
    cdef myResonance resonance
    
    cdef str name

    self.pySpinSystem = pySpinSystem()

    self.pySpinSystem.spinSystemNumber = self.spinSystemNumber
    
    self.pySpinSystem.ccpCode = self.ccpCode
  
    self.pySpinSystem.ccpnSeqCode = self.ccpnSeqCode
    
    self.pySpinSystem.isJoker = self.isJoker
    
    self.pySpinSystem.tentativeCcpCodes = self.tentativeCcpCodes
  
    self.pySpinSystem.tentativeSeqCodes = self.tentativeSeqCodes
    
    self.pySpinSystem.allowedResidues = self.allowedResidues

    
    
    for res in self.solutions :
      
      self.pySpinSystem.solutions.append(res.pyResidue)
      
    for name,  resonance in self.resonanceDict.items() :
    
      self.pySpinSystem.resonanceDict[name] = resonance.pyResonance 