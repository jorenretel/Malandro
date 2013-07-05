
cdef class aResidue :

  cdef aResidue previousResidue, nextResidue
  
  cdef public object ccpnResidue
  
  cdef public str ccpCode
  
  cdef public mySpinSystem userDefinedSolution
  
  cdef myChain chain
  
  cdef public list solutions
  
  cdef list atoms
  
  cdef dict atomsByName, atomsByAtomSiteName, atomsByCcpnChemAtom, linkDict, intraDict
  
  cdef mySpinSystem currentSpinSystemAssigned
  
  cdef int seqCode
  
  def __init__(self, chain, ccpnResidue):


    self.chain = chain    # Parent link

    self.atoms = []
    
    self.solutions = []
    
    self.linkDict = {}
    
    self.intraDict = {}
    
    self.atomsByName = {}
    
    self.atomsByCcpnChemAtom = {}
    
    self.atomsByAtomSiteName = {}
    
    self.userDefinedSolution = None
    
    if ccpnResidue :
    
      self.ccpnResidue = ccpnResidue
  
      self.ccpCode = ccpnResidue.ccpCode
      
      self.seqCode = ccpnResidue.seqCode
      
      self.setupAtoms()

  def setupAtoms(self):
    
    for atom in self.ccpnResidue.atoms :
      
      newatom = anAtom(self,atom)
      self.atoms.append(newatom)
      self.atomsByName[atom.chemAtom.name] = newatom
      self.atomsByCcpnChemAtom[atom.chemAtom] = newatom
      
    self.groupAtomsByAtomSite()
      
  cdef void groupAtomsByAtomSite(self) :
    
    cdef anAtom atom
    
    atomsByAtomSiteName = self.atomsByAtomSiteName
    atomsByName = self.atomsByName

    atomsByAtomSiteName['CA'] = [atomsByName['CA']]
    atomsByAtomSiteName['CO'] = [atomsByName['C']]
    
    HAs = []
    HBs = []
    CXs = []
    Calis = []
    
    for atom in self.atoms :
      
      elementSymbol = atom.ccpnAtom.chemAtom.elementSymbol
      
      atomsByAtomSiteName[elementSymbol] = atomsByAtomSiteName.get(elementSymbol, []) + [atom]
      
      if elementSymbol == 'C' :
        
        CXs.append(atom)
    
    if 'CB' in atomsByName :

      atomsByAtomSiteName['CB'] = [atomsByName['CB']]
      
    for name in ['HA','HA1','HA2','HA3'] :
      
      if name in atomsByName :
      
        HAs.append(atomsByName[name])
        
    for name in ['HB','HB1','HB2','HB3'] :
      
      if name in atomsByName :
      
        HBs.append(atomsByName[name])    
    
    if self.ccpCode == 'Phe' or self.ccpCode == 'Tyr' :
      
      Calis = [atomsByName[name] for name in ['CA','CB']]
      atomsByAtomSiteName['Caro'] = [atomsByName[name] for name in ['CG','CD1','CD2','CE1','CE2','CZ']]
      
    elif self.ccpCode == 'Trp' :
      
      Calis = [atomsByName[name] for name in ['CA','CB']]
      atomsByAtomSiteName['Caro'] = [atomsByName[name] for name in ['CG','CD1','CD2','CE2','CE3','CZ2', 'CZ3','CH2']]
      
    elif self.ccpCode == 'His' :
      
      Calis = [atomsByName[name] for name in ['CA','CB']]
      atomsByAtomSiteName['Caro'] = [atomsByName[name] for name in ['CG','CD2','CE1']]
      
    elif self.ccpCode == 'Glu' or self.ccpCode == 'Gln' :
      
      Calis = [atomsByName[name] for name in ['CA','CB','CG']]
      atomsByAtomSiteName['CO'].append(atomsByName['CD'])
      
    elif self.ccpCode == 'Asp' or self.ccpCode == 'Asn' :
      
      Calis = [atomsByName[name] for name in ['CA','CB']]
      atomsByAtomSiteName['CO'].append(atomsByName['CG'])
      
    else :
      
      for atom in atomsByAtomSiteName['C'] :
        
        if atom.atomName != 'C' :
          
          Calis.append(atom)
          
    atomsByAtomSiteName['HA'] = HAs
    atomsByAtomSiteName['HB'] = HBs
    atomsByAtomSiteName['CX'] = CXs
    atomsByAtomSiteName['Cali'] = Calis
          
  cdef list getAtomsForAtomSite(self, object atomSite) :
    
    cdef anAtom atom
    atomSiteName = atomSite.name
    
    return self.atomsByAtomSiteName.get(atomSiteName, [])
      
  cdef void addToLinkDict(self,mySpinSystem spinSys1, mySpinSystem spinSys2, list realPeaks, list simulatedPeaks, list notFoundSimulatedPeaks, list scores) :
  
    cdef spinSystemLink linkObject
    
    if (spinSys1.spinSystemNumber*10000+spinSys2.spinSystemNumber) in self.linkDict :

      linkObject = self.linkDict[(spinSys1.spinSystemNumber*10000+spinSys2.spinSystemNumber)]
      
    else :
      
      linkObject = spinSystemLink()
      linkObject.spinSystem1 = spinSys1
      linkObject.spinSystem2 = spinSys2
      self.linkDict[(spinSys1.spinSystemNumber*10000+spinSys2.spinSystemNumber)] = linkObject
      
    for realPeak, simPeak, score in zip(realPeaks,simulatedPeaks,scores) :
      
      newPeakLink = peakLink(realPeak,simPeak,score)
      linkObject.peakLinks.append(newPeakLink)

    linkObject.simulatedPeaks.extend(simulatedPeaks)
    linkObject.realPeaks.extend(realPeaks)
    linkObject.notFoundSimulatedPeaks.extend(notFoundSimulatedPeaks)
        
  cdef void addToIntraDict(self,mySpinSystem spinSys, list realPeaks, list simulatedPeaks, list notFoundSimulatedPeaks, list scores) :
  
    cdef spinSystemLink linkObject
    
    if spinSys.spinSystemNumber in self.intraDict :

      linkObject = self.intraDict[spinSys.spinSystemNumber]
      
    else :
      
      linkObject = spinSystemLink()
      linkObject.spinSystem1 = spinSys
      linkObject.spinSystem2 = spinSys
      self.intraDict[spinSys.spinSystemNumber] = linkObject

    for realPeak, simPeak, score in zip(realPeaks,simulatedPeaks,scores) :
      
      newPeakLink = peakLink(realPeak,simPeak,score)
      linkObject.peakLinks.append(newPeakLink)
      
    linkObject.simulatedPeaks.extend(simulatedPeaks)
    linkObject.realPeaks.extend(realPeaks)
    linkObject.notFoundSimulatedPeaks.extend(notFoundSimulatedPeaks)
    
  cdef addPeakToLinkDict(self,mySpinSystem spinSys1, mySpinSystem spinSys2, aPeak realPeak, simulatedPeak simPeak, list resonances, double score, bint isIntra) :
    
    cdef spinSystemLink linkObject
    
    if isIntra is True:                   #Not correct, in the way the spectra are matched now, it might by accident happen to for instance an Ala-Ala pair
      
      if spinSys1.spinSystemNumber in self.intraDict :
  
        linkObject = self.intraDict[spinSys1.spinSystemNumber]
        
      else :
        
        linkObject = spinSystemLink()
        linkObject.spinSystem1 = spinSys1
        linkObject.spinSystem2 = spinSys2
        self.intraDict[spinSys1.spinSystemNumber] = linkObject
    
    else :
      
      if (spinSys1.spinSystemNumber*10000+spinSys2.spinSystemNumber) in self.linkDict :
  
        linkObject = self.linkDict[(spinSys1.spinSystemNumber*10000+spinSys2.spinSystemNumber)]
        
      else :
        
        linkObject = spinSystemLink()
        linkObject.spinSystem1 = spinSys1
        linkObject.spinSystem2 = spinSys2
        self.linkDict[(spinSys1.spinSystemNumber*10000+spinSys2.spinSystemNumber)] = linkObject
      
    newPeakLink = peakLink(realPeak,simPeak,resonances, score)  
      
    if realPeak is not None :

      linkObject.peakLinks.append(newPeakLink)
      
    else :
      
      linkObject.notFoundPeakLinks.append(newPeakLink)
      
  cdef spinSystemLink getFromLinkDict(self, mySpinSystem spinSystem1, mySpinSystem spinSystem2) :
    
    #cdef bool joker1
    #cdef bool joker2
    
    #joker1 = spinSystem1.isJoker
    #joker2 = spinSystem2.isJoker
    
    if spinSystem1.isJoker or spinSystem2.isJoker:
      
      return emptyLink

    cdef spinSystemLink link
    
    cdef int hashCode
    
    hashCode = spinSystem1.spinSystemNumber*10000+spinSystem2.spinSystemNumber
    
    link = <spinSystemLink>self.linkDict[hashCode]
    
    return link #link.realPeaks, link.score
      
  def __reduce__(self) :

    return (generalFactory, (type(self),), self.__getstate__())
    
  def __getstate__(self) :

    return (self.chain, self.seqCode, self.ccpCode, self.atoms, self.solutions, self.linkDict, self.intraDict, self.previousResidue, self.nextResidue, self.userDefinedSolution)
  
  def __setstate__(self, state) :

    self.chain, self.seqCode, self.ccpCode, self.atoms, self.solutions, self.linkDict, self.intraDict, self.previousResidue, self.nextResidue, userDefinedSolution = state

  def connectToProject(self) :
    
    cdef anAtom atom
    
    self.ccpnResidue = self.chain.ccpnChain.findFirstResidue(seqCode=self.seqCode)
    
    if not self.ccpnResidue :
      
      print 'Error: It seems like you changed the length of the chain since you ran the calculations. I can not find residue %s .' %str(self.seqCode)
      return
    
    elif not self.ccpCode == self.ccpnResidue.ccpCode :
      
      print 'Error: It seems like you changed the residue type of residue %s from %s to %s since you ran the calculations' %(str(self.seqCode), self.ccpCode, self.ccpnResidue.ccpCode)
      
    for atom in self.atoms :
      
      atom.connectToProject()

  def getSolutions(self) :
    
    return self.solutions
  
  def getLink(self,spinSystem1,spinSystem2) :
    
    return self.getFromLinkDict(spinSystem1,spinSystem2)
  
  def getIntraLink(self, mySpinSystem spinSystem) :
    
    if spinSystem.isJoker :
      
      return emptyLink

    cdef spinSystemLink link
    
    cdef int hashCode
    
    hashCode = spinSystem.spinSystemNumber
    
    link = self.intraDict[hashCode]
    
    return link

  def getSeqCode(self):
    
    return self.seqCode
  
  def getCcpCode(self):
    
    return self.ccpCode
