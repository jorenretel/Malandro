
cdef class aResidue :

  cdef aResidue previousResidue, nextResidue
  
  cdef object ccpnResidue, pyResidue
  
  cdef str ccpCode
  
  cdef myChain chain
  
  cdef list atoms, solutions
  
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
    
    if ccpnResidue :
    
      self.ccpnResidue = ccpnResidue
  
      self.ccpCode = ccpnResidue.ccpCode
      
      self.seqCode = ccpnResidue.seqCode
      
      self.setupAtoms()
    
  def __getstate__(self):
    state = dict(self.__dict__)
    if 'ccpnResidue' in state :
      del state['ccpnResidue']
    return state

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
    
  cdef addPeakToLinkDict(self,mySpinSystem spinSys1, mySpinSystem spinSys2, aPeak realPeak, aPeak simulatedPeak, list resonances, double score) :
    
    cdef spinSystemLink linkObject
    
    if spinSys1 is spinSys2 :                   #Not correct, in the way the spectra are matched now, it might by accident happen to for instance an Ala-Ala pair
      
      if spinSys.spinSystemNumber in self.intraDict :
  
        linkObject = self.intraDict[spinSys.spinSystemNumber]
        
      else :
        
        linkObject = spinSystemLink()
        linkObject.spinSystem1 = spinSys
        linkObject.spinSystem2 = spinSys
        self.intraDict[spinSys.spinSystemNumber] = linkObject
    
    else :
      
      if (spinSys1.spinSystemNumber*10000+spinSys2.spinSystemNumber) in self.linkDict :
  
        linkObject = self.linkDict[(spinSys1.spinSystemNumber*10000+spinSys2.spinSystemNumber)]
        
      else :
        
        linkObject = spinSystemLink()
        linkObject.spinSystem1 = spinSys1
        linkObject.spinSystem2 = spinSys2
        self.linkDict[(spinSys1.spinSystemNumber*10000+spinSys2.spinSystemNumber)] = linkObject
      
    newPeakLink = peakLink(realPeak,simPeak,score)  
      
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
    
  cdef void createPythonStyleObject(self):
    
      self.pyResidue = pyResidue()
      
      self.pyResidue.seqCode = self.seqCode
      
      self.pyResidue.ccpCode = self.ccpCode