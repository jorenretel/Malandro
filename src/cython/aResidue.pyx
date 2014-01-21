
cdef class Residue :
  '''Residue, is used as a wrapper around the ccpn residue object.
  
    ccpnResidue:        the ccpnResidue that is wrapped.
    chain:              the molecular chain this residue belongs to.
    atoms:              list of atoms inside this residue.
    seqCode:            int descibing position in molecular chain.
    ccpCode:            three-letter code indicating the amino acid type
                        of the residue, 'Asp' for instance.
    previousResidue:    n-terminal neighbor of residue.
    nextResidue:        c-terminal neighbor of residue.
    linkDict:           dict that has the link objects that describes
                        the evidence that supports the existence of a
                        sequential link between every combination
                        of two spin-systems that could be assigned to
                        this residue and its c-terminal neighbor.
                        Keys to this dict are a combination of the serials
                        of these two spin-systems. The reason that this
                        'matrix' is an instance variable instead of having one
                        global matrix somewhere else is because, depending on
                        the amino acid types and labelling of the two residues,
                        the sequential peak-pattern differs and therefor also
                        the matching between expected and real peaks.
    intraDict:          same as linkDict, only here a link objects are stored
                        that are not really links between two spin-systems but
                        rather holds important information like the
                        intra-residual peaks for each spin-system
                        that could be assigned to this residue.
    atomsByName:        dict of {atomName:atom} in this residue.
    atomsByAtomSiteName: dict of {atomSite.name:atom}, where atomSite
                        is used in the graphs that descibe nmr experiments
                        in the ccpn data model.
    currentSpinSystemAssigned: can be used to store the current spin-system
                        assignment of the residue during the monte-carlo run.
    userDefinedSolution: meant to be used temporarely to store a spin-system
                        assignment that is made by the user while puzzling.
    solutions:          Usually multiple runs are performed. Here the different
                        spinSystems that were assigned in the different runs can
                        be stored.
                        
  '''

  cdef Residue previousResidue, nextResidue
  cdef public object ccpnResidue
  cdef public str ccpCode
  cdef public SpinSystem currentSpinSystemAssigned, userDefinedSolution
  cdef Chain chain
  cdef public list solutions
  cdef list atoms
  cdef dict atomsByAtomSiteName, linkDict, intraDict, atomsByName  #, atomsByCcpnChemAtom
  cdef int seqCode
  
  def __init__(self, chain, ccpnResidue):
    '''Init residue (wrapper). Also triggers the setup of the atoms
       within the residue.
       kwargs:  chain:        the chain the residue belongs to.
                              Not the ccpn chain object but the
                              chain object used in this macro.
                ccpnResidue:  the ccpn residue that is wrapped.

    '''

    self.chain = chain    # Parent link

    self.atoms = []
    
    self.solutions = []
    
    self.linkDict = {}
    
    self.intraDict = {}
    
    self.atomsByName = {}
    
    #self.atomsByCcpnChemAtom = {}  #Not used
    
    self.atomsByAtomSiteName = {}
    
    self.userDefinedSolution = None
    
    if ccpnResidue :
    
      self.ccpnResidue = ccpnResidue
  
      self.ccpCode = ccpnResidue.ccpCode
      
      self.seqCode = ccpnResidue.seqCode
      
      self.setupAtoms()

  def setupAtoms(self):
    '''Creates all the atom instances within this residue.
    
    '''
    for atom in self.ccpnResidue.atoms :
      
      newatom = Atom(self,atom)
      self.atoms.append(newatom)
      self.atomsByName[atom.chemAtom.name] = newatom
      #self.atomsByCcpnChemAtom[atom.chemAtom] = newatom       #Not used
      
    self.groupAtomsByAtomSite()
      
  cdef void groupAtomsByAtomSite(self) :
    '''Sort all atoms in the residue by atomSites. These atomSites
       are used in the desciption of experiments in the ccpn
       data model. The atomSite is in principle just a rough indicator
       of a specific group of atoms that can be pulsed on or measured
       in a dimension of an nmr experiment.
       A dict is generated and stored under the instance
       variable atomsByAtomSiteName: {atomSite.name : [atom,...]}.
       This dict is later on used to combine the experiment graph
       with the actual atoms in the molecule.
    
    '''
    cdef Atom atom
    
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
    '''Get the atoms in this residue that correspond to a
       ccpn atomSite.
       kwargs:  atomSite: the ccpn atomSite that is used to
                          describe a group of atoms in the
                          experiment graph.
       returns:           list of atoms.
    '''
    
    cdef Atom atom
    atomSiteName = atomSite.name
    
    return self.atomsByAtomSiteName.get(atomSiteName, [])
  
  #Not used, should maybe be removed.     
  cdef void addToLinkDict(self,SpinSystem spinSys1, SpinSystem spinSys2, list realPeaks, list simulatedPeaks, list notFoundSimulatedPeaks, list scores) :
    '''Generate a new spinSystemLink, if not present yet, and add
       it to the linkDict.
       kwargs:  spinSys1:       spinSystem that could be assigned
                                to this residue.
                spinSys2:       spinSystem that could be assigned
                                to the c-terminal neighboring residue.
                realPeaks:      list of peaks that have been found in
                                the spectra.
                simulatedPeaks: simulated peaks that correspond to
                                realPeaks.
                notFoundSimulatedPeaks: all other simulated peaks that
                                where not found in the spectra.
                scores:         list with scores that indicate how well
                                a peak fits the simulated peak. I.e. how
                                well do the peaks dimensions fit the
                                frequencies of the resonances thought to
                                give rise to them.
       N.B. This function is not in use in the rest of the code.
       Everything works directly with link.addPeak.
       
    '''
    cdef SpinSystemLink linkObject
    
    if (spinSys1.spinSystemNumber*10000+spinSys2.spinSystemNumber) in self.linkDict :

      linkObject = self.linkDict[(spinSys1.spinSystemNumber*10000+spinSys2.spinSystemNumber)]
      
    else :
      
      linkObject = SpinSystemLink()
      linkObject.spinSystem1 = spinSys1
      linkObject.spinSystem2 = spinSys2
      self.linkDict[(spinSys1.spinSystemNumber*10000+spinSys2.spinSystemNumber)] = linkObject
      
    for realPeak, simPeak, score in zip(realPeaks,simulatedPeaks,scores) :
      
      newPeakLink = PeakLink(realPeak,simPeak,score)
      linkObject.peakLinks.append(newPeakLink)

    linkObject.simulatedPeaks.extend(simulatedPeaks)
    linkObject.realPeaks.extend(realPeaks)
    linkObject.notFoundSimulatedPeaks.extend(notFoundSimulatedPeaks)
  
  #Not used, should maybe be removed.      
  cdef void addToIntraDict(self,SpinSystem spinSys, list realPeaks, list simulatedPeaks, list notFoundSimulatedPeaks, list scores) :
    '''Generate a new spinSystemLink, if not present yet, for
       intra-residual peaks and add it to the intraDict.
       kwargs:  spinSys:        spinSystem that could be assigned
                                to this residue.
                realPeaks:      list of peaks that have been found in
                                the spectra.
                simulatedPeaks: simulated peaks that correspond to
                                realPeaks.
                notFoundSimulatedPeaks: all other simulated peaks that
                                where not found in the spectra.
                scores:         list with scores that indicate how well
                                a peak fits the simulated peak. I.e. how
                                well do the peaks dimensions fit the
                                frequencies of the resonances thought to
                                give rise to them.
       N.B. This function is not in use in the rest of the code.
       Everything works directly with link.addPeak.
       
    '''
    
    cdef SpinSystemLink linkObject
    
    if spinSys.spinSystemNumber in self.intraDict :

      linkObject = self.intraDict[spinSys.spinSystemNumber]
      
    else :
      
      linkObject = SpinSystemLink()
      linkObject.spinSystem1 = spinSys
      linkObject.spinSystem2 = spinSys
      self.intraDict[spinSys.spinSystemNumber] = linkObject

    for realPeak, simPeak, score in zip(realPeaks,simulatedPeaks,scores) :
      
      newPeakLink = PeakLink(realPeak,simPeak,score)
      linkObject.peakLinks.append(newPeakLink)
      
    linkObject.simulatedPeaks.extend(simulatedPeaks)
    linkObject.realPeaks.extend(realPeaks)
    linkObject.notFoundSimulatedPeaks.extend(notFoundSimulatedPeaks)
    
  #Not used, should maybe be removed.
  cdef addPeakToLinkDict(self,SpinSystem spinSys1, SpinSystem spinSys2,Peak realPeak, SimulatedPeak simPeak, list resonances, double score, bint isIntra) :
    '''Add one peak to the linkDict'''
    
    cdef SpinSystemLink linkObject
    
    if isIntra is True:                   #Not correct, in the way the spectra are matched now, it might by accident happen to for instance an Ala-Ala pair
      
      if spinSys1.spinSystemNumber in self.intraDict :
  
        linkObject = self.intraDict[spinSys1.spinSystemNumber]
        
      else :
        
        linkObject = SpinSystemLink()
        linkObject.residue1 = self
        linkObject.residue2 = self
        linkObject.spinSystem1 = spinSys1
        linkObject.spinSystem2 = spinSys2
        self.intraDict[spinSys1.spinSystemNumber] = linkObject
    
    else :
      
      if (spinSys1.spinSystemNumber*10000+spinSys2.spinSystemNumber) in self.linkDict :
  
        linkObject = self.linkDict[(spinSys1.spinSystemNumber*10000+spinSys2.spinSystemNumber)]
        
      else :
        
        linkObject = SpinSystemLink()
        
        linkObject.residue1 = self
        linkObject.residue2 = self.nextResidue
        linkObject.spinSystem1 = spinSys1
        linkObject.spinSystem2 = spinSys2
        self.linkDict[(spinSys1.spinSystemNumber*10000+spinSys2.spinSystemNumber)] = linkObject
      
     
    if realPeak is not None :
      
      newPeakLink = PeakLink(realPeak,simPeak,resonances, score) 
      linkObject.peakLinks.append(newPeakLink)
      
    else :
      
      linkObject.notFoundSimulatedPeaks.append(simPeak)
      #linkObject.notFoundPeakLinks.append(newPeakLink)
       
  cdef SpinSystemLink getFromLinkDict(self, SpinSystem spinSystem1, SpinSystem spinSystem2) :
    '''Get a link object from the linkDict for this residue,
       only for sequential links.
       kwargs:  spinSystem1:  spinSystem corresponding to this residue
                spinSystem2:  spinSystem corresponding to c-terminal
                              neighboring residue.
       returns: link object describing the evidence that spinSystem1 and
                spinSystem2 are sequential neighbors that can be assigned
                to this residue and its c-terminal neighboring residue.
                
    '''
    #cdef bool joker1
    #cdef bool joker2
    
    #joker1 = spinSystem1.isJoker
    #joker2 = spinSystem2.isJoker
    
    if spinSystem1.isJoker or spinSystem2.isJoker:
      
      return emptyLink

    cdef SpinSystemLink link
    
    cdef int hashCode
    
    hashCode = spinSystem1.spinSystemNumber*10000+spinSystem2.spinSystemNumber
    
    link = <SpinSystemLink>self.linkDict[hashCode]
    
    return link #link.realPeaks, link.score
      
  def __reduce__(self) :

    return (generalFactory, (type(self),), self.__getstate__())
    
  def __getstate__(self) :

    return (self.chain, self.seqCode, self.ccpCode, self.atoms, self.solutions, self.linkDict, self.intraDict, self.previousResidue, self.nextResidue, self.userDefinedSolution)
  
  def __setstate__(self, state) :

    self.chain, self.seqCode, self.ccpCode, self.atoms, self.solutions, self.linkDict, self.intraDict, self.previousResidue, self.nextResidue, userDefinedSolution = state

  def connectToProject(self) :
    '''After unpycling this method can be ran to set the
       connection between this residue and its equivalent
       in the ccpn data model.
    
    '''
    
    cdef Atom atom
    
    self.ccpnResidue = self.chain.ccpnChain.findFirstResidue(seqCode=self.seqCode)
    
    if not self.ccpnResidue :
      
      print 'Error: It seems like you changed the length of the chain since you ran the calculations. I can not find residue %s .' %str(self.seqCode)
      return
    
    elif not self.ccpCode == self.ccpnResidue.ccpCode :
      
      print 'Error: It seems like you changed the residue type of residue %s from %s to %s since you ran the calculations' %(str(self.seqCode), self.ccpCode, self.ccpnResidue.ccpCode)
      
    for atom in self.atoms :
      
      atom.connectToProject()

  def getSolutions(self) :
    '''Returns spinSystems assigned in different runs.'''
    return self.solutions
  
  def getLink(self,spinSystem1,spinSystem2) :
    '''Returns link between spinSystem1 and spinSystem2.
       Same as getFromLinkDict but then accessable from
       python.
       
    '''
    
    return self.getFromLinkDict(spinSystem1,spinSystem2)
  
  def getIntraLink(self, SpinSystem spinSystem) :
    '''Returns a link object that contains the mapping
       between expected and real intra-residual peaks
       for 'spinSystem' being assigned to this residue.
       
    '''
    
    if spinSystem.isJoker :
      
      return emptyLink

    cdef SpinSystemLink link
    
    cdef int hashCode
    
    hashCode = spinSystem.spinSystemNumber
    
    link = self.intraDict[hashCode]
    
    return link

  def getSeqCode(self):
    '''Returns the chain-position identifier'''
    
    return self.seqCode
  
  def getCcpCode(self):
    '''Returns three-letter of amino acid type.'''
    return self.ccpCode
  
  def getCcpnResidue(self) :
    '''Returns residue from ccpn data model corresponding
       to this residue.
    '''
    
    return self.ccpnResidue
  
  def getPreviousResidue(self) :
    '''Returns n-terminal neighboring residue.'''
    
    return self.previousResidue
  
  def getNextResidue(self) :
    '''Returns n-terminal neighboring residue.'''
    
    return self.nextResidue

