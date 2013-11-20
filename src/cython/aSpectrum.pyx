
cdef class Spectrum:
  
  cdef public str name
  
  cdef myDataModel DataModel
  
  cdef object ccpnSpectrum, ccpnPeakList, labellingScheme, pySpectrum
 
  cdef list simulatedPeakMatrix, intraResidualSimulatedPeakMatrix, peaks
  
  cdef int symmetry, experimentSerial, serial
  
  cdef frozenset molLabelFractions
  
  def __init__(self, temporary_spectrum_object):
    
    self.DataModel = None
    self.ccpnSpectrum = temporary_spectrum_object.ccpnSpectrum
    self.experimentSerial = self.ccpnSpectrum.experiment.serial
    self.serial = self.ccpnSpectrum.serial
    self.name = temporary_spectrum_object.ccpnSpectrum.name
    self.labellingScheme = temporary_spectrum_object.labellingScheme
    self.symmetry = 1
    
    if self.labellingScheme is True :
      
      self.molLabelFractions = self.ccpnSpectrum.experiment.findFirstLabeledMixture().getMolLabelFractions()
      
    else :
      
      self.molLabelFractions = None
    
    self.ccpnPeakList = temporary_spectrum_object.peakList
    
    self.simulatedPeakMatrix = []
    self.intraResidualSimulatedPeakMatrix = []

    self.peaks = []
    
    self.setupPeaks()
    
  def setupPeaks(self):

    peaks = self.ccpnPeakList.peaks

    for peak in peaks:
      
      newpeak =Peak(self,peak)

      self.peaks.append(newpeak)

  def match(self) :
    
    self.matchIntraResidual()
    self.matchSequential()

  cdef void simulate(self) :
    
    cdef myDataModel DataModel
    
    cdef list residues
    
    cdef Spectrum spectrum
    
    cdef object ccpnSpectrum
    
    cdef object scheme
    
    cdef list simulatedPeakMatrix
    
    cdef Residue resA
    
    cdef Residue resB
    
    cdef object refExperiment
    
    cdef object PT
    
    cdef object expGraph
    
    cdef dict recordedExpMeasuments
    
    cdef list expTransfers
    
    #cdef object refExpDim
    
    cdef list atomSitePathWay
    
    cdef object expStep
    
    cdef object atomSite
    
    cdef Atom atom
    
    cdef list atomPathWays
    
    cdef list atomPathWay
    
    cdef list peaks
    
    cdef list isotopeCodes
    
    cdef double minIsoFrac
    
    cdef SimulatedPeakContrib contrib
    
    cdef tuple dimNumbers
    cdef tuple stepNumbers
    cdef int dimNumber
    cdef int stepNumber
    
    DataModel = self.DataModel
    
    minIsoFrac = DataModel.minIsoFrac
    
    residues = DataModel.chain.residues 
    
    ccpnSpectrum = self.ccpnSpectrum
    
    scheme = self.labellingScheme
    
    simulatedPeakMatrix = self.simulatedPeakMatrix
    intraResidualSimulatedPeakMatrix = self.intraResidualSimulatedPeakMatrix
  
    refExperiment = ccpnSpectrum.experiment.refExperiment
    
    #molLabelFractions = self.molLabelFractions  #ccpnSpectrum.experiment.findFirstLabeledMixture().getMolLabelFractions()
    
    PT =  refExperiment.nmrExpPrototype
    
    # Starting off with walking over the first expGraph to see which steps are recorded in which dimensions and which isotopes are on the magnetization
    # transfer pathway. I do take projection spectra in account. I also asume experiments here where not two different isotopes are independently recorded
    # in the same step or where different expGraphs have a different way of mapping dimensions to steps in the pulse sequence (And I think if you were
    # doing this, you are are most likely processing your data in such a way that you end up with two spectra instead of one.
    
    expGraph = PT.findFirstExpGraph()    
    
    expTransfers = expGraph.sortedExpTransfers()
    
    #expSteps = sorted(expGraph.sortedExpSteps(), key=lambda x: x.stepNumber)


    
    expSteps = list(zip(*sorted([(expStep.stepNumber, expStep) for expStep in expGraph.sortedExpSteps()]))[1])
    
    outAndBack = self.isOutAndBack(expSteps)
    
    if outAndBack :

      expSteps = expSteps[:int(len(expSteps)/2.0+0.5)]
    
    atomSitePathWay = self.createAtomSitesList(expSteps)

    dimStepDict = self.mapExpStepToDimension(expSteps, refExperiment.sortedRefExpDims())
    dimNumbers, stepNumbers = zip(*sorted(dimStepDict.items()))
    isotopeCodes = [atomSite.isotopeCode for atomSite in atomSitePathWay]
    measuredIsotopeCodes = [isotopeCodes[stepNumber-1] for stepNumber in stepNumbers]
    
    
    #Now do get the visited atomSites and the transfers between these sites for each expGraph.
    
    expGraphs = PT.getExpGraphs()
    atomSiteAndTransferPathways = []
      
    for expGraph in expGraphs :
      
      expTransfers = expGraph.sortedExpTransfers()
      
      expSteps = list(zip(*sorted([(expStep.stepNumber, expStep) for expStep in expGraph.sortedExpSteps()]))[1])

        
      if outAndBack :

        expSteps = expSteps[:int(len(expSteps)/2.0+0.5)]

      # Making a list of atomSites that are visited.
      
      atomSitePathWay = self.createAtomSitesList(expSteps)
        
      # And a list with expTransfers that connect the atomSites
      
      transferPathWay = self.createTransferList(atomSitePathWay,expTransfers)
      
      atomSiteAndTransferPathways.append( (atomSitePathWay,transferPathWay) )
      
      

    # Going through each sequential pair in the sequence
    
    intraResidualAtomSetDictB = {}
    
    lastResidue = residues[-1]

    for resA, resB in zip(residues,residues[1:]) :
      
      isLast = resB is lastResidue
      
      atomPathWays = []
      
      for atomSitePathWay, transferPathWay in atomSiteAndTransferPathways :
      
        atomGroups = []
   
        for atomSite in atomSitePathWay :
          
          atomGroups.append(resA.getAtomsForAtomSite(atomSite) + resB.getAtomsForAtomSite(atomSite))

        atomPathWays.extend( self.walkExperimentTree([], atomGroups, transferPathWay,0) )
        

      self.cacheLabellingInfo(atomPathWays)

      
      # In the next few lines the atomPathways are grouped by atomSet. This takes care of two things :
      # 1. Different pathways that end up having the same atom(sets) in the measured dimensions (i.e give rise to the same peak) are grouped together.
      # 2. Also multiple pathways that in practice give rise to the same peak because the measured atom is in fast exchange (for instance the three protons attached to the CB of Ala)
      #    are grouped together, so only one peak is generated. The only time that this could actually cause a underestimation of the amount of peaks is in Phe and Tyr
      #    if the ring two sides of the ring are not equivalent. This can later on cause some peaks not being picked up in the matching procedure with the real spectra when
      #    the spin system used for matching has the ring set to 'non-equivalent' and therefor the CD1/CD2 and CE1/CE2 would be assigned to different resonances.
      # At the same time they get sorted into sequential peaks, intra-residual peaks from the first residue and second residue.
      sequentialAtomSetDict = {}
      intraResidualAtomSetDictA = intraResidualAtomSetDictB
      intraResidualAtomSetDictB = {}

      for atomPathWay in atomPathWays :
        
        atomSets = []
        measuredResidues = set()
        
        for stepNumber in stepNumbers :
          
          atom = atomPathWay[stepNumber-1]
          measuredResidues.add(atom.residue)
          atomSets.append(atom.ccpnAtom.atomSet)
        
        atomSetTuple = tuple(atomSets)
        
        # Sequential Peak.
        if len(measuredResidues) > 1 :

          sequentialAtomSetDict[atomSetTuple] = sequentialAtomSetDict.get(atomSetTuple,[]) + [atomPathWay]
        
        # All dimensional contributions to peak are from the first residue.  
        elif measuredResidues.pop() is resA :
            
          intraResidualAtomSetDictA[atomSetTuple] = intraResidualAtomSetDictA.get(atomSetTuple,[]) + [atomPathWay]
        
        # All dimensional contributions to peak are from the second residue
        # but the magnetization transferpathway includes atoms from the first residue.
        elif len(set([atom.residue for atom in atomPathWay])) > 1 or isLast :
            
          intraResidualAtomSetDictB[atomSetTuple] = intraResidualAtomSetDictB.get(atomSetTuple,[]) + [atomPathWay]  
            
            
        
      # Create a new simulatedPeak if the colabelling is sufficient.

      intraResidualPeaksA = []
      sequentialPeaks = []
      
      peaksToMake = [(sequentialAtomSetDict,sequentialPeaks),(intraResidualAtomSetDictA,intraResidualPeaksA)]
      
      if isLast :
        intraResidualPeaksB = []
        peaksToMake = [(sequentialAtomSetDict,sequentialPeaks),(intraResidualAtomSetDictA,intraResidualPeaksA),(intraResidualAtomSetDictB,intraResidualPeaksB)]
      
      for atomSetDict, peakList in peaksToMake :
       
        for atomSetTuple, atomPathWayList in atomSetDict.items() :

          colabelling = 0.0
          
          for atomPathWay in atomPathWayList :

            colabelling += self.getCoLabellingFractionNew(atomPathWay, isotopeCodes)
            
          colabelling = colabelling / len(atomPathWayList)          # Average colabelling over all pathways giving rise to the same peak.

          if not colabelling > minIsoFrac :
            
            continue
            
          firstAtomPathWay = atomPathWayList[0]

          newPeak = SimulatedPeak()
          newPeak.colabelling = colabelling
          newPeak.spectrum = self

          for dimNumber, stepNumber, atomSet, isotopeCode in zip(dimNumbers, stepNumbers, atomSetTuple, measuredIsotopeCodes) :

            contrib = SimulatedPeakContrib()

            atom = firstAtomPathWay[stepNumber-1]

            residue = atom.residue
            
            contrib.residue = residue
            
            contrib.ccpCode = residue.ccpCode
            
            contrib.atomName = atomSet.name
            
            contrib.isotopeCode = isotopeCode #isotopeCodes[stepNumber-1]
            
            if atom.residue is resA :
            
              contrib.firstOrSecondResidue = 1
              
            elif atom.residue is resB :
            
              contrib.firstOrSecondResidue = 2  
        
            contrib.dimNumber = dimNumber

            newPeak.simulatedPeakContribs.append(contrib)
            
          peakList.append(newPeak)
      
      simulatedPeakMatrix.append(sequentialPeaks)
      intraResidualSimulatedPeakMatrix.append(intraResidualPeaksA)
      if isLast :
        intraResidualSimulatedPeakMatrix.append(intraResidualPeaksB)
         
  cdef object transferIsPossible(self, Atom atomA, Atom atomB, object expTransfer) :
    
    if not expTransfer :
      
      if atomA is atomB :
        
        return True
      
      else :
        
        return False
    
    transferType = expTransfer.transferType
    
    if not expTransfer.transferToSelf and atomA is atomB :
      
      return False
    
    if (transferType == 'onebond' or transferType == 'Jcoupling')  :
      
      # Over the amide bond
      #if (atomA.atomName == 'N' and atomB.atomName == 'C' and resIDA == 2 and resIDB == 1) or (atomA.atomName == 'C' and atomB.atomName == 'N' and resIDA == 1 and resIDB == 2) :
      if (atomA.atomName == 'N' and atomB.atomName == 'C' and atomA.residue is atomB.residue.nextResidue) or (atomA.atomName == 'C' and atomB.atomName == 'N' and atomB.residue is atomA.residue.nextResidue) :
    
        return True
        
      # Within the same residue and directly bonded
      elif atomA.residue is atomB.residue and atomA.ccpnAtom.chemAtom.getChemBonds().intersection(atomB.ccpnAtom.chemAtom.getChemBonds()):
        
        return True
        
      # N CAi-1 transfer
      elif transferType == 'Jcoupling' and (atomA.atomName == 'N' and atomB.atomName == 'CA' and atomA.residue is atomB.residue.nextResidue) or (atomA.atomName == 'CA' and atomB.atomName == 'N' and atomB.residue is atomA.residue.nextResidue) :
      
        return True
      
      else :
        
        return False
      
    else :
      
      return True
 
  cdef object isOutAndBack(self, list expSteps) :
    
    cdef list measurements
    
    measurements = [expStep.expMeasurement for expStep in expSteps]
    
    if measurements == measurements[::-1] :                     # Checking whether lists are symetric by turning it around
      
      return True
    
    else :
      
      return False
    
  cdef dict mapExpStepToDimension(self, list expSteps, list refExpDims) :
    
    dimStepDict = {}
    
    for refExpDim in refExpDims :
      
      i = 0
      
      for expStep in expSteps[::-1] :
        
        if expStep.expMeasurement is refExpDim.findFirstRefExpDimRef().expMeasurement :
          
          dimStepDict[refExpDim.dim] = expStep.stepNumber
          
          break
    
    return dimStepDict
  
  cdef dict createTransferDict(self, list expTransfers) :
    
    transferDict = {}
    
    for expTransfer in expTransfers :
      
      atomSites = expTransfer.sortedAtomSites()
      
      atomSite1, atomSite2 = atomSites[0], atomSites[1]
      
      transferDict[(atomSite1,atomSite2)] = expTransfer
      transferDict[(atomSite2,atomSite1)] = expTransfer
      
    return transferDict
  
  cdef list createTransferList(self, list atomSites, list expTransfers) :
    
    transferDict = self.createTransferDict(expTransfers)
    
    transferPathWay = []
    
    for atomSiteA, atomSiteB in zip(atomSites,atomSites[1:]) :
      
      if (atomSiteA,atomSiteB) in transferDict :

        transferPathWay.append( transferDict[(atomSiteA,atomSiteB)] )
        
      else :
        
        transferPathWay.append(None)
        
    return transferPathWay
  
  cdef list createAtomSitesList(self, list expSteps) :
    
    return [expStep.expMeasurement.findFirstAtomSite() for expStep in expSteps]

  cdef list walkExperimentTree(self,list pastAtoms, list atomGroups, list expTransfers, int depth) :
    
    if depth == len(atomGroups) :
      
      #if self.atomsBelongToSameResidue(pastAtoms) :
        
      #  return []
      
      return [pastAtoms]
    
    listOfLists = []
    
    expTransfer = expTransfers[depth-1]

    atomGroup = atomGroups[depth]

    depth = depth+1

    for atomB in atomGroup :

      if not ( pastAtoms and not self.transferIsPossible(pastAtoms[-1], atomB, expTransfer) ):

        extendedPathway = self.walkExperimentTree(pastAtoms + [atomB],atomGroups,expTransfers,depth)

        listOfLists += extendedPathway

    return listOfLists
      
  cdef void cacheLabellingInfo(self, list atomPathWays) :
    
    cdef Atom atom
    cdef set importantAtoms
    
    molLabelFractions = self.molLabelFractions
    
    importantAtoms = set()
    
    for atomPathWay in atomPathWays :   

      for atom in atomPathWay :
        
        importantAtoms.add(atom)
        
    for atom in importantAtoms :    
      
      atomName = atom.ccpnAtom.name
      
      subType = atom.ccpnAtom.chemAtom.subType
    
      molResidue = atom.residue.ccpnResidue.molResidue
        
      atom.labelInfoTemp = {}
          
      for mlf in molLabelFractions:                                                             
      
        molLabel = mlf.molLabel
        
        resLabel = molLabel.findFirstResLabel(resId=molResidue.serial)
        resLabelFractions = resLabel.resLabelFractions
        
        for rlf in resLabelFractions :
          
          isotopomers = rlf.isotopomers
          
          for isotopomer in isotopomers :
            
            atom.labelInfoTemp[isotopomer] = {}
            
            isotopeDict = atom.labelInfoTemp[isotopomer]
            
            atomLabels = isotopomer.findAllAtomLabels(name=atomName, subType=subType)
            
            atomLabelWeightSum = float(sum([atomLabel.weight for atomLabel in atomLabels]))
            
            for atomLabel in atomLabels :
              
              isotopeDict[atomLabel.isotopeCode] = atomLabel.weight / atomLabelWeightSum
  
  cdef double getCoLabellingFractionNew(self,list atoms, list isotopeCodes):
    
    cdef object scheme
    cdef object ccpnSpectrum
    cdef Residue x
    cdef list ccpnResidues
    
    cdef double colabelling
    
    
    scheme = self.labellingScheme
    
    if scheme is None :
      
      return 1.0
    
    elif scheme is True :
      
      colabelling = self.calculateCoLabellingFromExperimentNew(atoms, isotopeCodes)
      
      return colabelling
    
  cdef double calculateCoLabellingFromExperimentNew(self, list atoms, list isotopeCodes) :
    
    cdef object mixture
    
    cdef double TotalCoLabellingFraction
    
    cdef dict residueDict
    
    cdef Atom atom
    
    cdef double molWeightSum
    
    cdef double resFactor
    
    cdef object molResidue
    
    cdef object ccpnAtom
    
    cdef frozenset molLabelFractions
    
    cdef object mlf
    
    cdef double molFactor
    
    cdef object molLabel
    
    cdef double colabellingOfAtomsInThismolLabelFraction
    
    cdef double colabellingOfAtomsWithinThisResidue
    
    cdef list ccpnAtomList
    
    cdef str atomName
    
    cdef object chemAtom
    
    cdef dict fracDict
    
    cdef double fraction
    
    cdef list subTypes
    
    cdef list atomNameList
    
    cdef double addedColabelling
    
    molLabelFractions = self.molLabelFractions
  
    #mixture = self.ccpnSpectrum.experiment.findFirstLabeledMixture()
      
    #molLabelFractions = mixture.molLabelFractions
    
    molWeightSum = sum([x.weight for x in molLabelFractions])                                                               # This is the sum of weight of all labelling patterns present in the sample
    
    residueDict = {}
    
    for atom, isotopeCode in zip(atoms, isotopeCodes) :                                                             # Group Atoms by residue
      
      molResidue = atom.residue.ccpnResidue.molResidue

      if molResidue in residueDict :
        
        residueDict[molResidue].add((atom,isotopeCode))
        
      else :
      
        residueDict[molResidue] = set([(atom,isotopeCode)])
        
    TotalCoLabellingFraction = 0.0    
         
    for mlf in molLabelFractions:                                                                                                                       # Now I want to loop over all different labelling patterns that are present in the sample (denoted with capital letter in the analysis GUI)
      
      molFactor = mlf.weight / molWeightSum
        
      molLabel = mlf.molLabel
      
      colabellingOfAtomsInThismolLabelFraction= 1.0
      
      for molResidue, atomIsotopeTuples in residueDict.items() :                                                                                         # Loop over involved residues
        
        resLabel = molLabel.findFirstResLabel(resId=molResidue.serial)                                                                                  # TODO: do something about this.
        resLabelFractions = resLabel.resLabelFractions
        
        rlfWeightSum = sum([x.weight for x in resLabelFractions])
        
        colabellingOfAtomsWithinThisResidue = 0
        
        for rlf in resLabelFractions :                                                                                                                         # Loop over the different resLabelFractions, bascically looping over isotopomers of one residue in a specific labelling pattern
          
          resFactor = rlf.weight / rlfWeightSum                                                                                                           # The relative weight of the different isotopomers

          fraction =self.getIsotopomerAtomSetFractions(rlf.isotopomers, atomIsotopeTuples)

          addedColabelling =  fraction * resFactor

          colabellingOfAtomsWithinThisResidue += addedColabelling
      
        colabellingOfAtomsInThismolLabelFraction *= colabellingOfAtomsWithinThisResidue

      TotalCoLabellingFraction += ( colabellingOfAtomsInThismolLabelFraction * molFactor )

    return TotalCoLabellingFraction 
         
  cdef double getIsotopomerAtomSetFractions(self,set isotopomers, set atomIsotopeTuples) :
    
    cdef double isoWeightSum

    cdef object isotopomer
    
    cdef double labellingFraction
    
    cdef Atom atom
    
    cdef str isotopeCode
    
    labellingFraction = 0.0
    
    isoWeightSum =  float(sum([isotopomer.weight for isotopomer in isotopomers]))
    
    for isotopomer in isotopomers :
      
      colabellingInThisIsotopomer = 1.0
      
      for atom, isotopeCode in atomIsotopeTuples :
      
        colabellingInThisIsotopomer *= atom.labelInfoTemp[isotopomer][isotopeCode]
        
      labellingFraction += colabellingInThisIsotopomer * isotopomer.weight / isoWeightSum
    
    return labellingFraction

  cdef void findPossiblePeakContributions(self,useDimenionalAssignments=False) :
    
    cdef list resonances
    cdef Resonance resonance
    cdef Peak peak
    cdef PeakDimension dim
    cdef Peak firstPeak
    cdef object atomSite
    cdef dict dimAtomsDict
    cdef double ppmValue
    cdef double tolerance
    
    dimAtomsDict = {}

    if self.peaks:
      
      firstPeak = self.peaks[0]
      
      for dim in firstPeak.dimensions :
        
        dataDim = dim.ccpnDim.dataDim
        
        atomSite = dataDim.expDim.refExpDim.findFirstRefExpDimRef().expMeasurement.findFirstAtomSite()
        
        resonances = self.DataModel.getResonancesForAtomSiteName(atomSite.name)
                
        tolerance = getAnalysisDataDim(dataDim).assignTolerance
        
        dimAtomsDict[dim.dimNumber] = (resonances, tolerance)
        
    
    for peak in self.peaks :

      for dim in peak.dimensions :
        
        resonances, tolerance = dimAtomsDict[dim.dimNumber]
        
        ppmValue = dim.ppmValue
        
        assignedContributions = dim.ccpnDim.peakDimContribs
        
        if assignedContributions and ( useDimenionalAssignments is True or peak.intraResidual ) :
          
          assignedResonances = set([contrib.resonance for contrib in assignedContributions])
          
          for resonance in resonances :
            
            if resonance.ccpnResonance in assignedResonances :
              
              dim.possibleContributions.append(resonance)
              resonance.addPeakToPeakDimsLib(peak,dim)
              
              
        else :
          
          for resonance in resonances :
          
            if abs(resonance.CS - ppmValue) <= tolerance :
              
              dim.possibleContributions.append(resonance)
              resonance.addPeakToPeakDimsLib(peak,dim)

  cdef void matchSequential(self):
    
    cdef myDataModel DataModel
    
    cdef list simulatedPeakMatrix, residues, simulatedPeakList, spinsystemsaa1, spinsystemsaa2, listWithPresentPeaks, listWithSimulatedPeaks, listWithNotFoundSimulatedPeaks, listWithScores, contributions, resonances, peakLists, peaksInWindow

    cdef Residue resA, resB
    
    cdef dict allSpinSystems
    
    cdef SpinSystem spinSys1, spinSys2, spinSystem
    
    cdef int presentPeaks, i
    
    cdef double symmetry
    
    cdef SimulatedPeak simulatedPeak
    
    cdef SimulatedPeakContrib contrib
    
    cdef PeakDimension dim
    
    cdef Resonance resonance
    
    cdef Peak peak
    
    cdef int hc
    
    cdef double bestScore
    
    cdef SpinSystemLink link
    
    hc = 10000 #self.hc

    DataModel = self.DataModel
    simulatedPeakMatrix = self.simulatedPeakMatrix
    
    allSpinSystems = DataModel.spinSystems
    
    residues = DataModel.chain.residues
    
    symmetry = float(self.symmetry)     # preventing int division

    for i,  simulatedPeakList in enumerate(simulatedPeakMatrix) :
      
      resA = residues[i]
      resB = residues[i+1]
      
      linkDict = resA.linkDict
      
      for link in linkDict.values() :
        
        spinSys1 = link.spinSystem1
        spinSys2 = link.spinSystem2
          
        spinSystems = [spinSys1,spinSys2]

        for simulatedPeak in simulatedPeakList :
          
          closestPeak = None
          bestScore = 0.0
          resonances = []
          
          contributions = simulatedPeak.simulatedPeakContribs

          for contrib in contributions :
            
            spinSystem = spinSystems[contrib.firstOrSecondResidue-1]
            
            resonance = spinSystem.getResonanceForAtomName(contrib.atomName)
            
            resonances.append(resonance)
            
          if resonances and not None in resonances :
            
            peakLists = [ resonance.getPeaksForSpectrumDim(self,contrib.dimNumber) for resonance, contrib in zip(resonances , contributions) ]

            peaksInWindow = commonElementInLists(peakLists)                                      # And check whether 1 or more peaks that fit in one dimension actually also fit in all other dimensions. In that case the peak is in the multidimensional tolerance window

            if peaksInWindow :
              
              peakScores = [scorePeak(peak.dimensions,resonances) for peak in peaksInWindow]
              
              bestScore, closestPeak = sorted(zip(peakScores,peaksInWindow))[-1]

              bestScore = min(1.0,bestScore) / symmetry * (len(set(resonances))**2.0)                      # Put a flat bottom (top) in. 
              
          link.addPeak(closestPeak, simulatedPeak, resonances, bestScore)

  cdef void matchIntraResidual(self):
    
    cdef myDataModel DataModel
    
    cdef list intraResidualSimulatedPeakMatrix, residues, simulatedPeakList, spinsystems, listWithPresentPeaks, listWithSimulatedPeaks, listWithNotFoundSimulatedPeaks, listWithScores, contributions, resonances, peakLists, peaksInWindow
    
    cdef Residue res
    
    cdef dict allSpinSystems
    
    cdef SpinSystem spinSys
    
    cdef SimulatedPeak simulatedPeak
    
    cdef SimulatedPeakContrib contrib
    
    cdef PeakDimension dim
    
    cdef int i
    
    cdef Resonance resonance
    
    cdef Peak peak
    
    cdef int hc
    
    cdef double bestScore
    
    cdef SpinSystemLink link
    
    hc = 10000 #self.hc

    DataModel = self.DataModel
    intraResidualSimulatedPeakMatrix = self.intraResidualSimulatedPeakMatrix
    
    allSpinSystems = DataModel.spinSystems
    
    residues = DataModel.chain.residues

    for res, simulatedPeakList in zip(residues,intraResidualSimulatedPeakMatrix) :
      
      linkDict = res.intraDict
      
      for link in linkDict.values() :
        
        spinSys = link.spinSystem1

        for simulatedPeak in simulatedPeakList :
          
          contributions = simulatedPeak.simulatedPeakContribs
          
          closestPeak = None
          bestScore = 0.0
          
          resonances = [spinSys.getResonanceForAtomName(contrib.atomName) for contrib in contributions]
          

          if resonances and not None in resonances :
            
            peakLists = [ resonance.getPeaksForSpectrumDim(self,contrib.dimNumber,True) for resonance, contrib in zip(resonances , contributions) ]
           
            peaksInWindow = commonElementInLists(peakLists)                                      # And check whether 1 or more peaks that fit in one dimension actually also fit in all other dimensions. In that case the peak is in the multidimensional tolerance window

            if peaksInWindow :
              
              peakScores = [scorePeak(peak.dimensions,resonances) for peak in peaksInWindow]
                  
              bestScore, closestPeak = sorted(zip(peakScores,peaksInWindow))[-1]
              
              bestScore = min(1.0,bestScore)
              
          link.addPeak(closestPeak, simulatedPeak, resonances, bestScore)

  cdef void determineSymmetry(self) :
    
    '''
    Checking the symmetry of the spectrum. This could be done in a number of ways. But an easy way seems to be to look at the simulated peaks.
    We only care about symmetry in the sequential peaks. In almost all cases the symmetry for intra-residual correlations are the same.
    The symmetry value means the amount of peaks that represents the same correlation between atoms. In a NCOCX the value is 1 since there is only
    1 peak that represents a correlations between 3 nuclei. In for istance  a 2D C-C through space correlation the value would be 2, because a peak
    representing a correlation between 2 nuclei shows up on both sides of the diagonal. We take the maximum value we find because a diagonal peak
    might produce a lower value than the other peaks in the spectrum.
    '''

    
    cdef int symmetry, maxSymmetry
    cdef set contribSetA, contribSetB
    cdef list peaksForOneResidue
    cdef SimulatedPeak peakA, peakB
    cdef SimulatedPeakContrib contrib
    
    maxSymmetry = 1
    
    for peaksForOneResidue in self.simulatedPeakMatrix :        # We only care about symmetry in sequential peaks
      
      for peakA in peaksForOneResidue :
        
        symmetry = 0
        
        contribSetA = set([(contrib.atomName, contrib.firstOrSecondResidue) for contrib in peakA.simulatedPeakContribs])
        
        for peakB in peaksForOneResidue :
          
          contribSetB = set([(contrib.atomName, contrib.firstOrSecondResidue) for contrib in peakB.simulatedPeakContribs])
          
          if contribSetA == contribSetB :
            
            symmetry += 1
        
        if symmetry > maxSymmetry :
          
          maxSymmetry = symmetry
        
        if len(contribSetA) == len(peakA.simulatedPeakContribs) :              # Not checking every single peak when we already found a non-diagonal one. I think this is justified.
          
          self.symmetry = maxSymmetry
          
          return
        
    self.symmetry = maxSymmetry
    return
   
  def getRandomSubSetofPeaks(self, fraction):
    ''' To avoid a lot of false-positives, it is possible to leave out a
        a sub-set of the peaks in each spectrum. These sub-sets can be changed
        every run so the results will show more variability.
        
        input: float fraction
        output: list of peaks
    '''
    
    randomPeakList = self.peaks[:]
    pyrandom.shuffle(randomPeakList)
    
    return randomPeakList[0:int(fraction*len(randomPeakList)+0.5)]
        
  def __reduce__(self) :

    return (generalFactory, (type(self),), self.__getstate__())
    
  def __getstate__(self) :

    return (self.name, self.peaks, self.experimentSerial, self.serial)
  
  def __setstate__(self, state) :

    self.name, self.peaks, self.experimentSerial, self.serial = state
    
  def connectToProject(self, nmrProject) :
    
    try :
      
      self.ccpnSpectrum = nmrProject.findFirstExperiment(serial=self.experimentSerial).findFirstDataSource(serial=self.serial)
      
    except AttributeError :
      
      print "Error: Cannot find spectrum %s, navigating to peaks in this spectrum will not be possible." % self.name
      
      return

    for peak in self.peaks :
      
      peak.connectToProject()
    
    
  def getName(self):
    
    return self.name
  
  def getCcpnSpectrum(self) :
    
    return self.ccpnSpectrum
