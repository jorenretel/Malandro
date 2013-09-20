
import math
from numpy.random import randint,random_sample
import random
import math
import cython

from cython.view cimport array as cvarray

#import cProfile

#from cpython cimport bool

#from libcpp cimport bool

from ccpnmr.analysis.core.ChemicalShiftBasic import getShiftsChainProbabilities

from ccpnmr.analysis.core.MoleculeBasic import getResidueCode

from ccp.util.LabeledMolecule import getIsotopomerSingleAtomFractions, getIsotopomerAtomPairFractions,  atomPairFrac,  atomPairFractions

from ccpnmr.analysis.core.Util import getAnalysisDataDim

from generalFactory import generalFactory

cdef extern from "math.h":

  double exp(double x)
  

from libc.stdlib cimport rand, RAND_MAX

cdef double randMax = float(RAND_MAX) 






cdef class Malandro :
  
  
  cdef object chain
  
  cdef public myDataModel DataModel
  
  cdef bint useAssignments
  
  cdef bint useTentative
  
  cdef int amountOfRepeats
  
  cdef int amountOfSteps
  
  cdef object project, shiftList
  
  cdef object nmrProject
  
  cdef double minIsoFrac
  
  cdef double minTypeScore, leavePeaksOutFraction
  
  cdef list selectedSpectra
  
  cdef str sourceName
  
  cdef int hc
  
  cdef object updateInfoText
  
  cdef object addEnergyPoint
  
  cdef list acceptanceConstantList
  
  cdef double score
  
  cdef object typeSpinSystems, reTypeSpinSystems
  
  cdef object useDimenionalAssignments
  
  
  
  def __init__(self):
    
    self.hc = 10000
    
  def getResults(self):
    
    return self.getResultsC()

  cdef object getResultsC(self):
    
    return self.DataModel #.pyDataModel

  def updateSettings(self,  connector):
    
    ''''
    This method is used to fetch all important parameters from the 'connector'
    class, that stores all connects the GUI to the algorithm. All important settings
    are fetched from the GUI by the connector and passed on to the algorithm.
    This to prevent mixing the GUI and the algorithm, which might be unhandy when
    a new GUI has to be written.
    
    '''
    
    self.chain = connector.chain
    
    self.useAssignments =connector.useAssignments
    
    self.useTentative = connector.useTentative
    
    self.amountOfRepeats = connector.amountOfRepeats
    
    self.amountOfSteps = connector.amountOfSteps
    
    self.shiftList = connector.shiftList
    
    self.nmrProject = connector.nmrProject
    
    self.project = connector.project
    
    self.minIsoFrac = connector.minIsoFrac
    
    self.leavePeaksOutFraction = connector.leavePeaksOutFraction
    
    self.minTypeScore = connector.minTypeScore
    
    self.selectedSpectra = connector.selectedSpectra
    
    self.sourceName = connector.sourceName                                # Not used
    
    self.updateInfoText = connector.updateInfoText
    
    self.acceptanceConstantList = connector.acceptanceConstantList
    
    self.addEnergyPoint = connector.addEnergyPoint
    
    self.typeSpinSystems = connector.typeSpinSystems
    
    self.reTypeSpinSystems = connector.reTypeSpinSystems
    
    self.useDimenionalAssignments = connector.useDimenionalAssignments
    
  def preCalculateDataModel(self)  :
    
    self.updateInfoText('Setting up model for calculations...')
    
    self.DataModel = myDataModel(self)
        
    self.updateInfoText('Setup-up all spectra...')
    
    self.DataModel.project = self.project
    self.DataModel.nmrProject = self.nmrProject
    
    self.DataModel.setupSpectra()

    self.DataModel.setupChain()

    self.createSpinSytemsAndResonances()
        
    self.updateInfoText('Simulating spectra...')
    
    self.DataModel.setupLinks()
    
    self.simulateSpectra()
    
    self.updateInfoText('Evaluating possible dimensional contributions to peak in real spectra...')
    
    self.calculateAllPeakContributions()
    
    self.updateInfoText('Matching simulated with real spectra...')
        
    self.matchSimulatedWithRealSpectra()

    self.createJokerSpinSystems()
    
    self.updateInfoText('Scoring links between spin systems...')
    
    #self.scoreAllLinks()
    
    #self.multiplyPeakScoresWithLinkScores()
    
    self.updateInfoText('Precalculations have finished...')

  cdef void doRandomAssignment(self):
    
    self.updateInfoText('Making a random assignment...')
    
    cdef myDataModel DataModel
    
    cdef dict assignedSpinSystems
    
    cdef dict tentativeSpinSystems
    
    cdef dict justTypedSpinSystems
    
    cdef dict allSpinSystems 
    
    cdef dict jokerSpinSystems 
    
    cdef dict dictio
    
    cdef Residue res
    
    cdef int i
    
    cdef Residue resimin1
    
    cdef Residue resiplus1
    
    cdef str ccpCode
    
    cdef int seqCode
    
    cdef list listWithSpinSystems
    
    cdef SpinSystem spinSystem
    
    cdef list listWithFittingSpinSystems
    
    cdef SpinSystem randomSpinSystem
    
    cdef bint useAssignments, useTentative

    useAssignments = self.useAssignments
    useTentative = self.useTentative
    
    DataModel = self.DataModel
    
    sample = random.choice

    assignedSpinSystems = makePrivateCopyOfDictContainingLists(DataModel.previouslyAssignedSpinSystems)
    
    tentativeSpinSystems =  makePrivateCopyOfDictContainingLists(DataModel.tentativeSpinSystems)
    
    justTypedSpinSystems = makePrivateCopyOfDictContainingLists(DataModel.justTypedSpinSystems)
    
    allSpinSystems = makePrivateCopyOfDictContainingLists(DataModel.spinSystems)
    
    jokerSpinSystems = makePrivateCopyOfDictContainingLists(DataModel.jokerSpinSystems)
    
    
    if useAssignments and useTentative:
      
      dictio = mergeDictionariesContainingLists([justTypedSpinSystems,jokerSpinSystems])  
      
    elif useAssignments :
    
      dictio = mergeDictionariesContainingLists([justTypedSpinSystems, tentativeSpinSystems,jokerSpinSystems])
                                              
    elif useTentative :
      
      dictio = mergeDictionariesContainingLists([assignedSpinSystems,justTypedSpinSystems,jokerSpinSystems])
                                                
    else :
      
      dictio = makePrivateCopyOfDictContainingLists(allSpinSystems)
    
    
    for res in DataModel.chain.residues :
      
      isAssigned = False
            
      ccpCode = res.ccpCode
      seqCode = res.seqCode
      
      if  useAssignments and (ccpCode in assignedSpinSystems) :                                                         # If wanted by the user: try to put the assigned spin system in this position
        
        listWithSpinSystems = assignedSpinSystems[ccpCode]
        
        for spinSystem in listWithSpinSystems :
          
          if spinSystem.ccpnSeqCode == seqCode :
             
            spinSystem.currentResidueAssignment = res
            
            res.currentSpinSystemAssigned = spinSystem
            
            isAssigned = True
            
            listWithSpinSystems.remove(spinSystem)
            
            break
            
            
      if not isAssigned and useTentative and (ccpCode in tentativeSpinSystems) :                            # If wanted by the user: try to put a tenitative assigned spin system on this position
        
        listWithSpinSystems = tentativeSpinSystems[ccpCode]
        
        listWithFittingSpinSystems = []
        
        for spinSystem in listWithSpinSystems :
          
          if not spinSystem.currentResidueAssignment and seqCode in spinSystem.tentativeSeqCodes :
            
            listWithFittingSpinSystems.append(spinSystem)
            
        if len(listWithFittingSpinSystems) > 0 :  
            
            
          randomSpinSystem =  sample(listWithFittingSpinSystems)                                        #listWithFittingSpinSystems[random.randint(0, (len(listWithFittingSpinSystems)-1))]                   # get a random element out of the list
            
          randomSpinSystem.currentResidueAssignment = res
          
          
          res.currentSpinSystemAssigned = randomSpinSystem 
          
          isAssigned = True
          
          listWithSpinSystems.remove(randomSpinSystem)
            
          

            
      if not isAssigned :                                                                                                                                                           # If no spin system is assigned to the residue by now, assign some unassigned spinsystem (or whatever random spinsystem, when the user-made assignment should be ignored)
        
        listWithSpinSystems = []
        
        listWithFittingSpinSystems = []
        
        if ccpCode in dictio :
          
          listWithSpinSystems = dictio[ccpCode]
          
          for spinSystem in listWithSpinSystems :
            
            if not spinSystem.currentResidueAssignment :
              
              listWithFittingSpinSystems.append(spinSystem)    
            
      
        if len(listWithFittingSpinSystems) > 0 :      
          
                      
          randomSpinSystem = sample( listWithFittingSpinSystems)                                                        #listWithFittingSpinSystems[random.randint(0, (len(listWithFittingSpinSystems)-1))] 
             
          randomSpinSystem.currentResidueAssignment = res
          
          res.currentSpinSystemAssigned = randomSpinSystem 
          
          isAssigned = True
          
          if randomSpinSystem in listWithFittingSpinSystems :
          
            listWithSpinSystems.remove(randomSpinSystem)
            
            
          
      if not isAssigned :
          
        print 'something went wrong during random assignment at the start of the procedure, not all residues have an assignment'
        print res.seqCode
        print res.ccpCode

  def runAnnealling(self):
    
    useAssignments = self.useAssignments
    useTentative = self.useTentative    
    amountOfStepsPerTemperature = self.amountOfSteps
    AcceptanceConstantList = self.acceptanceConstantList
    DataModel = self.DataModel
            
    allSpinSystems = DataModel.spinSystems

    allSpinSystemsWithoutAssigned = DataModel.allSpinSystemsWithoutAssigned

        
    if useAssignments :                                                 # If the already made assignments are used, the dictionary will just include unassigned spinsystems and Joker spinsystems
      relevantSpinSystems = allSpinSystemsWithoutAssigned
    else :                                                                       # If the already assigned are chosen to still be permutated, the dictionary will also include the ready assinged spinsystems
      relevantSpinSystems = allSpinSystems
      
    listWithSpinSystems = list(set([val for subl in relevantSpinSystems.values() for val in subl]))
    
    for x,AcceptanceConstant in enumerate(AcceptanceConstantList) :
      
      self.annealingSub(AcceptanceConstant,amountOfStepsPerTemperature,listWithSpinSystems)
      
      self.scoreInitialAssignment()
      
      self.addEnergyPoint(self.score,x)
 
  def startMonteCarlo(self):
    
    cdef int repeat
    
    cdef Residue res
    
    cdef str info
    
    repeat = self.amountOfRepeats
    
    info = 'Running annealing number %s out of ' + str(repeat) + '...'
    
    self.setupSpinSystemExchange()

    for x in range(repeat) :
      
      self.cleanAssignments()
      
      self.leaveOutSubSetOfPeaks(self.leavePeaksOutFraction)
      
      self.scoreAllLinks()

      self.doRandomAssignment()

      self.setupPeakInformationForRandomAssignment()

      self.scoreInitialAssignment()

      self.updateInfoText(info %str(x+1))

      self.runAnnealling()
      
      i = 1
      matches = 0

      for res in self.DataModel.chain.residues :
        
        res.solutions.append(res.currentSpinSystemAssigned)
        
        res.currentSpinSystemAssigned.solutions.append(res.seqCode)
        
        
        if res.ccpnResidue.findFirstResonanceGroup() :

          if res.currentSpinSystemAssigned.spinSystemNumber == res.ccpnResidue.findFirstResonanceGroup().serial :

            matches = matches + 1
        i = i + 1
    
    self.updateInfoText('Done')

  cdef void cleanAssignments(self) :
    
    cdef myDataModel DataModel
    
    cdef Chain chain
    
    cdef list residues
    
    cdef list spectra
    
    cdef Spectrum spectrum
    
    cdef dict mySpinSystems
    
    cdef dict tentativeSpinSystems
    
    cdef str key
    
    cdef list spinSystemList
    
    cdef SpinSystem spinSys
    
    cdef Peak peak
    
    cdef Residue res
    
    
    
    DataModel = self.DataModel
    
    chain = DataModel.chain
    
    residues = chain.residues
    
    spectra = DataModel.spectra
    
    mySpinSystems = DataModel.spinSystems
    
    tentativeSpinSystems = DataModel.tentativeSpinSystems
    
    for key,  spinSystemList in mySpinSystems.items() :                                 # de-assigning residues from spinsystems
      
      for spinSys in spinSystemList :

        spinSys.currentResidueAssignment = None
        
        
    for key,  spinSystemList in tentativeSpinSystems.items() :                                 # also de-assigning residues from tentative assigned spinsystems
      
      for spinSys in spinSystemList :

        spinSys.currentResidueAssignment = None
          
        
    for res in residues :                                                                       # de-assigning spinsystems from residues
      
      res.currentSpinSystemAssigned = None
      
      
    for spectrum in spectra :                                                               # Setting all peak-degeneracies back to 0
      
      for peak in spectrum.peaks :
        
        peak.degeneracy = 0

  cdef void setupPeakInformationForRandomAssignment(self):
    
    cdef list residues
    
    cdef Residue res
      
    cdef SpinSystemLink link
    
    cdef Peak peak
  
    cdef PeakLink pl
    
    residues = self.DataModel.chain.residues
      
    for res in residues :
      
      nextRes = res.nextResidue
      
      link = res.getFromLinkDict(res.currentSpinSystemAssigned, nextRes.currentSpinSystemAssigned)
 
      for pl in link.activePeakLinks :
        
        pl.peak.degeneracy += 1
      
        #for peak in link.realPeaks :
        
        #  peak.degeneracy += 1

  cdef void createJokerSpinSystems(self):
    
    '''
    When there are less spin systems that are typed to a
    certain amino acid type than there are residues of that
    type in the sequence, some Jokers have to be introduced. 
    
    '''
    
    cdef myDataModel DataModel
    
    cdef str ccpCode
    
    cdef int amountOfAssignedSpinsystems, NTypedSpinSystems, NResiduesOfThisType, short, i, x
        
    cdef SpinSystem spinSys, newSpinSystem
    
    DataModel = self.DataModel    
    
    
    i = 1
    
    for ccpCode, NResiduesOfThisType in DataModel.chain.residueTypeFrequencyDict.items() : #DataModel.spinSystems.keys() :

      NResiduesAssigned = len(set([spinSys.ccpnSeqCode for spinSys in DataModel.previouslyAssignedSpinSystems.get(ccpCode, [])]))
      NTypedSpinSystems = len(DataModel.justTypedSpinSystems.get(ccpCode,[]))

      
      #short = amountOfResiduesOfThisType- (amountOfAssignedSpinsystems + amountOfTypedSpinSystems + amountOfTentativeSpinSystemsWithOnlyOneCcpCode)
      short = NResiduesOfThisType - (NResiduesAssigned + NTypedSpinSystems)
      
      
      if short < 0 :
        
        string = 'You seem to have more ' + ccpCode + ' spin systems than residues of this type are in the sequence'

        short = 0
      
      for x in range(short) :
        
        newSpinSystem = SpinSystem(DataModel=DataModel,ccpCode=ccpCode)
        
        newSpinSystem.spinSystemNumber = i * 1000000
        
        i = i + 1
        
        DataModel.addToDictWithLists(DataModel.spinSystems,ccpCode,newSpinSystem)
        DataModel.addToDictWithLists(DataModel.jokerSpinSystems,ccpCode,newSpinSystem)
        DataModel.addToDictWithLists(DataModel.allSpinSystemsWithoutAssigned,ccpCode,newSpinSystem)

        

        #if ccpCode in DataModel.spinSystems :
        #
        #  DataModel.spinSystems[ccpCode].append(newSpinSystem)
        #  
        #else :
        #  
        #  DataModel.spinSystems[ccpCode] = [newSpinSystem]
        #  
        #  
        #if ccpCode in DataModel.jokerSpinSystems:
        #  
        #  DataModel.jokerSpinSystems[ccpCode].append(newSpinSystem)
        #  
        #else :
        #  
        #  DataModel.jokerSpinSystems[ccpCode] = [newSpinSystem]
        #  
        #  
        #if ccpCode in DataModel.allSpinSystemsWithoutAssigned:
        #  
        #  DataModel.allSpinSystemsWithoutAssigned[ccpCode].append(newSpinSystem)
        #  
        #else :
        #  
        #  DataModel.allSpinSystemsWithoutAssigned[ccpCode] = [newSpinSystem]

    string = str(i-1) + ' joker spinsystems are used in this calculation.'    

  cdef void simulateSpectra(self) :
    
    cdef myDataModel DataModel
    
    cdef Spectrum spectrum
    
    DataModel = self.DataModel
    
    for spectrum in DataModel.spectra:                 
      
      self.updateInfoText('Simulating ' + spectrum.name)
      
      spectrum.simulate()
      
      spectrum.determineSymmetry()

  cdef void createSpinSytemsAndResonances(self):
    
    self.DataModel.setupSpinSystems(minTypeScore=self.minTypeScore)

  cdef void scoreAllLinks(self):
    
    cdef list residues 
    cdef Residue res
    cdef dict linkDict
    cdef SpinSystemLink linkObject

    DataModel = self.DataModel

    residues = DataModel.chain.residues
    
    for res in residues :
      
      linkDict = res.linkDict
      
      for linkObject in linkDict.values() :
        
        linkObject.determineScore()
             
  def leaveOutSubSetOfPeaks(self, fraction):
    
    cdef Spectrum spectrum
    cdef Residue residue
    cdef SpinSystemLink spinSystemLink
    cdef PeakLink peakLink
    cdef set peaksToLeaveOut
    cdef list peaksToLeaveOutList
    
    DataModel = self.DataModel
    residues = DataModel.chain.residues
    spectra = DataModel.spectra
    
    if fraction == 0.0 :
      
      for residue in residues :
      
        for spinSystemLink in residue.linkDict.values():
          
          spinSystemLink.activePeakLinks = spinSystemLink.peakLinks
      
    else :
    
      peaksToLeaveOutList = []
      
      for spectrum in spectra:
        
        peaksToLeaveOutList.extend(spectrum.getRandomSubSetofPeaks(fraction))
        
        
      peaksToLeaveOut = set(peaksToLeaveOutList)
      
      for residue in residues :
        
        for spinSystemLink in residue.linkDict.values():
          
          activePeakLinks = []
          
          for peakLink in spinSystemLink.peakLinks:
            
            if not peakLink.peak in peaksToLeaveOut :
              
              activePeakLinks.append(peakLink)
          
          spinSystemLink.activePeakLinks = activePeakLinks
      
  cdef void matchSimulatedWithRealSpectra(self):
    
    self.updateInfoText('Matching real to simulated spectra.....')
    
    cdef Spectrum spectrum
    
    DataModel = self.DataModel
    
    for spectrum in DataModel.spectra :
      
      self.updateInfoText('Match simulated with real spectrum: ' + spectrum.name)
      
      spectrum.match()

  cdef void calculateAllPeakContributions(self):
        
    cdef Spectrum spectrum
    
    cdef str info

    for spectrum in self.DataModel.spectra:                 # Determine for each dimension of every peak in all (used) spectra, which resonances can contribute to the peak
              
      info = 'Evaluating dimensional contributions to peaks: ' + spectrum.name
      
      self.updateInfoText(info)
      
      spectrum.findPossiblePeakContributions(self.useDimenionalAssignments)

  cdef void setupSpinSystemExchange(self) :
    
    cdef SpinSystem spinSystem
    
    spinSystems = self.DataModel.spinSystems
    
    allSpinSystems = []
    
    for spinSystemList in spinSystems.values() :
      
      allSpinSystems += spinSystemList
      
    allSpinSystems = list(set(allSpinSystems))
    
    for spinSystem in allSpinSystems :
        
      spinSystem.setupAllowedResidues(self.useAssignments, self.useTentative)
      spinSystem.setupAllowedResidueView()

    for spinSystem in allSpinSystems :
    
      spinSystem.setupExchangeSpinSystems(self.useAssignments)

  cdef double getAtomLabellingFraction(self,str molType, str ccpCode, str atomName, object labellingScheme):
    """
    Get the fraction of labelling for a given atom in a given amino acid.

    .. describe:: Input
  
    Nmr.Resonance
  
    .. describe:: Output
  
    Float
    """


    if labellingScheme == None :
      fraction = 1.0 # When no specific labelling scheme is chosen.
  

    else :

      chemCompLabel = labellingScheme.findFirstChemCompLabel(ccpCode=ccpCode,
                                                             molType=molType)
 
      if not chemCompLabel:
        natAbun = self.NatAbun #resonance.root.findFirstLabelingScheme(name='NatAbun')
 
        if natAbun:
          chemCompLabel = natAbun.findFirstChemCompLabel(ccpCode=ccpCode,
                                                         molType=molType)
 
      if chemCompLabel:
        isotopomers = chemCompLabel.isotopomers
        if atomName[0] == 'H':
          isotope = '1H'
        elif atomName[0] == 'C':
          isotope = '13C'
        elif atomName[0] == 'N':
          isotope = '15N'
        elif atomName[0] == 'P':
          isotope = '31P'
 
        fracDict = getIsotopomerSingleAtomFractions(isotopomers,
                                                      atomName, 1)

        fraction = fracDict.get(isotope, 1.0)
 

          
    return fraction

  cdef double getAtomPairLabellingFraction_intra(self,str molType, str ccpCode, str atomNameA, str atomNameB, object labellingScheme):
    """
    Get the fraction of a pair of atoms both being labelled
    given a labelling scheme. Considers individual isotopomers if
    the resonances are within the same residue.

    .. describe:: Input
  
 
  
    .. describe:: Output
  
    Float
    """
  
    if labellingScheme == None :
      fraction = 1.0 # When no specific labelling scheme is chosen.  




    else:
      findFirstChemCompLabel = labellingScheme.findFirstChemCompLabel
  
      chemCompLabel = findFirstChemCompLabel(ccpCode=ccpCode,
                                             molType=molType)
 
      if not chemCompLabel:
        natAbun = self.NatAbun #resonanceA.root.findFirstLabelingScheme(name='NatAbun')
 
        if natAbun:
          chemCompLabel = natAbun.findFirstChemCompLabel(ccpCode=ccpCode,
                                                         molType=molType)

 
      if chemCompLabel:
        isotopomers  = chemCompLabel.isotopomers
        if atomNameA[0] == 'H':
          isotope = '1H'
        elif atomNameA[0] == 'C':
          isotope = '13C'
        elif atomNameA[0] == 'N':
          isotope = '15N'
        elif atomNameA[0] == 'P':
          isotope = '31P'

        if atomNameB[0] == 'H':
          isotope_2 = '1H'
        elif atomNameB[0] == 'C':
          isotope_2 = '13C'
        elif atomNameB[0] == 'N':
          isotope_2 = '15N'
        elif atomNameB[0] == 'P':
          isotope_2 = '31P'

        atomNames = (atomNameA, atomNameB)
        isotopes = (isotope, isotope_2)

        pairDict  = getIsotopomerAtomPairFractions(isotopomers, atomNames, (1,1))
        fraction = pairDict.get(isotopes, 1.0)

      
    return fraction

  cdef double getAtomPairLabellingFraction_inter(self,str molType, str ccpCodeA, str ccpCodeB, str atomNameA, str atomNameB, object labellingScheme):
    fractionA = self.getAtomLabellingFraction(molType, ccpCodeA, atomNameA, labellingScheme)
    fractionB = self.getAtomLabellingFraction(molType, ccpCodeB, atomNameB, labellingScheme)
    fraction = fractionA * fractionB
      
    return fraction

  cdef double getIsotopomerSingleAtomFractionsForAtom(self, set isotopomers, Atom atom, str isotopeCode) :                # Not used at the moment
    
    cdef double isoWeightSum

    cdef object isotopomer
    
    cdef double labellingFraction
    
    labellingFraction = 0.0
    
    isoWeightSum =  float(sum([isotopomer.weight for isotopomer in isotopomers]))

    for isotopomer in isotopomers :
      
      labellingFraction += atom.labelInfoTemp[isotopomer][isotopeCode] * isotopomer.weight / isoWeightSum
    
    return labellingFraction 

  @cython.wraparound(False)
  @cython.boundscheck(False)
  @cython.cdivision(True)
  @cython.nonecheck(False)
  cdef void annealingSub(self, double AcceptanceConstant,int amountOfStepsPerTemperature,list listWithSpinSystems):
     
    cdef int improvements, equals, worse, r, seqCodeA, seqCodeB, deltaLinkScore, lengthOfListWithSpinSystems
    
    cdef double score, DeltaScore
    
    cdef list exchangeSpinSystems, oldPeaks, newPeaks, peakSet
    
    cdef SpinSystem A, B, Am1, Ap1, Bm1, Bp1
  
    cdef Residue currentResidueA, currentResidueB, previousResA, previousResB, nextResA, nextResB
        
    cdef SpinSystemLink l1, l2, l3, l4, l5, l6, l7, l8
    
    cdef Peak peak
    
    cdef PeakLink pl
    
    cdef int path

    lengthOfListWithSpinSystems = len(listWithSpinSystems)
    score = self.score
    improvements = 0
    equals = 0
    worse = 0
    
    for aTry in xrange(amountOfStepsPerTemperature) :
      
      r = int(rand()/(randMax+1)*lengthOfListWithSpinSystems)
      A = <SpinSystem>listWithSpinSystems[r]

      exchangeSpinSystems = A.exchangeSpinSystems
      
      if exchangeSpinSystems :
        
        r = int(rand()/(randMax+1)*len(exchangeSpinSystems)) #int(rand()/(RAND_MAX*len(exchangeSpinSystems)))
        
        B = <SpinSystem>exchangeSpinSystems[r]
        
        
      else :
      
        continue
      
      currentResidueA = A.currentResidueAssignment
      currentResidueB = B.currentResidueAssignment
      
      if not currentResidueA is None and not currentResidueB is None :             # Both spin systems are assigned to a residue
        
        path = 0
        
        seqCodeA = currentResidueA.seqCode
        seqCodeB = currentResidueB.seqCode
        
        if not ( A.allowedResidueView[seqCodeB] and B.allowedResidueView[seqCodeA] ) :
                    
          continue
        
        previousResA = currentResidueA.previousResidue
        previousResB = currentResidueB.previousResidue
        
        nextResA = currentResidueA.nextResidue
        nextResB = currentResidueB.nextResidue
        
        oldLinkScore = 0
        newLinkScore = 0
       
        if currentResidueA is previousResB :                                                                            # A and B happen to be a sequential pair AB
  
          Am1 = previousResA.currentSpinSystemAssigned
          Bp1 = nextResB.currentSpinSystemAssigned
          
          l1 = currentResidueA.getFromLinkDict(A,B)
          l2 = previousResA.getFromLinkDict(Am1,A)
          l3 = currentResidueB.getFromLinkDict(B,Bp1)
          l4 = currentResidueA.getFromLinkDict(B,A)
          l5 = previousResA.getFromLinkDict(Am1,B)
          l6 = currentResidueB.getFromLinkDict(A,Bp1)
          
          #deltaLinkScore = l4.score + l5.score + l6.score - (l1.score + l2.score + l3.score)
          oldPeaks = l1.activePeakLinks + l2.activePeakLinks + l3.activePeakLinks
          newPeaks = l4.activePeakLinks + l5.activePeakLinks + l6.activePeakLinks
          peakSet = oldPeaks + newPeaks
          
        elif currentResidueB is previousResA :                                                                                    # sequential pair BA
  
          Bm1 = previousResB.currentSpinSystemAssigned
          Ap1 = nextResA.currentSpinSystemAssigned
          
          l1 = currentResidueB.getFromLinkDict(B,A)
          l2 = previousResB.getFromLinkDict(Bm1,B)
          l3 = currentResidueA.getFromLinkDict(A,Ap1)
          l4 = currentResidueB.getFromLinkDict(A,B)
          l5 = previousResB.getFromLinkDict(Bm1,A)
          l6 = currentResidueA.getFromLinkDict(B,Ap1)
          
          #deltaLinkScore = l4.score + l5.score + l6.score - (l1.score + l2.score + l3.score)
          oldPeaks = l1.activePeakLinks + l2.activePeakLinks + l3.activePeakLinks
          newPeaks = l4.activePeakLinks + l5.activePeakLinks + l6.activePeakLinks
          peakSet = oldPeaks + newPeaks

        else :                                                                                                            # A and B are not sequential
        
          Am1 = previousResA.currentSpinSystemAssigned
          Ap1 = nextResA.currentSpinSystemAssigned
          Bm1 = previousResB.currentSpinSystemAssigned
          Bp1 = nextResB.currentSpinSystemAssigned

          l1 = previousResA.getFromLinkDict(Am1,A)
          l2 = currentResidueA.getFromLinkDict(A,Ap1)
          l3 = previousResB.getFromLinkDict(Bm1,B)
          l4 = currentResidueB.getFromLinkDict(B,Bp1)
          l5 = previousResA.getFromLinkDict(Am1,B)
          l6 = currentResidueA.getFromLinkDict(B,Ap1)
          l7 = previousResB.getFromLinkDict(Bm1,A)
          l8 = currentResidueB.getFromLinkDict(A,Bp1)
          
          #deltaLinkScore = l5.score + l6.score + l7.score + l8.score - (l1.score + l2.score + l3.score + l4.score)
          oldPeaks = l1.activePeakLinks + l2.activePeakLinks + l3.activePeakLinks + l4.activePeakLinks
          newPeaks = l5.activePeakLinks + l6.activePeakLinks + l7.activePeakLinks + l8.activePeakLinks
          peakSet = oldPeaks + newPeaks
          
      elif not currentResidueA is None :                                                                      # spin system B is not assigned to any residue
        
        path = 1
        
        seqCodeA = currentResidueA.seqCode
        
        if not B.allowedResidueView[seqCodeA] :
          
          continue
          
        previousResA = currentResidueA.previousResidue
        nextResA = currentResidueA.nextResidue
  
        Am1 = previousResA.currentSpinSystemAssigned
        Ap1 = nextResA.currentSpinSystemAssigned
        
        l1  = previousResA.getFromLinkDict(Am1,A)
        l2  = currentResidueA.getFromLinkDict(A,Ap1)
        l3  = previousResA.getFromLinkDict(Am1,B)
        l4  = currentResidueA.getFromLinkDict(B,Ap1)
        
        #deltaLinkScore = l3.score + l4.score - (l1.score + l2.score)
        oldPeaks = l1.activePeakLinks + l2.activePeakLinks
        newPeaks = l3.activePeakLinks + l4.activePeakLinks
        peakSet = oldPeaks+newPeaks
  
      elif not currentResidueB is None :                                                                      # spin system A is not assigned to any residue
        
        path = 2
        
        seqCodeB = currentResidueB.seqCode
        
        if not A.allowedResidueView[seqCodeB] :          
          continue
          
        previousResB = currentResidueB.previousResidue
        nextResB = currentResidueB.nextResidue
  
        Bm1 = previousResB.currentSpinSystemAssigned
        Bp1 = nextResB.currentSpinSystemAssigned
        
        l1  = previousResB.getFromLinkDict(Bm1,B)
        l2  = currentResidueB.getFromLinkDict(B,Bp1)
        l3  = previousResB.getFromLinkDict(Bm1,A)
        l4  = currentResidueB.getFromLinkDict(A,Bp1)

        #deltaLinkScore = l3.score + l4.score - (l1.score + l2.score)
        oldPeaks = l1.activePeakLinks + l2.activePeakLinks
        newPeaks = l3.activePeakLinks + l4.activePeakLinks
        peakSet = oldPeaks+newPeaks

      else :
        
        continue
      
      DeltaScore = CcalcDeltaPeakScore(peakSet,oldPeaks,newPeaks) #+ deltaLinkScore
      
      if DeltaScore >= 0 or exp(AcceptanceConstant*DeltaScore) > rand()/randMax :
        
        score += DeltaScore
        
        for pl in peakSet :
          
          pl.peak.degeneracy = pl.peak.degeneracyTemp
          
        B.currentResidueAssignment = currentResidueA
        A.currentResidueAssignment = currentResidueB
          
        if path == 0 :

          currentResidueA.currentSpinSystemAssigned = B
          currentResidueB.currentSpinSystemAssigned = A
          
        elif path == 1 :

          currentResidueA.currentSpinSystemAssigned = B
          
        else :

          currentResidueB.currentSpinSystemAssigned = A
          
    self.score = score

  cdef void scoreInitialAssignment(self) :
    
    cdef list residues
    
    cdef Residue res
    
    cdef Residue nextRes
    
    cdef double score
    
    cdef dict linkDict
    
    cdef int keyA
    
    cdef int keyB
    
    cdef int key
    
    cdef SpinSystemLink link
    
    cdef list peaks
    
    cdef Peak peak
    
    cdef double peakScore
    
    cdef PeakLink pl
    
    residues = self.DataModel.chain.residues
    
    score = 0.0
    
    for res in residues :
      
      nextRes = res.nextResidue
        
      link = res.getFromLinkDict(res.currentSpinSystemAssigned, nextRes.currentSpinSystemAssigned)
      
      score += sum([1.0/pl.peak.degeneracy * pl.preMultipliedScore for pl in link.activePeakLinks]) #+ link.score  
          
    self.score = score      
          
  cdef multiplyPeakScoresWithLinkScores(self) : # Not used, already happens during spinSystemLink.determineScore()
    
    cdef list residues, spinSystemLinks, peakLinks
    cdef Residue res
    cdef dict linkDict
    cdef SpinSystemLink spinSystemLink
    cdef int spinSystemLinkScore
    cdef PeakLink peakLink
    
    residues = self.DataModel.chain.residues
    
    for res in residues :
      
      linkDict = res.linkDict
      
      spinSystemLinks = linkDict.values()
      
      for spinSystemLink in spinSystemLinks :
        
        spinSystemLinkScore = spinSystemLink.score
        
        peakLinks = spinSystemLink.peakLinks
        
        for peakLink in peakLinks:
          
          peakLink.score *= spinSystemLinkScore
 