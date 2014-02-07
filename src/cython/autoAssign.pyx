
import math
#from numpy.random import randint,random_sample
import random as pyrandom
import math
import cython
from cython.view cimport array as cvarray
from time import time
from ccpnmr.analysis.core.ChemicalShiftBasic import getShiftsChainProbabilities
from ccpnmr.analysis.core.MoleculeBasic import getResidueCode
from ccp.util.LabeledMolecule import getIsotopomerSingleAtomFractions, getIsotopomerAtomPairFractions,  atomPairFrac,  atomPairFractions
from ccpnmr.analysis.core.Util import getAnalysisDataDim
from generalFactory import generalFactory

cdef extern from "math.h":

  double exp(double x)
  
from libc.stdlib cimport rand, srand, RAND_MAX

cdef double randMax = float(RAND_MAX)

acceptanceConstants = [0.0, 0.01, 0.015, 0.022, 0.033,
                       0.050, 0.075, 0.113, 0.170, 0.256,
                       0.384, 0.576, 0.864, 1.297, 1.946,
                       2.919, 4.378, 6.568, 9.852, 14.77,
                       22.16, 33.25]

cdef class Malandro :
  
  cdef public myDataModel DataModel
  cdef bint useAssignments, useTentative, typeSpinSystems, reTypeSpinSystems, useDimenionalAssignments
  cdef object project, shiftList, nmrProject, chain
  cdef double minTypeScore, minIsoFrac
  cdef list selectedSpectra
  cdef int hc
  cdef double score
  cdef object textObservers, energyObservers
  
  def __init__(self):
    
    self.hc = 10000
    self.energyObservers = []
    self.textObservers = []
    
  def getResults(self):
    
    return self.getResultsC()

  cdef object getResultsC(self):
    
    return self.DataModel

  def updateSettings(self,  connector):
    ''''This method is used to fetch all important parameters from the 'connector'
    class, that stores all connects the GUI to the algorithm. All important settings
    are fetched from the GUI by the connector and passed on to the algorithm.
    This to prevent mixing the GUI and the algorithm, which might be unhandy when
    a new GUI has to be written.
    
    '''
    
    self.chain = connector.chain
    self.useAssignments =connector.useAssignments
    self.useTentative = connector.useTentative
    self.shiftList = connector.shiftList
    self.nmrProject = connector.nmrProject
    self.project = connector.project          #Not really used, nmrProject is sufficient
    self.minIsoFrac = connector.minIsoFrac
    self.minTypeScore = connector.minTypeScore
    self.selectedSpectra = connector.selectedSpectra
    self.typeSpinSystems = connector.typeSpinSystems
    self.reTypeSpinSystems = connector.reTypeSpinSystems
    self.useDimenionalAssignments = connector.useDimenionalAssignments
    
  def preCalculateDataModel(self):
    '''Creates the data model for the calculations and sequentially
       sets up everything that is needed to run the annealing later,
       like simulating the spectra, finding out which resonances could
       contribute to which peak dimensions. Matching the simulated
       with the real spectra and introducing joker spin systems.
       
    '''
    sendText = self.notifyTextObservers
    
    sendText('Setting up model for calculations...')
    self.DataModel = myDataModel()
    sendText('Setup-up all spectra...')
    self.DataModel.project = self.project
    self.DataModel.nmrProject = self.nmrProject
    self.DataModel.setupSpectra(self.selectedSpectra)
    self.DataModel.setupChain(self.chain)
    self.DataModel.setupSpinSystems(resonanceGroups=self.nmrProject.resonanceGroups, shiftList=self.shiftList,
                                    useAssignments=self.useAssignments, useTentative=self.useTentative,
                                    useType=not self.reTypeSpinSystems, includeUntypedSpinSystems=self.typeSpinSystems,
                                    minTypeScore=self.minTypeScore, makeJokers=True)
    sendText('Simulating spectra...')
    self.DataModel.setupLinks()
    self.simulateSpectra(minIsoFrac=self.minIsoFrac)
    sendText('Evaluating possible dimensional contributions to peak in real spectra...')
    self.calculateAllPeakContributions()
    sendText('Matching simulated with real spectra...')
    self.matchSimulatedWithRealSpectra()
    self.setupSpinSystemExchange()
    #self.createJokerSpinSystems()
    #sendText('Scoring links between spin systems...')
    #self.scoreAllLinks()
    #self.multiplyPeakScoresWithLinkScores()
    sendText('Precalculations have finished...')

  def doRandomAssignment(self):
    '''Performs a random assignment that is consistent in
       terms of amino acid types.
       
    '''
    self.notifyTextObservers('Making a random assignment...')
    
    cdef myDataModel DataModel
    cdef Residue residue
    cdef SpinSystem spinSystem
    
    choice = pyrandom.choice
    shuffle = pyrandom.shuffle
    DataModel = self.DataModel
    
    usedSpinSystems = set()
    residues = DataModel.chain.residues[:]
    shuffle(residues)
    
    for residue in residues:
      
      unUsedPossibilities = residue.allowedSpinSystems - usedSpinSystems
      if not unUsedPossibilities:
        print 'Something went wrong during the random assignment with %s%s' %(residue.seqCode, residue.ccpCode)
        print '%s%s' %(2,'swd')
        continue
      spinSystem = choice(list(unUsedPossibilities))
      spinSystem.currentResidueAssignment = residue
      residue.currentSpinSystemAssigned = spinSystem
      usedSpinSystems.add(spinSystem)

  def runAnnealling(self, stepsPerTemperature=10000, acceptanceConstants=acceptanceConstants):
    '''Runs one full annealling by going through a list
       of temperature constants and calling annealingSub
       with it.
       
    '''

    cdef SpinSystem spinSystem
    
    # Only spin systems that can change assignment during the monte carlo procedure
    # are relevant, all others will stay were they are anyway.
    relevantSpinSystems = [spinSystem for spinSystem in self.DataModel.getSpinSystemSet() if spinSystem.exchangeSpinSystems]

    # For every annealing I create a new instance of the mersenne twister random number generator.
    # As a seed I use a number from the linear congruential generator present in the c standard
    # library.
    rng = Random()
    
    if relevantSpinSystems:
      # doing one annealing
      for x,acceptanceConstant in enumerate(acceptanceConstants) :
        
        rng.seed(rand())
        self.annealingSub(acceptanceConstant,stepsPerTemperature,relevantSpinSystems,rng)
        self.scoreInitialAssignment()
        self.notifyEnergyObservers(self.score,x+1)
    
    else:
      self.notifyTextObservers('Nothing to do.')      #todo, more explanation to user here.
    
    # storing the result
  
  def storeResults(self):
    
    cdef Residue residue
    cdef SpinSystem spinSystem
    
    for residue in self.DataModel.chain.residues :
      residue.solutions.append(residue.currentSpinSystemAssigned)
      
    for spinSystem in self.DataModel.getSpinSystemSet():
      spinSystem.solutions.append(spinSystem.currentResidueAssignment.seqCode)
    
    self.DataModel.energies.append(self.score*-1)
    
      
  def startMonteCarlo(self, amountOfRuns=1, stepsPerTemperature=10000, acceptanceConstants=acceptanceConstants, fractionOfPeaksLeftOut=0.0):
    '''Run the optimization for a number of times
       as defined in amountOfRuns.
    
    '''
    
    info = 'Running annealing number %s out of ' + str(amountOfRuns) + '...'    
    srand(int(time()*1000000%10000000))        # Seeding the linear congruential pseudo-random number generator

    for run in range(amountOfRuns) :
      
      self.cleanAssignments()
      self.leaveOutSubSetOfPeaks(fractionOfPeaksLeftOut)
      self.scoreAllLinks()
      self.doRandomAssignment()
      self.setupPeakInformationForRandomAssignment()
      self.notifyTextObservers(info %str(run+1))
      self.scoreInitialAssignment()
      self.notifyEnergyObservers(self.score,0)
      self.runAnnealling(stepsPerTemperature=stepsPerTemperature,acceptanceConstants=acceptanceConstants)
      self.storeResults()
 
    self.notifyTextObservers('Done')

  cdef void cleanAssignments(self):
    '''Cleans all assignment information generated
       by previous runs of the annealing.
       
    '''
    
    cdef Spectrum spectrum
    cdef SpinSystem spinSystem
    cdef Peak peak
    cdef Residue residue

    # de-assigning residues from spin systems
    for spinSystem in self.DataModel.getSpinSystemSet():
      spinSystem.currentResidueAssignment = None
      
    # de-assigning spin systems from residues
    for residue in self.DataModel.chain.residues :
      residue.currentSpinSystemAssigned = None
      
    # Setting all peak-degeneracies back to 0
    for spectrum in self.DataModel.spectra :                                                               
      for peak in spectrum.peaks :
        peak.degeneracy = 0

  cdef void setupPeakInformationForRandomAssignment(self):
    '''Sets the degeneracy property of each peak to
       the correct value. This value will later on be
       updated whenever a change is made during the
       annealing.
       
    '''
    
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

  def simulateSpectra(self, minIsoFrac=0.0) :
    '''Simulates all spectra, in the sense that
       a list with simulated peaks is created.
    
    '''
    
    cdef myDataModel DataModel
    cdef Spectrum spectrum
    DataModel = self.DataModel
    
    for spectrum in DataModel.spectra:
      self.notifyTextObservers('Simulating ' + spectrum.name)
      spectrum.simulate(minIsoFrac=minIsoFrac)
      spectrum.determineSymmetry()

  cdef void scoreAllLinks(self):
    '''Scores every link in the model.
    '''
    
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
    '''Can be ran before the anneling to leave out
       a random subset of the peaks. This can be used to
       generate more heterogeneous results over multiple
       runs.
       
    '''
    
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
    '''Combine lists of simulated peaks with resonances
       of each two spin
    '''
    
    cdef Spectrum spectrum
    
    DataModel = self.DataModel
    
    for spectrum in DataModel.spectra :
      
      self.notifyTextObservers('Match simulated with real spectrum: ' + spectrum.name)
      
      spectrum.match()

  cdef void calculateAllPeakContributions(self):
        
    cdef Spectrum spectrum
    cdef str info

    for spectrum in self.DataModel.spectra:                 # Determine for each dimension of every peak in all (used) spectra, which resonances can contribute to the peak
              
      info = 'Evaluating dimensional contributions to peaks: ' + spectrum.name
      self.notifyTextObservers(info)
      
      spectrum.findPossiblePeakContributions(self.useDimenionalAssignments)

  cdef void setupSpinSystemExchange(self):
    
    cdef SpinSystem spinSystem
    spinSystems = self.DataModel.getSpinSystemSet()
    for spinSystem in spinSystems :
      spinSystem.setupExchangeSpinSystems()

  @cython.wraparound(False)
  @cython.boundscheck(False)
  @cython.cdivision(True)
  @cython.nonecheck(False)
  cdef void annealingSub(self, double AcceptanceConstant,int stepsPerTemperature,list listWithSpinSystems, Random rng):
     
    cdef int improvements, equals, worse, r, seqCodeA, seqCodeB, deltaLinkScore, lengthOfListWithSpinSystems, path
    cdef double score, DeltaScore
    cdef list exchangeSpinSystems, oldPeaks, newPeaks, peakSet
    cdef SpinSystem A, B, Am1, Ap1, Bm1, Bp1
    cdef Residue currentResidueA, currentResidueB, previousResA, previousResB, nextResA, nextResB
    cdef SpinSystemLink l1, l2, l3, l4, l5, l6, l7, l8
    cdef Peak peak
    cdef PeakLink pl

    lengthOfListWithSpinSystems = len(listWithSpinSystems)
    score = self.score
    improvements = 0
    equals = 0
    worse = 0
    
    for aTry in xrange(stepsPerTemperature) :
      
      r = rng.cy_randrange(0,lengthOfListWithSpinSystems,1)
      #r = int(rand()/(randMax+1)*lengthOfListWithSpinSystems)
      A = <SpinSystem>listWithSpinSystems[r]
      exchangeSpinSystems = A.exchangeSpinSystems
      r = rng.cy_randrange(0,len(exchangeSpinSystems),1)
      B = <SpinSystem>exchangeSpinSystems[r]

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
      
      if DeltaScore >= 0 or exp(AcceptanceConstant*DeltaScore) > rng.random() :
        
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

  cdef void scoreInitialAssignment(self):
    
    cdef list residues
    cdef Residue res
    cdef Residue nextRes
    cdef double score
    cdef SpinSystemLink link
    cdef PeakLink pl
    
    residues = self.DataModel.chain.residues
    score = 0.0
    
    for res in residues :
      
      nextRes = res.nextResidue
      link = res.getFromLinkDict(res.currentSpinSystemAssigned, nextRes.currentSpinSystemAssigned)
      score += sum([1.0/pl.peak.degeneracy * pl.preMultipliedScore for pl in link.activePeakLinks]) #+ link.score  
          
    self.score = score
    
  def registerTextObserver(self,observerFunction):
    '''Register functions that get text information
       during the procedure.
       kwargs:  observerFunction:  function that is called
                                to send a text message.
                                This function receives
                                one argument which is
                                a string.
    '''
    
    self.textObservers.append(observerFunction)
    
  def unRegisterTextObserver(self,observerFunction):
    '''Unregister notifier function.'''
    
    self.textObservers.remove(observerFunction)
    
  def notifyTextObservers(self, text):
    '''Call all listener functions.'''
    for observerFunction in self.textObservers:
      
      observerFunction(text)
    
  def registerEnergyObserver(self,observerFunction):
    '''Register functions that get information
       about the energy during the annealling.
       kwargs:  observerFunction: function that is called
                                to send a text message.
                                This function receives
                                two arguments, the energy
                                and a which point in the
                                annealing it corresponds to.
    '''
    
    self.energyObservers.append(observerFunction)
    
  def unRegisterEnergyObserver(self,observerFunction):
    '''Unregister notifier function.'''
    
    self.energyObservers.remove(observerFunction)
    
  def notifyEnergyObservers(self, energy, pointNumber):
    '''Call all listener functions.'''
    for observerFunction in self.energyObservers:
      
      observerFunction(energy, pointNumber)
  
          