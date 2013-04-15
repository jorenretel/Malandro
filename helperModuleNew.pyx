# cython: profile=True

import math
import random
from numpy.random import randint,random_sample



import Tkinter
import re
import os
import random
import math
import time

import cProfile

import pickle


from memops.gui.MessageReporter import showWarning

from ccpnmr.analysis.core.ChemicalShiftBasic import getShiftsChainProbabilities

from ccpnmr.analysis.core.MoleculeBasic import getResidueCode
from ccpnmr.analysis.core.Util import stringFromExperimentSpectrum, getSpectrumPosContourColor

from ccp.util.LabeledMolecule import getIsotopomerSingleAtomFractions, getIsotopomerAtomPairFractions,  atomPairFrac,  atomPairFractions

from ccpnmr.analysis.core.Util import getAnalysisDataDim

from pythonStyleClasses import *

cdef extern from "math.h":

  double exp(double x)
  
  

  

cdef class autoAssign :
  
  
  cdef object chain
  
  cdef public myDataModel DataModel
  
  cdef object useAssignments
  
  cdef object useTentative
  
  cdef int amountOfRepeats
  
  cdef int amountOfSteps
  
  cdef object shiftList
  
  cdef object nmrProject
  
  cdef double minIsoFrac
  
  cdef list selectedSpectra
  
  cdef object project
  
  cdef str sourceName
  
  cdef int hc
  
  cdef object updateInfoText
  
  cdef object addEnergyPoint
  
  cdef list acceptanceConstantList
  
  cdef double score
  
  cdef object typeSpinSystems
  
  cdef object useDimenionalAssignments
  
  
  
  def __init__(self):
    
    self.hc = 10000
    
    
  def getResults(self):
    
    return self.getResultsC()


  cdef object getResultsC(self):
    
    return self.DataModel.pyDataModel


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
    
    self.selectedSpectra = connector.selectedSpectra
    
    self.sourceName = connector.sourceName
    
    self.updateInfoText = connector.updateInfoText
    
    self.acceptanceConstantList = connector.acceptanceConstantList
    
    self.addEnergyPoint = connector.addEnergyPoint
    
    self.typeSpinSystems = connector.typeSpinSystems
    
    self.useDimenionalAssignments = connector.useDimenionalAssignments
    
    
  def preCalculateDataModel(self)  :
    
    self.updateInfoText('Setting up model for calculations...')
    
    self.DataModel = myDataModel(self)
        
    self.updateInfoText('Setup-up all spectra...')
    
    self.DataModel.setupSpectra()

    self.DataModel.setupChain()

    self.createSpinSytemsAndResonances()
        
    self.updateInfoText('Simulating spectra...')
    
    self.simulateSpectra()
    
    self.updateInfoText('Evaluating possible dimensional contributions to peak in real spectra...')
    
    self.calculateAllPeakContributions()
    
    self.updateInfoText('Matching simulated with real spectra...')
        
    self.matchSimulatedWithRealSpectra()

    self.createJokerSpinSystems()
    
    self.updateInfoText('Scoring links between spin systems...')
    
    self.scoreAllLinks()
    
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
    
    cdef aResidue res
    
    cdef int i
    
    cdef aResidue resimin1
    
    cdef aResidue resiplus1
    
    cdef str ccpCode
    
    cdef int seqCode
    
    cdef list listWithSpinSystems
    
    cdef mySpinSystem spinSystem
    
    cdef list listWithFittingSpinSystems
    
    cdef mySpinSystem randomSpinSystem
    
    
    
    
    
    
    

    useAssignments = self.useAssignments
    useTentative = self.useTentative
    
    DataModel = self.DataModel
    
    sample = random.choice

    assignedSpinSystems = self.makePrivateCopyOfDictContainingLists(DataModel.previouslyAssignedSpinSystems)
    
    tentativeSpinSystems =  self.makePrivateCopyOfDictContainingLists(DataModel.tentativeSpinSystems)
    
    justTypedSpinSystems = self.makePrivateCopyOfDictContainingLists(DataModel.justTypedSpinSystems)
    
    allSpinSystems = self.makePrivateCopyOfDictContainingLists(DataModel.mySpinSystems)
    
    jokerSpinSystems = self.makePrivateCopyOfDictContainingLists(DataModel.jokerSpinSystems)
    
    
    if useAssignments and useTentative:
      
      dictio = self.mergeDictionariesContainingLists([justTypedSpinSystems,jokerSpinSystems])  
      
    elif useAssignments :
    
      dictio = self.mergeDictionariesContainingLists([justTypedSpinSystems, tentativeSpinSystems,jokerSpinSystems])
                                              
    elif useTentative :
      
      dictio = self.mergeDictionariesContainingLists([justTypedSpinSystems,jokerSpinSystems])
                                                
    else :
      
      dictio = self.makePrivateCopyOfDictContainingLists(allSpinSystems)                                                            
    
    
    
    
    
    i = 0
    
    for res in DataModel.myChain.residues :
      #print i + 1
      
      isAssigned = False
      
      if i > 0 :
        resimin1 = DataModel.myChain.residues[i-1]
      else :
        resimin1 = None         
      if i < (len(DataModel.myChain.residues)-1) :
        resiplus1 = DataModel.myChain.residues[i+1]
      else :
        resiplus1 = None
      i = i + 1  
        
      
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
          
          listWithFittingSpinSystems = dictio[ccpCode]
          
              
        if len(listWithFittingSpinSystems) > 0 :      
          
                      
          randomSpinSystem = sample( listWithFittingSpinSystems)                                                        #listWithFittingSpinSystems[random.randint(0, (len(listWithFittingSpinSystems)-1))] 
             
          randomSpinSystem.currentResidueAssignment = res
          
          res.currentSpinSystemAssigned = randomSpinSystem 
          
          isAssigned = True
          
          if randomSpinSystem in listWithFittingSpinSystems :
          
            listWithFittingSpinSystems.remove(randomSpinSystem)
            
            
          
      if not isAssigned :
          
        print 'something went wrong during random assignment at the start of the procedure, not all residues have an assignment'
              

  def runAnnealling(self):
    

    useAssignments = self.useAssignments
    useTentative = self.useTentative    

    amountOfStepsPerTemperature = self.amountOfSteps
        
    DataModel = self.DataModel
            
    exp = math.exp
    cutoff = random.random
    sample = random.choice
    randint = random.randint
        
    AcceptanceConstantList = self.acceptanceConstantList
    
    justTypedSpinSystems =DataModel.justTypedSpinSystems
    allSpinSystems = DataModel.mySpinSystems
    tentativeSpinSystems = DataModel.tentativeSpinSystems
    allSpinSystemsWithoutAssigned = DataModel.allSpinSystemsWithoutAssigned

    
    if useAssignments :                                                 # If the already made assignments are used, the dictionary will just include unassigned spinsystems and Joker spinsystems
      dict = allSpinSystemsWithoutAssigned
    else :                                                                       # If the already assigned are chosen to still be permutated, the dictionary will also include the ready assinged spinsystems
      dict = allSpinSystems
      
    listWithSpinSystems = []  
      
    for key in dict :
      
      listWithSpinSystems = listWithSpinSystems + (dict[key])
      
    listWithSpinSystems = list(set(listWithSpinSystems))
      
    score = 0
    
    for x,AcceptanceConstant in enumerate(AcceptanceConstantList) :
      
      self.annealingSub(AcceptanceConstant,amountOfStepsPerTemperature,listWithSpinSystems)
      
      self.addEnergyPoint(self.score,x)


  def testB(self,string):
    
    self.updateInfoText(string)
    
    
  def startMonteCarlo(self):
    
    cdef int repeat
    
    cdef aResidue res
    
    cdef str info
    
    repeat = self.amountOfRepeats
    
    info = 'Running annealing number %s out of ' + str(repeat) + '...'
    
    self.setupSpinSystemExchange()
    
    for x in range(repeat) :
      
      self.cleanAssignments()
    
      self.doRandomAssignment()
      
      self.setupPeakInformationForRandomAssignment()
      
      self.scoreInitialAssignment()
      
      self.updateInfoText(info %str(x+1))
      
      self.runAnnealling()
      
      i = 1
      matches = 0

      for res in self.DataModel.myChain.residues :
        
        res.solutions.append(res.currentSpinSystemAssigned)
        
        res.currentSpinSystemAssigned.solutions.append(res)
        
        
        if res.ccpnResidue.findFirstResonanceGroup() :

          if res.currentSpinSystemAssigned.spinSystemNumber == res.ccpnResidue.findFirstResonanceGroup().serial :

            matches = matches + 1
        i = i + 1

    self.convertToPythonStyleDataModel()
    
    self.updateInfoText('Done')
    
      
  cdef void cleanAssignments(self) :
    
    cdef myDataModel DataModel
    
    cdef myChain chain
    
    cdef list residues
    
    cdef list spectra
    
    cdef aSpectrum spectrum
    
    cdef dict mySpinSystems
    
    cdef dict tentativeSpinSystems
    
    cdef str key
    
    cdef list spinSystemList
    
    cdef mySpinSystem spinSys
    
    cdef aPeak peak
    
    cdef aResidue res
    
    
    
    DataModel = self.DataModel
    
    chain = DataModel.myChain
    
    residues = chain.residues
    
    spectra = DataModel.spectra
    
    mySpinSystems = DataModel.mySpinSystems
    
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
    
    cdef myDataModel DataModel
    
    cdef myChain chain
    
    cdef list residues
    
    cdef list resWithoutLast
    
    cdef aResidue residue
    
    cdef mySpinSystem spinSysI
    
    cdef mySpinSystem spinSysIplus1
    
    cdef spinSystemLink link
    
    cdef aPeak peak
    
    cdef int hc
    
    hc = self.hc
    
    DataModel = self.DataModel
    
    chain = DataModel.myChain
    
    residues = chain.residues
    
    
    resWithoutLast =  residues[0:len(residues)-1]
    
    for residue in resWithoutLast :
      
      spinSysI = residue.currentSpinSystemAssigned
      spinSysIplus1 = residue.nextResidue.currentSpinSystemAssigned
      
      if not spinSysI.isJoker and not spinSysIplus1.isJoker :
      
        link = residue.linkDict[spinSysI.spinSystemNumber*hc+spinSysIplus1.spinSystemNumber]
      
        for peak in link.realPeaks :
        
          peak.degeneracy += 1
        

  cdef dict makePrivateCopyOfDictContainingLists(self,  dict dictio):                                                   # Make a copy of the spinsystemDict that can be modified without changing the original dictionary
    
    cdef dict copiedDict
    cdef list copyOfList
    
    copiedDict = {}
    
    for key in dictio :
      
      copyOfList = dictio[key][:]                                 #a copy of each list is made by slicing the list, trick
      copiedDict[key] = copyOfList
      
    return copiedDict
      
      
  cdef dict mergeDictionariesContainingLists(self, list dictionaries):
    
    
    cdef list dicts
    cdef dict dictio
    cdef dict newDict
    cdef list newList
    
    
    keys = []
    
    dicts = []
      
    for dictio in dictionaries :
      
      dicts.append(self.makePrivateCopyOfDictContainingLists(dictio))
      
      keys = keys + dictio.keys()
      
    keys = set(keys)
    
    newDict = {}
    
    for key in keys :
      
      newList = []
      
      for dictio in dicts :
        
        if key in dictio :
          
          newList += dictio[key]
          
  
      newDict[key] = list(set(newList))
      
    return newDict  
          
          
  cdef void createJokerSpinSystems(self):
    
    '''
    When there are less spin systems that are typed to a
    certain amino acid type than there are residues of that
    type in the sequence, some Jokers have to be introduced. 
    
    '''
    
    cdef myDataModel DataModel
    
    cdef str key
    
    cdef int amountOfAssignedSpinsystems
    
    cdef int amountOfTypedSpinSystems
    
    cdef int amountOfTentativeSpinSystemsWithOnlyOneCcpCode
    
    cdef mySpinSystem spinsys
    
    cdef int amountOfResiduesOfThisType
    
    cdef int short
    
    cdef int i
    
    cdef mySpinSystem newSpinSystem
    
    DataModel = self.DataModel    
    
    
    i = 1
    
    for key in DataModel.mySpinSystems.keys() :

      
      amountOfAssignedSpinsystems = 0
      
      if key in DataModel.previouslyAssignedSpinSystems :
      
        amountOfAssignedSpinsystems = len(DataModel.previouslyAssignedSpinSystems[key])
        
        
      amountOfTypedSpinSystems = 0
      
      
      if key in DataModel.justTypedSpinSystems :
        
        amountOfTypedSpinSystems = len(DataModel.justTypedSpinSystems[key])
        
      
      amountOfTentativeSpinSystemsWithOnlyOneCcpCode = 0
      
#      if key in DataModel.tentativeSpinSystems :                                                         # Commented this out since it gives problems when there is more than one tentative spin system pointing to the same residue and nothing else. In this case too little Jokers would be added. Rather have a Joker too much than too little.
#      
#        for spinsys in DataModel.tentativeSpinSystems[key] :
#          
#          if len(spinsys.tentativeCcpCodes) == 1 :
#            
#            amountOfTentativeSpinSystemsWithOnlyOneCcpCode += 1
          
  
      
      if key in DataModel.residueTypeFrequencyDict :
        
        amountOfResiduesOfThisType = DataModel.residueTypeFrequencyDict[key]
        
      else :
      
        amountOfResiduesOfThisType = 0  
      
      short = amountOfResiduesOfThisType- (amountOfAssignedSpinsystems + amountOfTypedSpinSystems + amountOfTentativeSpinSystemsWithOnlyOneCcpCode)
      
      
      
      if short < 0 :
        
        string = 'You seem to have more ' + key + ' spin systems than residues of this type are in the sequence'

        short = 0
      
      for x in range(short) :
        

        
        newSpinSystem = mySpinSystem()
        
        newSpinSystem.isJoker = True
        
        newSpinSystem.ccpCode = key
        
        newSpinSystem.spinSystemNumber = i * 1000000
        
        i = i + 1
        

        
        if key in DataModel.mySpinSystems :
        
          DataModel.mySpinSystems[key].append(newSpinSystem)
          
        else :
          
          DataModel.mySpinSystems[key] = [newSpinSystem]
          
          
        if key in DataModel.jokerSpinSystems:
          
          DataModel.jokerSpinSystems[key].append(newSpinSystem)
          
        else :
          
          DataModel.jokerSpinSystems[key] = [newSpinSystem]
          
          
        if key in DataModel.allSpinSystemsWithoutAssigned:
          
          DataModel.allSpinSystemsWithoutAssigned[key].append(newSpinSystem)
          
        else :
          
          DataModel.allSpinSystemsWithoutAssigned[key] = [newSpinSystem]
          



        
        
        
      
        
    
    string = str(i-1) + ' joker spinsystems are used in this calculation.'    

  cdef void simulateSpectra(self) :
    
    cdef myDataModel DataModel
    
    cdef list residues
    
    cdef aSpectrum spectrum
    
    cdef object ccpnSpectrum
    
    cdef str refExpName
    
    cdef object scheme
    
    cdef list simulatedPeakMatrix
    
    cdef int i
    
    cdef aResidue resA
    
    cdef aResidue resB
    
    cdef set simpleAtomSites
    
    cdef object refExperiment
    
    cdef object PT
    
    cdef object expGraph
    
    cdef dict recordedExpMeasuments
    
    cdef list expTransfers
    
    cdef object refExpDim
    
    cdef list expDimensions
    
    cdef list expMeasurementPathWay
    
    cdef list atomSitePathWay
    
    cdef object expStep
    
    #cdef int i
    
    cdef int ii
    
    
    #cdef aResidue resA
    
    #cdef aResidue resB
    
    cdef aResidue res
    
    cdef list atomPathWayResA
    
    cdef list atomPathWayResB
    
    cdef list atomPathWayResX
    
    cdef object atomSite
    
    cdef list atomsForThisAtomSite
    
    cdef str atomSiteName
    
    cdef anAtom atom
    
    cdef list resPathWays
    
    cdef list resPathWaysTemp
    
    cdef list resPathWay
    
    cdef int resID 
    
    cdef list atomGroupPathWays
    
    cdef list atomPathWays
    
    cdef list atomPathWay
    
    cdef list atomGroup
    
    cdef list tempList
    
    cdef object isPossible
    
    cdef object transfer
    
    cdef anAtom atomA
    
    cdef anAtom atomB
    
    cdef frozenset bondSetA
    
    cdef frozenset bondSetB
    
    cdef int resIDA
    
    cdef int resIDB
    
    cdef list peaks
    
    cdef double minIsoFrac
    
    cdef simulatedPeakContrib contrib
    
    
    
    
    minIsoFrac = self.minIsoFrac
    
    
    
    
    simpleAtomSites = set(['H', 'HA','HB','N', 'CA', 'CB'])              #TODO: actually H stands for all H's here, not just the backbone amide H. And check HB out, actually not one atom either. Gly actually has two HAs.
    
    carbons = set(['C', 'Cali','Cmet', 'Caro'])
    
    DataModel = self.DataModel
    
    residues = DataModel.myChain.residues 
    
    for spectrum in DataModel.spectra:                 
      
          
      ccpnSpectrum = spectrum.ccpnSpectrum
    
      refExpName = spectrum.ccpnSpectrum.experiment.refExperiment.name
      
      
      self.updateInfoText('Simulating ' + spectrum.name)
      
      
      scheme = spectrum.labellingScheme
      
      simulatedPeakMatrix = spectrum.simulatedPeakMatrix
    
      refExperiment = ccpnSpectrum.experiment.refExperiment
      
      molLabelFractions = ccpnSpectrum.experiment.findFirstLabeledMixture().molLabelFractions
      
      PT =  refExperiment.nmrExpPrototype
    
      expGraph = PT.findFirstExpGraph()    
      
      recordedExpMeasuments = {}
      
      expTransfers = expGraph.sortedExpTransfers()
      
      
      # Looking which expMeasurements where actually done
      
      for refExpDim in refExperiment.sortedRefExpDims() :
        
        recordedExpMeasuments[refExpDim.findFirstRefExpDimRef().expMeasurement] = refExpDim.dim
        
        
      
      expDimensions = []
      expMeasurementPathWay = []
      
      atomSitePathWay = []
      
      # Making a list of atomSites that are visited during the experiment and another list of the same length that marks whether they were recorded and if so, to which dimension
      
      for expStep in expGraph.sortedExpSteps() :                                                   # TODO: This is not good with in and out experiments, rather take atomSites from ExpTransfers
        
        atomSitePathWay.append(expStep.expMeasurement.findFirstAtomSite())
        
        if expStep.expMeasurement in recordedExpMeasuments.keys() :
          
          expDimensions.append(recordedExpMeasuments[expStep.expMeasurement])
          
        else :
          
          expDimensions.append(None)
          
      
    
      # Create lists with all combinations of 1 and 2 (1 indicating first residue, 2 indicating second) of a the amount of different atomSites visited. Will check out later which combinations actually work.
      resPathWays = self.createPermutations(len(expGraph.sortedExpTransfers())+1)
      
      resPathWaysTemp = []
      for resPathWay in resPathWays :
        
        if len(set(resPathWay)) == 2 :                            # We are not interested in peaks that do not connect the two residues in any way.
          
          resPathWaysTemp.append(resPathWay)
          
      resPathWays =  resPathWaysTemp   
          
      for i, resA in enumerate(residues[:-1]) :
    
        peaks = []
          
        resB = residues[i+1]
        
        atomPathWayResA = []
        atomPathWayResB = []
        
        possibleAtomPathWays = []
        
        # Transforming atomSides into atoms for both residues in the sequential pair of resA and resB, thereby filling the lists atomPathWayResA and atomPathWayResB
     
        for atomPathWayResX,  res in zip([atomPathWayResA, atomPathWayResB], [resA, resB]) :
        
          for atomSite in atomSitePathWay :
            
            atomsForThisAtomSite = []
          
            atomSiteName = atomSite.name
            
            if atomSiteName in simpleAtomSites and atomSiteName in res.atomsByName :
              
              atomsForThisAtomSite.append(res.atomsByName[atomSiteName])
              
            elif atomSiteName == 'CO' and 'C' in res.atomsByName :
              
              atomsForThisAtomSite.append(res.atomsByName['C'])
              
            elif atomSiteName in carbons :                                                      #TODO, make distinction between these
              
              for atom in res.atoms :
                
                if atom.ccpnAtom.chemAtom.elementSymbol == 'C':
                  
                  atomsForThisAtomSite.append(atom)
                  
                  
            atomPathWayResX.append(atomsForThisAtomSite)
          
          
        
        # Going through all generated possibilities of a sequence of atomSites being located in either resA or resB. For an experiment with three steps (like HNCA), this would be [[1,1,1],[1,1,2]......[2,2,2]]
        for resPathWay in resPathWays :
          
          atomPathWays = [[]]
          
          atomGroupPathWays = []

          #Lining up the different atoms that can be visited in a nested list
          for ii,  resID in enumerate(resPathWay) :               # looping through sequence of 1s and 2s
          
            if resID == 1 :
              
              atomGroupPathWays.append(atomPathWayResA[ii])
              
            elif resID == 2 :
              
              atomGroupPathWays.append(atomPathWayResB[ii])
          

          #Making all possible routes between different atoms that could be there (when not considering the type of transfers yet).    
          for atomGroup in atomGroupPathWays :
        
            tempList = []
        
            for atom in atomGroup :
          
              for atomPathWay in atomPathWays :
                
                tempList.append(atomPathWay + [atom])
                
            atomPathWays = tempList
            

          # Now we're going to check which of these pathways can actually really happen when we take the type of transfer (for instance 'one bond') into account.  
          for atomPathWay in atomPathWays :
            
            isPossible = True

            # Therefore we'll look whether every transfer is possible
            for ii, transfer in enumerate(expTransfers):
              
              atomA = atomPathWay[ii]
              atomB = atomPathWay[ii+1]
              
              resIDA = resPathWay[ii]
              resIDB = resPathWay[ii+1]
              
              
              if (transfer.transferType == 'onebond' or transfer.transferType == 'Jcoupling')  :
                
                
                # Over the amide bond
                if (atomA.atomName == 'N' and atomB.atomName == 'C' and resIDA == 2 and resIDB == 1) or (atomA.atomName == 'C' and atomB.atomName == 'N' and resIDA == 1 and resIDB == 2) :
              
                  pass
                  
                # Within the same residue and directly bonded
                elif resIDA == resIDB and atomA.ccpnAtom.chemAtom.getChemBonds().intersection(atomB.ccpnAtom.chemAtom.getChemBonds()):
                  
                  pass
                  
                # In all other cases the transfer pathway is not possible  
                else :
                  
                  isPossible = False
                  
                  break
                  


            if isPossible :
           
              possibleAtomPathWays.append(atomPathWay)
              
              
              
        importantAtoms = set()
        for atomPathWay in  possibleAtomPathWays :   

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
                
                atomLabels = isotopomer.findAllAtomLabels(name=atomName, subType=subType)
                
                atom.labelInfoTemp[isotopomer] = atomLabels
                
                    
                    
        for atomPathWay in possibleAtomPathWays :
          
          colabelling = self.getCoLabellingFractionNew(spectrum,atomPathWay)   


          if colabelling > minIsoFrac :
              
            newPeak = simulatedPeak()
            newPeak.colabelling = colabelling
            newPeak.spectrum = spectrum
            
            for dimension,  atom in zip(expDimensions, atomPathWay) :
              
              if dimension :                                                              # If None, the dimension was not recorded
            
                contrib = simulatedPeakContrib()
            
                contrib.ccpCode = atom.residue.ccpCode
                
                contrib.atomName = atom.atomName
                
                if atom.residue == resA :
                
                  contrib.firstOrSecondResidue = 1
                  
                elif atom.residue == resB :
                
                  contrib.firstOrSecondResidue = 2  
            
                contrib.dimNumber = dimension
            
                newPeak.simulatedPeakContribs.append(contrib)
            
            peaks.append(newPeak)
            
        simulatedPeakMatrix.append(peaks)    

     
  cdef list createPermutations(self, int length):
    
    '''This method creates all permutations of a given length of
       ones and twos. Output for length=3, will be a nested list like
       [ [1,1,1], [2,1,1], [2,2,1] .... etc. ]. Used as a hypothetical
       series of magnetization transfers between two residues with a
       certain amount of transfer steps (length). Later it will be
       checked whether the transfers are actually possible given the
       experiment.'''
       
    cdef list lists
    cdef list newLists
    cdef list newList1
    cdef list newList2
    cdef list oneList
    cdef int i
    
    lists = [[]]

    for i in range(length):
      
      newLists = []
      
      for oneList in lists :
        
        newList1 = oneList[:]
        newList2 = oneList[:]
        newList1.append(1)
        newList2.append(2)
        
        newLists.append(newList1)
        newLists.append(newList2)
        
      lists = newLists
      
    return lists
  
    
  cdef void createSpinSytemsAndResonances(self):
    
    
    print 'Creating an initial setup for the calculations......'


    cdef int cc
    
    cdef int a
    
    cdef object shiftList
    cdef object ccpnChain
    
    cdef mySpinSystem newSpinSystem
    
    cdef object resonanceGroup
    
    cc = 0
    
    shiftList = self.shiftList
    ccpnChain = self.DataModel.myChain.ccpnChain
    
    
    
    mySpinSystems = self.DataModel.mySpinSystems
      
    previouslyAssignedSpinSystems = self.DataModel.previouslyAssignedSpinSystems
        
    justTypedSpinSystems = self.DataModel.justTypedSpinSystems 
    
    tentativeSpinSystems = self.DataModel.tentativeSpinSystems 
    
    untypedSpinSystems = self.DataModel.untypedSpinSystems
    
    allSpinSystemsWithoutAssigned = self.DataModel.allSpinSystemsWithoutAssigned
    
      
    for resonanceGroup in self.nmrProject.resonanceGroups :                                                                     # taking all spinsystems in the project
    
      SpinSysType = None
          
      
      if resonanceGroup.resonances : 

        
        if resonanceGroup.residue and resonanceGroup.residue.chain == ccpnChain :                                  # SpinSystem is assigned to a residue in the selected chain
        
          newSpinSystem = mySpinSystem()                                                                                                      # creating a spinsystemobject for myself, to store some things in without modifying the project

          newSpinSystem.spinSystemNumber = resonanceGroup.serial
          
          newSpinSystem.ccpnResonanceGroup = resonanceGroup
        
          SpinSysType = 'assigned'
          
          newSpinSystem.ccpCode = getResidueCode(resonanceGroup.residue.molResidue)

          newSpinSystem.ccpnSeqCode = int(resonanceGroup.residue.seqCode)

          cc = cc +1

          
          
          
          
          
        elif resonanceGroup.residueProbs :                                                                                                      # SpinSystem has one or more tentative assignments. Got this piece out of EditSpinSystem.py in popups.
          
          newSpinSystem = mySpinSystem()                                                                                                      # creating a spinsystemobject for myself, to store some things in without modifying the project

          newSpinSystem.spinSystemNumber = resonanceGroup.serial
          
          newSpinSystem.ccpnResonanceGroup = resonanceGroup
          
          
          SpinSysType = 'tentative'
          

          

          cc = cc +1
          for residueProb in resonanceGroup.residueProbs:
            if not residueProb.weight:
              continue
              
            residue = residueProb.possibility
            
            seq = residue.seqCode
            resCode = residue.ccpCode #getResidueCode(residue)
            
            if residue.chain == ccpnChain :

              newSpinSystem.tentativeCcpCodes.append(resCode)                                                                   # Only consider the tentative assignments to residues in the selected chain.                                                         
              newSpinSystem.tentativeSeqCodes.append(seq)

        elif resonanceGroup.ccpCode :                                                                                                              # Residue is just Typed
        
          newSpinSystem = mySpinSystem()                                                                                                      # creating a spinsystemobject for myself, to store some things in without modifying the project

          newSpinSystem.spinSystemNumber = resonanceGroup.serial
          
          newSpinSystem.ccpnResonanceGroup = resonanceGroup
        
          SpinSysType = 'justTyped'
        
          cc = cc +1
          
          newSpinSystem.ccpCode = resonanceGroup.ccpCode
          
          
          
        elif self.typeSpinSystems :                                                                                                                                                      # For spin systems that are not typed at all, I want to type them to one or more amino acid types here on the fly, later the user can decide whether these are actually used. 
        
        
          newSpinSystem = mySpinSystem()                                                                                                      # creating a spinsystemobject for myself, to store some things in without modifying the project
        
          newSpinSystem.spinSystemNumber = resonanceGroup.serial
          
          newSpinSystem.ccpnResonanceGroup = resonanceGroup
        
          SpinSysType = 'unTyped'
        
        
          shifts = []
          for resonance in resonanceGroup.resonances:
            if resonance.isotopeCode in ('1H',  '13C',  '15N') :
              shift = resonance.findFirstShift(parentList=shiftList)
              if shift:
                shifts.append(shift)
        
          scores = getShiftsChainProbabilities(shifts, ccpnChain)
          total = sum(scores.values())
          
          print '++++++++++++++++++++++++++'
          
          print scores
          
          
          scoreList = []
          
          scoreDict = {}
          
          if total:
              
            for ccpCode, score in scores.items() :
              
              if score > float(self.DataModel.residueTypeFrequencyDict[ccpCode])/len(self.DataModel.myChain.residues) :
                
                scoreDict[ccpCode] = score
              
        
          
          
          print ' ------------------- '
          print 'serial: '
          print resonanceGroup.serial
          
          print scoreDict
          
          

          
          
          
          
          
          
          
          

        a = 0 
          
        for resonance in resonanceGroup.resonances :
          

          a = a +1
          
          if len(resonance.assignNames) > 0 :
            
          
            newResonance = myResonance()
          
            newResonance.mySpinSystem = newSpinSystem
            
            newResonance.ccpnResonance = resonance
            
            assignName = resonance.assignNames[0]
            
            if assignName == 'CE*':                                 #taking out of some chirality stuff that makes everything harder
              assignName = 'CE1'
            elif assignName == 'CD*':
              assignName = 'CD1'
            elif assignName == 'CG*':
              assignName = 'CG1'
            
            newResonance.atomName = assignName
            
            
            newResonance.isotope = resonance.isotopeCode
            
            newResonance.CS = resonance.findFirstShift().value
            
            newSpinSystem.resonanceDict[assignName] = newResonance                  # adding the resonance to the dictionary in the spinsystem
            
            if assignName == 'N' :                                                                              # Also adding it to a bunch of lists from quick retrieval
              
              self.DataModel.Nresonances.append(newResonance)
              
            elif assignName == 'CA' :
              
              self.DataModel.CAresonances.append(newResonance)
              
            elif assignName == 'CB' :
            
              self.DataModel.CBresonances.append(newResonance)  
              
            elif assignName == 'C' :
              
              self.DataModel.COresonances.append(newResonance)
              
            elif assignName == 'H' :
              
              self.DataModel.Hresonances.append(newResonance)

            if resonance.isotopeCode == '13C' :
              
              self.DataModel.CXresonances.append(newResonance)
              
            elif resonance.isotopeCode == '1H' :
            
              self.DataModel.HXresonances.append(newResonance)
            
            elif resonance.isotopeCode == '15N' :
            
              self.DataModel.NXresonances.append(newResonance)
              
            
            
            
            
            
            
        if SpinSysType== 'assigned' :
          
          if newSpinSystem.ccpCode in mySpinSystems :       
            
            mySpinSystems[newSpinSystem.ccpCode].append(newSpinSystem)
            
          else :                                                                                                                          
            
            mySpinSystems[newSpinSystem.ccpCode] = [newSpinSystem]
            
            
          if newSpinSystem.ccpCode in previouslyAssignedSpinSystems :       
            
            previouslyAssignedSpinSystems[newSpinSystem.ccpCode].append(newSpinSystem)
            
          else :                                                                                                                          
            
            previouslyAssignedSpinSystems[newSpinSystem.ccpCode] = [newSpinSystem]
            
            
          
        elif SpinSysType== 'tentative' and len(newSpinSystem.tentativeSeqCodes) > 0 :
          
          for ccpCode in list(set(newSpinSystem.tentativeCcpCodes))  :

            if ccpCode in mySpinSystems :       
              
              mySpinSystems[ccpCode].append(newSpinSystem)
              
            else :                                                                                                                          
              
              mySpinSystems[ccpCode] = [newSpinSystem]
              
              
            if newSpinSystem.ccpCode in allSpinSystemsWithoutAssigned :       
              
              allSpinSystemsWithoutAssigned[newSpinSystem.ccpCode].append(newSpinSystem)
              
            else :                                                                                                                          
              
              allSpinSystemsWithoutAssigned[newSpinSystem.ccpCode] = [newSpinSystem]
              
              
            if ccpCode in tentativeSpinSystems :       
              
              tentativeSpinSystems[ccpCode].append(newSpinSystem)
              
            else :                                                                                                                          
              
              tentativeSpinSystems[ccpCode] = [newSpinSystem]
              
              
              

            
          
        elif SpinSysType== 'justTyped' :
          
          if newSpinSystem.ccpCode in mySpinSystems :       
            
            mySpinSystems[newSpinSystem.ccpCode].append(newSpinSystem)
            
          else :                                                                                                                          
            
            mySpinSystems[newSpinSystem.ccpCode] = [newSpinSystem]
            
            
          if newSpinSystem.ccpCode in allSpinSystemsWithoutAssigned :       
            
            allSpinSystemsWithoutAssigned[newSpinSystem.ccpCode].append(newSpinSystem)
            
          else :                                                                                                                          
            
            allSpinSystemsWithoutAssigned[newSpinSystem.ccpCode] = [newSpinSystem]
            
            
          if newSpinSystem.ccpCode in justTypedSpinSystems :       
            
            justTypedSpinSystems[newSpinSystem.ccpCode].append(newSpinSystem)
            
          else :                                                                                                                          
            
            justTypedSpinSystems[newSpinSystem.ccpCode] = [newSpinSystem]
            
            
            

          
        elif SpinSysType== 'unTyped' :
          
          newSpinSystem.aminoAcidProbs = scoreDict
            
          for ccpCode in scoreDict.keys() :

            if ccpCode in mySpinSystems :       
              
              mySpinSystems[ccpCode].append(newSpinSystem)
              
            else :                                                                                                                          
              
              mySpinSystems[ccpCode] = [newSpinSystem]
              
              
            if newSpinSystem.ccpCode in allSpinSystemsWithoutAssigned :       
              
              allSpinSystemsWithoutAssigned[newSpinSystem.ccpCode].append(newSpinSystem)
              
            else :                                                                                                                          
              
              allSpinSystemsWithoutAssigned[newSpinSystem.ccpCode] = [newSpinSystem]
              

  cdef void scoreAllLinks(self):
    
    cdef list residues 
    
    cdef aResidue res
    
    cdef dict linkDict
    
    cdef spinSystemLink linkObject
    
    cdef int NfoundPeaks
    
    cdef simulatedPeak peak
    
    cdef simulatedPeakContrib contrib
    
    cdef str atomName
    
    cdef list atomsUsed
    
    cdef int NatomsUsed
    
    cdef int hc
    
    hc = self.hc
    
      
    DataModel = self.DataModel

    residues = DataModel.myChain.residues
    
    for res in residues :
      
      linkDict = res.linkDict
      

      
      
      for key,  linkObject in linkDict.items() :
        
        
        NfoundPeaks = len(linkObject.realPeaks)                                                                                          # These are the simulated peaks where a real peak was found for in the spectra during the matching procedure

        atomsUsed = []
        
        for peak in linkObject.simulatedPeaks :                                                                                                              # These are the simulated variant of the 'realPeaks'.
            
          for contrib in peak.simulatedPeakContribs :
                
            atomName = contrib.atomName
            
            if atomName not in atomsUsed :
                
                atomsUsed.append(atomName)
                
        NatomsUsed = len(atomsUsed)
        
        linkObject.score = NfoundPeaks + NatomsUsed                     #- NfoundPeaksThatShouldNotHaveBeenThere


  cdef void matchSimulatedWithRealSpectra(self):
    
    self.updateInfoText('Matching real to simulated spectra.....')
    
    cdef aSpectrum spectrum
    
    DataModel = self.DataModel
    
    for spectrum in DataModel.spectra :
      
      self.updateInfoText('Match simulated with real spectrum: ' + spectrum.name)
    
      self.matchSpectrum(spectrum)


  cdef void matchSpectrum(self, aSpectrum spectrum):
    
    cdef myDataModel DataModel
    
    cdef list simulatedPeakMatrix
    
    cdef list residues
    
    cdef list simulatedPeakList
    
    cdef aResidue resA
    
    cdef aResidue resB
    
    cdef dict allSpinSystems
    
    cdef list spinsystemsaa1
    
    cdef list spinsystemsaa2
    
    cdef mySpinSystem spinSys1
    
    cdef mySpinSystem spinSys2
    
    cdef int presentPeaks 
    
    cdef list listWithPresentPeaks
    cdef list listWithSimulatedPeaks
    cdef list listWithNotFoundSimulatedPeaks
    
    cdef simulatedPeak simulatedPeak
    
    cdef list contributions
    
    cdef list resonances
    
    cdef simulatedPeakContrib contrib
    
    cdef list peakLists
    
    cdef list peaksInWindow
    
    cdef list rootOfSquaresList
    
    cdef list deltaCSsquaredList
    
    cdef aDimension dim
    
    cdef spinSystemLink linkObject
    
    cdef int i
    
    cdef double smallest
    
    cdef double x
    
    cdef int smallestIndex
    
    cdef myResonance resonance
    
    cdef aPeak peak
    
    cdef int hc
    

    
    
    
    
    
    
    hc = self.hc

    DataModel = self.DataModel
    simulatedPeakMatrix = spectrum.simulatedPeakMatrix
    
    allSpinSystems = DataModel.mySpinSystems
    
    residues = DataModel.myChain.residues

    for i,  simulatedPeakList in enumerate(simulatedPeakMatrix) :
      
      #print 'pppppppppppppppppppppppppppp'
      #print i
      #print simulatedPeakList
      resA = residues[i]
      resB = residues[i+1]
      
      ccpCodeA = resA.ccpCode
      ccpCodeB = resB.ccpCode
      

      
      spinsystemsaa1 = allSpinSystems[ccpCodeA]
      spinsystemsaa2 = allSpinSystems[ccpCodeB]
      
      
      
      count = 0
      
      for spinSys1 in spinsystemsaa1 :
        
        count = count + 1
        
        for spinSys2 in spinsystemsaa2 :
          
          presentPeaks = 0
          listWithPresentPeaks = []
          listWithSimulatedPeaks = []
          listWithNotFoundSimulatedPeaks = []
  
      
          for simulatedPeak in simulatedPeakList :
            
            contributions = simulatedPeak.simulatedPeakContribs
            

            
            resonances = []
            
            for contrib in contributions :
              
              
              if contrib.firstOrSecondResidue == 1 and contrib.atomName in spinSys1.resonanceDict :
                
                resonances.append(spinSys1.resonanceDict[contrib.atomName])
                
                  
              elif contrib.firstOrSecondResidue == 2 and contrib.atomName in spinSys2.resonanceDict :
                
                resonances.append(spinSys2.resonanceDict[contrib.atomName])
                  

              
            
            
            
            
            if len(resonances) == len(contributions) :                                                      # If the involved spinSystems have all atoms assigned that contribute to the simualated peak
              
              peakLists = []
              
              for resonance,  contrib in zip(resonances , contributions) :                          # Gather for each dimension in the spectrum the peaks that are on the right frequency (as determined in 'findPeakContributions')
                

                if spectrum.name in resonance.peakDimsLib and contrib.dimNumber in resonance.peakDimsLib[spectrum.name] :

                  peakLists.append( resonance.peakDimsLib[spectrum.name][contrib.dimNumber] )
                  
              
              if len(peakLists) == len(resonances) :                                                            # If there are 1 or more peaks for each dimension, we can go on...
                
                peaksInWindow = self.commonElementInLists(peakLists)                                      # And check whether 1 or more peaks that fit in one dimension actually also fit in all other dimensions. In that case the peak is in the multidimensional tolerance window
                
              else :
                
                peaksInWindow = []

              
              if peaksInWindow :
                
                presentPeaks = presentPeaks + 1
                
                if len(peaksInWindow) == 1 :                                                                            # If there is only one peak that falls within the window, we don't have to choose.
                
                  listWithPresentPeaks.append(peaksInWindow[0])
                  listWithSimulatedPeaks.append(simulatedPeak)
                
                else :                                                                                                                 # There might be more than on peak within the tolerances in all dimensions. We will choose the closest one 
                  
                  rootOfSquaresList = []
                  smallestIndex = 0
                  
                  for peak in peaksInWindow :
                    
                    deltaCSsquaredList = []

                    for dim in peak.dimensions :
                      
                      for resonance,  contrib in zip(resonances , contributions) :
                      
                        if dim.dimNumber == contrib.dimNumber :
                        
                          deltaCSsquaredList.append((dim.ppmValue - resonance.CS)**2)
                          
                    if len(deltaCSsquaredList) == len(contributions) :     
                    
                      rootOfSquaresList.append((sum(deltaCSsquaredList))**0.5)
                      
                    else :
                      
                      print 'something went wrong with the peak dimensions'
                      
                    #i = 0
                    smallest = 10000  
                    
                    for i, x in enumerate(rootOfSquaresList) :
                      
                      if x < smallest :
                        
                        smallest = x
                        smallestIndex = i
                        
                      #i += 1
                      
                      
                  
                  listWithPresentPeaks.append(peaksInWindow[smallestIndex])                                             # The peak with the smallest deviation is added to the list of found peaks
                  
                  listWithSimulatedPeaks.append(simulatedPeak)
             
              else :
               
                listWithNotFoundSimulatedPeaks.append(simulatedPeak)   
                  
          if (spinSys1.spinSystemNumber*hc+spinSys2.spinSystemNumber) in resA.linkDict :
     
            linkObject = resA.linkDict[(spinSys1.spinSystemNumber*hc+spinSys2.spinSystemNumber)] 
            
            linkObject.score = linkObject.score + presentPeaks
            linkObject.simulatedPeaks = linkObject.simulatedPeaks + listWithSimulatedPeaks
            linkObject.realPeaks = linkObject.realPeaks + listWithPresentPeaks
            linkObject.notFoundSimulatedPeaks += listWithNotFoundSimulatedPeaks
              
          else :                                                                                                                        # If there is no linkObject yet, make one
            
            linkObject = spinSystemLink()
            
            linkObject.score = presentPeaks
            linkObject.simulatedPeaks = listWithSimulatedPeaks
            linkObject.realPeaks = listWithPresentPeaks
            linkObject.notFoundSimulatedPeaks += listWithNotFoundSimulatedPeaks
            linkObject.spinSystem1 = spinSys1
            linkObject.spinSystem2 = spinSys2
            
            resA.linkDict[(spinSys1.spinSystemNumber*hc+spinSys2.spinSystemNumber)] = linkObject  
              
     
  cdef list commonElementInLists(self, list listOfLists):
    
    cdef int length
    cdef list intersect
    
    length= len(listOfLists)
     
    if length > 0 :
      
      intersect = list(set(listOfLists[0]).intersection(*listOfLists))
            
    return intersect
             

  cdef void calculateAllPeakContributions(self):
        
    
    cdef aSpectrum spectrum
    
    cdef object ccpnSpectrum
    
    cdef str refExpName
    
    cdef str info

    for spectrum in self.DataModel.spectra:                 # Determine for each dimension of every peak in all (used) spectra, which resonances can contribute to the peak
      
      spectrum.setupPeaks()   
              
      info = 'Evaluating dimensional contributions to peaks: ' + spectrum.name
      
      self.updateInfoText(info)
      
      self.findPossiblePeakContributions(spectrum)
    

  cdef void findPossiblePeakContributions(self,aSpectrum  spectrum) :
    
  
    
    cdef list resonances
    cdef myResonance resonance
    cdef aPeak peak
    cdef aDimension dim
    cdef dict entryForSpectrum
    cdef list listWithPeaks
    cdef dict newlib
    cdef aPeak firstPeak
    cdef object atomSite
    cdef str atomSiteName
    cdef str isotope
    cdef dict dimAtomsDict
    cdef double ppmValue
    cdef object possible
    
    cdef object onlyIntra
    
    cdef int numberOfDimensions
    
    cdef int serial
    
    
    dimAtomsDict = {}

    
    if spectrum.peaks:
      
      firstPeak = spectrum.peaks[0]
      
      for dim in firstPeak.dimensions :
        
        atomSite = dim.ccpnDim.dataDim.expDim.refExpDim.findFirstRefExpDimRef().expMeasurement.findFirstAtomSite()
        
        atomSiteName = atomSite.name
        
        isotope = atomSite.isotopeCode
        
        if atomSiteName == 'H' :
          
          dimAtomsDict[dim.dimNumber] = self.DataModel.Hresonances
        
        elif atomSiteName == 'N' :
          
          dimAtomsDict[dim.dimNumber] = self.DataModel.Nresonances
          
        elif atomSiteName == 'CA' :
          
          dimAtomsDict[dim.dimNumber] = self.DataModel.CAresonances
          
        elif atomSiteName == 'CB' :
          
          dimAtomsDict[dim.dimNumber] = self.DataModel.CBresonances
          
        elif atomSiteName == 'CO' :
          
          dimAtomsDict[dim.dimNumber] = self.DataModel.COresonances
          
        elif isotope == '1H' :
          
          dimAtomsDict[dim.dimNumber] = self.DataModel.HXresonances
          
        elif isotope == '15N' :
          
          dimAtomsDict[dim.dimNumber] = self.DataModel.NXresonances
          
        elif isotope == '13C' :
          
          dimAtomsDict[dim.dimNumber] = self.DataModel.CXresonances
          
        else :
          
          dimAtomsDict[dim.dimNumber] = 0
        
          print 'Isotope not implemented yet'  
    
    
    
    for peak in spectrum.peaks :
      
      onlyIntra = True
      
      ssDict = {}
      
      for dim in peak.dimensions :                                      # Here I quickly check whether the peak does not have a completely intra-residual assignment (for instance in PDSD's this might be the case), and the peak should therefor be ignored, since it does not contain any sequential information.
        
        for contrib in dim.ccpnDim.peakDimContribs :
          
          if contrib.resonance and contrib.resonance.resonanceGroup :
          
            if contrib.resonance.resonanceGroup.serial in ssDict :
            
              ssDict[contrib.resonance.resonanceGroup.serial] +=1
              
            else :  
            
              ssDict[contrib.resonance.resonanceGroup.serial] = 1
            
            
      if ssDict :
        
        for serial, numberOfDimensions in ssDict.items() :
          
          if numberOfDimensions != len(peak.dimensions) :               # A peak is considered a pure auto-peak only if all dimensions are assigned to the same spin system. 
            
            onlyIntra = False
            
      else :
        
        onlyIntra = False
            
            
      if not onlyIntra :                                                # Only proceed when the peak is not an auto-peak.   
                  
        for dim in peak.dimensions :
          
          resonances = dimAtomsDict[dim.dimNumber]
          
          ppmValue = dim.ppmValue
          
          for resonance in resonances :
          
            possible = False
            
            if self.useDimenionalAssignments and len(dim.ccpnDim.peakDimContribs) > 0 :                                         # Here I say that if there are already resonances assigned to this peak dimension, these should be the only resonances taken into account.
              
              
              
              for contrib in dim.ccpnDim.peakDimContribs :
                
                if contrib.resonance == resonance.ccpnResonance :
                  
                  possible = True
            
            
          
            elif abs(resonance.CS - ppmValue) <= getAnalysisDataDim(dim.ccpnDim.dataDim).assignTolerance :                    # If there are no peak resonances assigned to this peak dimension, every resonance within the tolerance qualifies for being a contribution
              
              possible = True
              
            if possible :
                              
              dim.possibleContributions.append(resonance)
              
              if spectrum.name in resonance.peakDimsLib :
                
                entryForSpectrum = resonance.peakDimsLib[spectrum.name]
                
                if dim.dimNumber in entryForSpectrum :
                  
                  listWithPeaks = entryForSpectrum[dim.dimNumber]
                  listWithPeaks.append(peak)
                
                else :
                  
                  entryForSpectrum[dim.dimNumber] = [peak]
                  
              else : 
                
                newlib = {}
                newlib[dim.dimNumber] = [peak]
                
                resonance.peakDimsLib[spectrum.name] = newlib
        

  cdef void setupSpinSystemExchange(self):
    
    cdef myDataModel DataModel
    
    cdef list residues
    cdef dict residuesByCcpCode

    cdef dict justTypedSpinSystems
    cdef dict allSpinSystems
    cdef dict tentativeSpinSystems
    cdef dict allSpinSystemsWithoutAssigned
    
    cdef object useAssignments
    cdef object useTentative
    
    cdef dict dictio
    
    cdef str key
    
    cdef list listWithSpinSystems
    
    cdef list spinSysList
    
    cdef mySpinSystem spinSys
    
    cdef mySpinSystem spinSys2
    
    cdef str ccpCode
    
    cdef int seqCode
    
    cdef aResidue res
    

    
    
      
      
    useAssignments = self.useAssignments
    useTentative = self.useTentative
    
    DataModel = self.DataModel
    
    residues = DataModel.myChain.residues

    residuesByCcpCode = self.makePrivateCopyOfDictContainingLists(DataModel.myChain.residuesByCcpCode)

    previouslyAssignedSpinSystems = DataModel.previouslyAssignedSpinSystems
    justTypedSpinSystems =DataModel.justTypedSpinSystems
    allSpinSystems = DataModel.mySpinSystems
    tentativeSpinSystems = DataModel.tentativeSpinSystems
    allSpinSystemsWithoutAssigned = DataModel.allSpinSystemsWithoutAssigned

    
    for key, spinSysList in allSpinSystems.items():                    # Cleaning up
      
      for spinSys in spinSysList :
        
        spinSys.allowedResidues = set()
        
        
    
    
    if useAssignments :                                                 # If the already made assignments are used, the dictionary will just include unassigned spinsystems and Joker spinsystems

      dictio = allSpinSystemsWithoutAssigned
      
      for key, spinSysList in previouslyAssignedSpinSystems.items() :
        
        for spinSys in spinSysList :
          
          spinSys.allowedResidues = set([spinSys.ccpnSeqCode])
          
          res = residues[spinSys.ccpnSeqCode-1]
          
          if spinSys.ccpCode in residuesByCcpCode :
            
            if res in residuesByCcpCode[spinSys.ccpCode] :
             
              residuesByCcpCode[spinSys.ccpCode].remove(res)                
    
    else :                                                                       # If the already assigned are chosen to still be permutated, the dictionary will also include the ready assinged spinsystems
      dictio = allSpinSystems
      
      
    listWithSpinSystems = []  
      
    for key in dictio :
      
      listWithSpinSystems = listWithSpinSystems + (dictio[key])
      
    listWithSpinSystems = list(set(listWithSpinSystems))
    
    
    for spinSys in listWithSpinSystems :
    
      if spinSys.ccpCode :  
          
        if spinSys.ccpCode in dictio :
            
          for spinSys2 in dictio[spinSys.ccpCode]  :
            
            if spinSys != spinSys2 and not (spinSys.isJoker and spinSys2.isJoker) :
                
              spinSys.exchangeSpinSystems.append(spinSys2)
            
          
      elif spinSys.tentativeCcpCodes :
            
        for ccpCode in spinSys.tentativeCcpCodes :
          
          if ccpCode in dictio :

            for spinSys2 in dictio[ccpCode] :
                  
              if useTentative and spinSys2.tentativeSeqCodes:      #not spinSys2.ccpCode : 
                      
                if spinSys != spinSys2 and not (spinSys.isJoker and spinSys2.isJoker) and list(set(spinSys.tentativeSeqCodes).intersection(spinSys2.tentativeSeqCodes)) :
                    
                  spinSys.exchangeSpinSystems.append(spinSys2)
                
              elif spinSys != spinSys2 and not (spinSys.isJoker and spinSys2.isJoker) :
                  
                spinSys.exchangeSpinSystems.append(spinSys2)
      
      elif spinSys.aminoAcidProbs :
        
        for ccpCode in spinSys.aminoAcidProbs.keys():
          
          if ccpCode in dictio :

            for spinSys2 in dictio[ccpCode] :
              
              if spinSys != spinSys2 and not (spinSys.isJoker and spinSys2.isJoker) :
                
                spinSys.exchangeSpinSystems.append(spinSys2)
      
              
      spinSys.exchangeSpinSystems = list(set(spinSys.exchangeSpinSystems))          
      
          
      if useTentative and spinSys.tentativeSeqCodes :
          
        for seqCode in spinSys.tentativeSeqCodes :
            
          spinSys.allowedResidues.add(seqCode)         
          
      elif spinSys.ccpCode :
        
        if spinSys.ccpCode in residuesByCcpCode :
        
          for res in residuesByCcpCode[spinSys.ccpCode] :
           
            spinSys.allowedResidues.add(res.seqCode)     
        
      elif spinSys.tentativeCcpCodes :                        # Only end up in here when not useTentative
          
        for ccpCode in spinSys.tentativeCcpCodes :
          
          if ccpCode in residuesByCcpCode :
          
            for res in residuesByCcpCode[ccpCode] :
              
              spinSys.allowedResidues.add(res.seqCode)
              
      elif spinSys.aminoAcidProbs :
        
        for ccpCode in spinSys.aminoAcidProbs.keys() :
          
          if ccpCode in residuesByCcpCode :
          
            for res in residuesByCcpCode[ccpCode] :
              
              spinSys.allowedResidues.add(res.seqCode)
        
        
          

      


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


  cdef double getCoLabellingFractionNew(self,aSpectrum spectrum,list atoms):
    
    cdef object scheme
    cdef object ccpnSpectrum
    cdef aResidue x
    cdef list ccpnResidues
    
    cdef double colabelling
    
    
    scheme = spectrum.labellingScheme
    
    if scheme is True :
      
      colabelling = self.calculateCoLabellingFromExperimentNew(spectrum, atoms)
      
      return colabelling
    
#    # TODO: look at this and fix
#    elif scheme :
#      
#      residueDict = {} 
#      
#      for atomName, res in zip(atomNames,  residues) :                                                                  
#        
#        if res in residueDict : 
#          
#          residueDict[res].append(atomName)                                                             # You can seemingly use objects of user-defined types as dictionary keys in python :-)
#          
#        else :
#          
#          residueDict[res] = [atomName]
#           
#      colabelling = 1     
#           
#      for res, atomNameList in residueDict.items() :
#        
#        numberOfAtoms = len(atomNameList)
#        
#        if numberOfAtoms == 1 :
#          
#          colabelling = colabelling * self.getAtomLabellingFraction('protein', res.ccpCode, atomNameList[0], scheme)
#          
#        elif numberOfAtoms == 2 :
#          
#          colabelling = colabelling * self.getAtomPairLabellingFraction_intra('protein', res.ccpCode, atomNameList[0], atomNameList[1], scheme)
#          
#        else :
#          
#          print 'something went wrong'
#      
#      print colabelling    
#      return colabelling    
          
          
          
  cdef double getCoLabellingFraction(self,aSpectrum spectrum,list residues, list atomNames):
    
    cdef object scheme
    cdef object ccpnSpectrum
    cdef aResidue x
    cdef list ccpnResidues
    
    cdef double colabelling
    
    
    scheme = spectrum.labellingScheme
    
    if scheme is True :
      
      ccpnSpectrum = spectrum.ccpnSpectrum
      ccpnResidues = [x.ccpnResidue for x in residues]
      
      colabelling = self.calculateCoLabellingFromExperiment(ccpnSpectrum, ccpnResidues, atomNames)
      
      return colabelling
    
    elif scheme :
      
      residueDict = {} 
      
      for atomName, res in zip(atomNames,  residues) :                                                                  
        
        if res in residueDict : 
          
          residueDict[res].append(atomName)                                                             # You can seemingly use objects of user-defined types as dictionary keys in python :-)
          
        else :
          
          residueDict[res] = [atomName]
           
      colabelling = 1     
           
      for res, atomNameList in residueDict.items() :
        
        numberOfAtoms = len(atomNameList)
        
        if numberOfAtoms == 1 :
          
          colabelling = colabelling * self.getAtomLabellingFraction('protein', res.ccpCode, atomNameList[0], scheme)
          
        elif numberOfAtoms == 2 :
          
          colabelling = colabelling * self.getAtomPairLabellingFraction_intra('protein', res.ccpCode, atomNameList[0], atomNameList[1], scheme)
          
        else :
          
          print 'something went wrong'
      
      #print colabelling    
      return colabelling    
               
      
  cdef double calculateCoLabellingFromExperimentNew(self, aSpectrum spectrum, list atoms) :
    
    #print '1'
    
    cdef object mixture
    
    cdef double TotalCoLabellingFraction
    
    cdef dict residueDict
    
    cdef anAtom atom
    
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
    
    cdef list isotopes
    
    cdef list atomNameList
    
    cdef double addedColabelling
    
    #mixtures = spectrum.ccpnSpectrum.experiment.labeledMixtures
    
   # print '2'
  
    mixture = spectrum.ccpnSpectrum.experiment.findFirstLabeledMixture()
      
    molLabelFractions = mixture.molLabelFractions
    
    molWeightSum = sum([x.weight for x in molLabelFractions])                                                               # This is the sum of weight of all labelling patterns present in the sample
    
   # print '4'
    
    
    residueDict = {} 
    
    for atom in atoms :                                                             # Group Atoms by residue
    
#      print atom.atomName
      
      molResidue = atom.residue.ccpnResidue.molResidue
      
      ccpnAtom = atom.ccpnAtom
      
      if molResidue in residueDict :
        
        residueDict[molResidue].append(atom)
        
      else :
      
        residueDict[molResidue] = [atom]  
        
      
         
         
    TotalCoLabellingFraction = 0.0   
  
  #  print '8'  
         
    for mlf in molLabelFractions:                                                                                                                       # Now I want to loop over all different labelling patterns that are present in the sample (denoted with capital letter in the analysis GUI)
      
      
 #     print '9'
      
      molFactor = mlf.weight / molWeightSum
        
      molLabel = mlf.molLabel
      
      colabellingOfAtomsInThismolLabelFraction= 1.0
      
 #     print '10'
      
      for molResidue, atomList in residueDict.items() :                                                                                         # Loop over involved residues
        
 #       print '11'
        
        #molResidue = mixture.labeledMolecule.molecule.findFirstMolResidue(serial=resID)
        resLabel = molLabel.findFirstResLabel(resId=molResidue.serial)
        resLabelFractions = resLabel.resLabelFractions
        
 #       print '12'
        
        rlfWeightSum = sum([x.weight for x in resLabelFractions])
        
  #      print '13'
        
        colabellingOfAtomsWithinThisResidue = 0
        
        for rlf in resLabelFractions :                                                                                                                         # Loop over the different resLabelFractions, bascically looping over isotopomers of one residue in a specific labelling pattern
          
          resFactor = rlf.weight / rlfWeightSum                                                                                                           # The relative weight of the different isotopomers
          
 #         print '14'
          if len(atomList) == 1 :
            
            #atomName = atomNameList[0]
            
            #chemAtom = molResidue.chemCompVar.findFirstChemAtom(name=atomName)
            
            
            atom = atomList[0]
            #chemAtom = atom.ccpnAtom.chemAtom
            
       #     print '15'
            
            
            #subType = chemAtom.subType
            
            atomName = atom.atomName
            
            #print atomName
            
            #fracDict = getIsotopomerSingleAtomFractions(rlf.isotopomers, atomName, subType)
            
            fracDict = self.getIsotopomerSingleAtomFractionsForChemAtom(rlf.isotopomers,  atom)
            
            fraction = fracDict.get(self.guessSpinHalfIsotopeFromAtomName(atomName), 1.0)
            
       #     print '16'
            
            
          else :                                                    # TODO:UPDATE
            #print '17'
            subTypes = []
            isotopes = []
            
            for atom in atomList :
              
              atomName = atom.atomName
              chemAtom = atom.ccpnAtom.chemAtom
              
              #print 'atomName....'
              #print atomName
              
              subTypes.append(chemAtom.subType)
              isotopes.append(self.guessSpinHalfIsotopeFromAtomName(atomName))
              
              #print '19'

            atomNameList = [atom.atomName for atom in atomList]
            #print '19A'
            #print rlf.isotopomers
            #print type(rlf.isotopomers)
            #print type(atomNameList)
            #print type(subTypes)
            fracDict =self.getIsotopomerAtomSetFractions(rlf.isotopomers, atomNameList, subTypes)
            #print '19B'
            fraction = fracDict.get(tuple(isotopes), 1.0)
            #print '19C'
            
          #print '19D' 
          
          #print fraction
          #print resFactor
          addedColabelling =  fraction * resFactor
          #print '19D1'
          colabellingOfAtomsWithinThisResidue = colabellingOfAtomsWithinThisResidue +  addedColabelling #+ (fraction * resFactor)
          #print '20'  
      
        colabellingOfAtomsInThismolLabelFraction = colabellingOfAtomsInThismolLabelFraction * colabellingOfAtomsWithinThisResidue
        
        #print '21'
            
      TotalCoLabellingFraction = TotalCoLabellingFraction + ( colabellingOfAtomsInThismolLabelFraction * molFactor )
      
      #print '22'
    #print TotalCoLabellingFraction
    return TotalCoLabellingFraction 


  cdef double calculateCoLabellingFromExperiment(self, object ccpnSpectrum, list ccpnResidues, list atomNames) :
    
    #print '1'
    
    cdef object mixture
    
    cdef double TotalCoLabellingFraction
    
    cdef dict residueDict
    
    mixtures = ccpnSpectrum.experiment.labeledMixtures
    
    #print '2'

    mixture = ccpnSpectrum.experiment.findFirstLabeledMixture()
      
    molLabelFractions = mixture.molLabelFractions
    
    molWeightSum = sum([x.weight for x in molLabelFractions])                                                               # This is the sum of weight of all labelling patterns present in the sample
    
    #print '4'
    residueDict = {} 
    
    for atomName, res in zip(atomNames,  ccpnResidues) :                                                                       # I want to group the atoms by residue
      
      #print '5'
      
      if res.molResidue.serial in residueDict : 
        
        #print '6'
        
        residueDict[res.molResidue.serial].append(atomName)
        
      else :
        
        #print '7'
        
        residueDict[res.molResidue.serial] = [atomName]
         
         
         
    TotalCoLabellingFraction = 0   
  
    #print '8'  
         
    for mlf in molLabelFractions:                                                                                                                       # Now I want to loop over all different labelling patterns that are present in the sample (denoted with capital letter in the analysis GUI)
      
      
      #print '9'
      
      molFactor = mlf.weight / molWeightSum
        
      molLabel = mlf.molLabel
      
      colabellingOfAtomsInThismolLabelFraction= 1
      
      #print '10'
      
      for resID, atomNameList in residueDict.items() :                                                                                         # Loop over involved residues
        
        #print '11'
        
        molResidue = mixture.labeledMolecule.molecule.findFirstMolResidue(serial=resID)
        resLabel = molLabel.findFirstResLabel(resId=resID)
        resLabelFractions = resLabel.resLabelFractions
        
        #print '12'
        
        rlfWeightSum = sum([x.weight for x in resLabelFractions])
        
        #print '13'
        
        colabellingOfAtomsWithinThisResidue = 0
        
        for rlf in resLabelFractions :                                                                                                                         # Loop over the different resLabelFractions, bascically looping over isotopomers of one residue in a specific labelling pattern
          
          resFactor = rlf.weight / rlfWeightSum                                                                                                           # The relative weight of the different isotopomers
          
          #print '14'
          if len(atomNameList) == 1 :
            
            atomName = atomNameList[0]
            
            chemAtom = molResidue.chemCompVar.findFirstChemAtom(name=atomName)
            
            #print '15'
            
            if not chemAtom :                                                                                                                   # This happens only when you try to calculate colabelling involving an atom that does not exist, like the backbone amide proton in proline
              
              return 0
            
            
            subType = chemAtom.subType
            
            fracDict = getIsotopomerSingleAtomFractions(rlf.isotopomers, atomName, subType)
            
            fraction = fracDict.get(self.guessSpinHalfIsotopeFromAtomName(atomName), 1.0)
            
            #print '16'
            
            
          else :
            print '17'
            subTypes = []
            isotopes = []
            
            for atomName in atomNameList :
              
              
              chemAtom = molResidue.chemCompVar.findFirstChemAtom(name=atomName)
              
              print '18'
              
              if not chemAtom :
                
                return 0
              
              
              subTypes.append(chemAtom.subType)
              isotopes.append(self.guessSpinHalfIsotopeFromAtomName(atomName))
              
              print '19'

            print '19A'
            print rlf.isotopomers
            print type(rlf.isotopomers)
            print type(atomNameList)
            print type(subTypes)
            fracDict =self.getIsotopomerAtomSetFractions(rlf.isotopomers, atomNameList, subTypes)
            print '19B'
            fraction = fracDict.get(tuple(isotopes), 1.0)
            print '19C'
            
       #   print '19D' 
       #   print fraction
      #    print resFactor
          addedColabelling =  fraction * resFactor
       #   print '19D1'
          colabellingOfAtomsWithinThisResidue = colabellingOfAtomsWithinThisResidue +  addedColabelling #+ (fraction * resFactor)
        #  print '20'  
      
        colabellingOfAtomsInThismolLabelFraction = colabellingOfAtomsInThismolLabelFraction * colabellingOfAtomsWithinThisResidue
        
        #print '21'
            
      TotalCoLabellingFraction = TotalCoLabellingFraction + ( colabellingOfAtomsInThismolLabelFraction * molFactor )
      
      #print '22'
 
    return TotalCoLabellingFraction 
            
            
  cdef str guessSpinHalfIsotopeFromAtomName(self, str atomName) :

    
    if atomName[0] == 'H':
      isotope = '1H'
    elif atomName[0] == 'C':
      isotope = '13C'
    elif atomName[0] == 'N':
      isotope = '15N'
    elif atomName[0] == 'P':
      isotope = '31P'
      
    return isotope
          

  cdef dict getIsotopomerAtomSetFractions(self,set isotopomers, list atomNames, list subTypes):
    """Descrn: Get the combined isotope proportions for a given set of named
               atoms within a set of isotopomers. Fractions normalised to 1.0
       Inputs: List of ChemCompLabel.Isotopomers, 2-Tuple of Words (ChemAtom.name), 2-Tuple of Words (ChemAtom.subType)
       Output: Dict of Tuple:Float - (IsotopeCode,IsotopeCode):fraction
    """
    
    #print '30'

    fractionDict = {}

    isoWeightSum = sum([x.weight for x in isotopomers])
  
    #atLabels   = []
    #sumWeights = []

#    atLabels   = [None, None]
#    sumWeights = [None, None]
   
    
    for isotopomer in isotopomers:

      i = 0
      atLabels   = []
      sumWeights = []
      for atomName in atomNames:
        atLabels.append(isotopomer.findAllAtomLabels(name=atomName, 
                                                   subType=subTypes[i]))
        sumWeights.append(sum([x.weight for x in atLabels[i]]))


#        atLabels[i] = isotopomer.findAllAtomLabels(name=atomNames[i], 
#                                                 subType=subTypes[i])
#        sumWeights[i] = sum([x.weight for x in atLabels[i]])



        i = i + 1
    
      # Done this way to guard against the divisor becoming zero
      factor  = isotopomer.weight / isoWeightSum
    
      divisor = 1

      for sumWeight in sumWeights :
        divisor = divisor * sumWeight


      if len(atLabels) == 2 :
    
        for atl0 in atLabels[0]:
          for atl1 in atLabels[1]:
      
            if atl0 is atl1:
              contrib = atl0.weight * factor / sumWeights[0]  
            else:
              contrib = atl0.weight * atl1.weight * factor / divisor
          
            key = (atl0.isotopeCode, atl1.isotopeCode)
            fractionDict[key] = fractionDict.get(key, 0.0) + contrib



      if len(atLabels) == 3 :
    
        for atl0 in atLabels[0]:
          for atl1 in atLabels[1]:
            for atl2 in atLabels[2]:
      
              if atl0 == atl1 and atl0 == atl2:
                contrib = atl0.weight * factor / sumWeights[0]

              elif atl0 is atl1 :
                contrib = atl0.weight * factor / sumWeights[0] * atl2.weight * factor / divisor

              elif atl0 is atl2 :
                contrib = atl0.weight * factor / sumWeights[0] * atl1.weight * factor / divisor

              elif atl1 is atl2 :
                contrib = atl1.weight * factor / sumWeights[0] * atl0.weight * factor / divisor

              else:
                contrib = atl0.weight * atl1.weight * atl2.weight * factor / divisor


              key = (atl0.isotopeCode, atl1.isotopeCode,atl2.isotopeCode)
              fractionDict[key] = fractionDict.get(key, 0.0) + contrib



      if len(atLabels) == 4 :
    
        for atl0 in atLabels[0]:
          for atl1 in atLabels[1]:
            for atl2 in atLabels[2]:
              for atl3 in atLabels[3]:  

                contrib = atl0.weight * atl1.weight * atl2.weight * atl3.weight * factor / divisor


                key = (atl0.isotopeCode, atl1.isotopeCode,atl2.isotopeCode,atl3.isotopeCode)
                fractionDict[key] = fractionDict.get(key, 0.0) + contrib
              





    return fractionDict
    

  cdef dict getIsotopomerSingleAtomFractionsForChemAtom(self,  set isotopomers, anAtom atom) :
    
    cdef dict fractionDict
    cdef double isoWeightSum
    
    cdef object x
    
    cdef object isotopomer
    
    cdef set atomLabels
    cdef frozenset allAtomLabels
    cdef object atomLabel
    cdef double atWeightSum
    
    cdef double atFactor
    
    cdef str isotopeCode
    
    cdef double contrib
    
    
    
    
    fractionDict = {}
    isoWeightSum =  sum([x.weight for x in isotopomers])

    for isotopomer in isotopomers :
      
      atomLabels = atom.labelInfoTemp[isotopomer]
  
      atWeightSum = sum([x.weight for x in atomLabels])
      atFactor    = isotopomer.weight  / isoWeightSum
      
      for atomLabel in atomLabels:
        
        isotopeCode = atomLabel.isotopeCode

        contrib     = atFactor * atomLabel.weight / atWeightSum

        fractionDict[isotopeCode] = fractionDict.get(isotopeCode,0.0) + contrib
  
    return fractionDict
      

  cdef void annealingSub(self, double AcceptanceConstant,int amountOfStepsPerTemperature,list listWithSpinSystems):
     
    cdef int improvements
    cdef int equals
    cdef int worse
    
    cdef mySpinSystem A
    cdef mySpinSystem B
    
    cdef aResidue currentResidueA
    cdef aResidue currentResidueB
    
    cdef aResidue previousResA
    cdef aResidue previousResB
    
    cdef aResidue nextResA
    cdef aResidue nextResB
    
    cdef list exchangeSpinSystems
    
  
    cdef list oldPeaks
    cdef list newPeaks
    cdef set PeakSet
    
    cdef int oldLinkScore
    cdef int newLinkScore
    
    cdef mySpinSystem spinSystemAm1
    cdef mySpinSystem spinSystemAp1 
    
    cdef mySpinSystem spinSystemBm1
    cdef mySpinSystem spinSystemBp1 
    
    cdef int spinSystemNumberA
    cdef int spinSystemNumberAm1
    cdef int spinSystemNumberAp1
    cdef int spinSystemNumberB
    cdef int spinSystemNumberBm1
    cdef int spinSystemNumberBp1
    
    cdef spinSystemLink linkAB
    cdef spinSystemLink linkBA
    cdef spinSystemLink linkAAm1
    cdef spinSystemLink linkBBm1
    cdef spinSystemLink linkAAp1
    cdef spinSystemLink linkBBp1
    cdef spinSystemLink linkABm1
    cdef spinSystemLink linkABp1
    cdef spinSystemLink linkBAm1
    cdef spinSystemLink linkBAp1
    
    cdef double cut
    
    cdef aPeak peak
    
    cdef set peakSet
    
    cdef int hc
    
    cdef double score
    
    cdef double DeltaScore
    
    #cdef dict linkDict
    
    hc = self.hc
    
    score = self.score
    
    improvements = 0
    equals = 0
    worse = 0
    
    #exp = math.exp
    cutoff = random_sample #random.random
    
    
    #sample = self.Csample #random.choice
    #randint = random.randint
    
    #calcDeltaPeakScore = self.CcalcDeltaPeakScore
    
    
  
    for aTry in xrange(amountOfStepsPerTemperature) :  
      
  
      #print 'been here'
      
      #print listWithSpinSystems
      #acceptThisOne = False
      
      A = listWithSpinSystems[randint(0, len(listWithSpinSystems))]
      
      exchangeSpinSystems = A.exchangeSpinSystems
      
      if exchangeSpinSystems :
        
        B = exchangeSpinSystems[randint(0, len(exchangeSpinSystems))]
        
      else :
      
        continue  
        
      #A = self.Csample(listWithSpinSystems)
      
      #print A.exchangeSpinSystems

     
      #B = self.Csample(A.exchangeSpinSystems)
      
  
      
  
      
      currentResidueA = A.currentResidueAssignment
      currentResidueB = B.currentResidueAssignment
      
      
      if currentResidueA and currentResidueB:             # Both spin systems are assigned to a residue
        
  
        
        if currentResidueA.seqCode in B.allowedResidues and currentResidueB.seqCode in A.allowedResidues :          # The switch is possible because aa types of residues fits with possible aa types of spinsystems
  
                     
          previousResA = currentResidueA.previousResidue
          previousResB = currentResidueB.previousResidue
          
          nextResA = currentResidueA.nextResidue
          nextResB = currentResidueB.nextResidue
          
          oldPeaks = []
          newPeaks = []
          
          oldLinkScore = 0
          newLinkScore = 0
         
          if previousResB and currentResidueA.seqCode == previousResB.seqCode :                                                                            # A and B happen to be a sequential pair AB
           
            AisJoker = A.isJoker
            BisJoker = B.isJoker
        
            if not AisJoker and not BisJoker :
              
              spinSystemNumberA = A.spinSystemNumber
              spinSystemNumberB = B.spinSystemNumber
              
              linkDict = currentResidueA.linkDict
            
              linkAB = linkDict[(spinSystemNumberA*hc+spinSystemNumberB)] 
              oldPeaks += (linkAB.realPeaks)
              oldLinkScore += linkAB.score
              
              linkBA = linkDict[(B.spinSystemNumber*hc+A.spinSystemNumber)] 
              newPeaks += (linkBA.realPeaks)
              newLinkScore += linkBA.score
            
              if previousResA: 
              
                spinSystemAm1 = previousResA.currentSpinSystemAssigned
                
                if not spinSystemAm1.isJoker :
  
                  spinSystemNumberAm1 = spinSystemAm1.spinSystemNumber
                  
                  linkDict = previousResA.linkDict
              
                  linkAAm1 = linkDict[(spinSystemNumberAm1*hc+spinSystemNumberA)]
                  oldPeaks += (linkAAm1.realPeaks)
                  oldLinkScore += linkAAm1.score
                  
                  linkBAm1 = linkDict[(spinSystemNumberAm1*hc+spinSystemNumberB)] 
                  newPeaks += (linkBAm1.realPeaks)
                  newLinkScore += linkBAm1.score
                  
              if  nextResB :
              
                spinSystemBp1 = nextResB.currentSpinSystemAssigned 
                
                if not spinSystemBp1.isJoker :
                  
                  spinSystemNumberBp1 = spinSystemBp1.spinSystemNumber
                  
                  linkDict = currentResidueB.linkDict
              
                  linkBBp1 = linkDict[(spinSystemNumberB*hc+spinSystemNumberBp1)] 
                  oldPeaks += (linkBBp1.realPeaks)
                  oldLinkScore += linkBBp1.score
  
                  linkABp1 = linkDict[(spinSystemNumberA*hc+spinSystemNumberBp1)] 
                  newPeaks += (linkABp1.realPeaks)
                  newLinkScore += linkABp1.score
              
  
  
            
            elif AisJoker :
              
              spinSystemNumberB = B.spinSystemNumber
              
              if previousResA : 
              
                spinSystemAm1 = previousResA.currentSpinSystemAssigned
                
                if not spinSystemAm1.isJoker :
  
                  spinSystemNumberAm1 = spinSystemAm1.spinSystemNumber
                  
                  linkBAm1 = previousResA.linkDict[(spinSystemNumberAm1*hc+spinSystemNumberB)] 
                  newPeaks += (linkBAm1.realPeaks)
                  newLinkScore += linkBAm1.score
                  
              if  nextResB :
              
                spinSystemBp1 = nextResB.currentSpinSystemAssigned 
                
                if not spinSystemBp1.isJoker :
                  
                  spinSystemNumberBp1 = spinSystemBp1.spinSystemNumber
              
                  linkBBp1 = currentResidueB.linkDict[(spinSystemNumberB*hc+spinSystemNumberBp1)] 
                  oldPeaks += (linkBBp1.realPeaks)
                  oldLinkScore += linkBBp1.score
                  
            elif B.isJoker :
              
              spinSystemNumberA = A.spinSystemNumber
            
              if previousResA : 
              
                spinSystemAm1 = previousResA.currentSpinSystemAssigned
                
                if not spinSystemAm1.isJoker :
  
                  spinSystemNumberAm1 = spinSystemAm1.spinSystemNumber
              
                  linkAAm1 = previousResA.linkDict[(spinSystemNumberAm1*hc+spinSystemNumberA)]
                  oldPeaks += (linkAAm1.realPeaks)
                  oldLinkScore += linkAAm1.score
                  
                  
              if  nextResB :
              
                spinSystemBp1 = nextResB.currentSpinSystemAssigned 
                
                if not spinSystemBp1.isJoker :
                  
                  spinSystemNumberBp1 = spinSystemBp1.spinSystemNumber
  
                  linkABp1 = currentResidueB.linkDict[(spinSystemNumberA*hc+spinSystemNumberBp1)] 
                  newPeaks += (linkABp1.realPeaks)
                  newLinkScore += linkABp1.score
              
  
  
  
              
              
              
              
              
              
              
          elif previousResA and currentResidueB.seqCode == previousResA.seqCode :                                                                                    # sequential pair BA
            
            AisJoker = A.isJoker
            BisJoker = B.isJoker
        
            if not AisJoker and not BisJoker :
              
              spinSystemNumberA = A.spinSystemNumber
              spinSystemNumberB = B.spinSystemNumber
              
              linkDict = currentResidueB.linkDict
              
              linkBA = linkDict[(spinSystemNumberB*hc+spinSystemNumberA)] 
              oldPeaks += (linkBA.realPeaks)
              oldLinkScore += linkBA.score
              
              linkAB = linkDict[(spinSystemNumberA*hc+spinSystemNumberB)] 
              newPeaks += (linkAB.realPeaks)
              newLinkScore += linkAB.score
              
              
              if previousResB : 
              
                spinSystemBm1 = previousResB.currentSpinSystemAssigned
                
                if not spinSystemBm1.isJoker :
                  
                  spinSystemNumberBm1 = spinSystemBm1.spinSystemNumber
                  
                  linkDict = previousResB.linkDict
                 
                  linkBBm1 = linkDict[(spinSystemNumberBm1*hc+spinSystemNumberB)]
                  oldPeaks += (linkBBm1.realPeaks)
                  oldLinkScore += linkBBm1.score
                  
                  linkABm1 = linkDict[(spinSystemNumberBm1*hc+spinSystemNumberA)] 
                  newPeaks += (linkABm1.realPeaks)
                  newLinkScore += linkABm1.score
                  
              if  nextResA :
              
                spinSystemAp1 = nextResA.currentSpinSystemAssigned 
                
                if not spinSystemAp1.isJoker :
                  
                  spinSystemNumberAp1 = spinSystemAp1.spinSystemNumber
                  
                  linkDict = currentResidueA.linkDict
           
                  linkAAp1 = linkDict[(spinSystemNumberA*hc+spinSystemNumberAp1)] 
                  oldPeaks += (linkAAp1.realPeaks)
                  oldLinkScore += linkAAp1.score
  
            
                  linkBAp1 = linkDict[(spinSystemNumberB*hc+spinSystemNumberAp1)] 
                  newPeaks += (linkBAp1.realPeaks)
                  newLinkScore += linkBAp1.score
                  
                  
            elif AisJoker :
  
              spinSystemNumberB = B.spinSystemNumber
  
              if previousResB : 
              
                spinSystemBm1 = previousResB.currentSpinSystemAssigned
                
                if not spinSystemBm1.isJoker :
                  
                  spinSystemNumberBm1 = spinSystemBm1.spinSystemNumber
                 
                  linkBBm1 = previousResB.linkDict[(spinSystemNumberBm1*hc+spinSystemNumberB)]
                  oldPeaks += (linkBBm1.realPeaks)
                  oldLinkScore += linkBBm1.score
                  
                  
              if  nextResA :
              
                spinSystemAp1 = nextResA.currentSpinSystemAssigned 
                
                if not spinSystemAp1.isJoker :
                  
                  spinSystemNumberAp1 = spinSystemAp1.spinSystemNumber
  
            
                  linkBAp1 = currentResidueA.linkDict[(spinSystemNumberB*hc+spinSystemNumberAp1)] 
                  newPeaks += (linkBAp1.realPeaks)
                  newLinkScore += linkBAp1.score
                  
                  
            elif BisJoker :
              
              spinSystemNumberA = A.spinSystemNumber
  
              if previousResB : 
              
                spinSystemBm1 = previousResB.currentSpinSystemAssigned
                
                if not spinSystemBm1.isJoker :
                  
                  spinSystemNumberBm1 = spinSystemBm1.spinSystemNumber
                  
                  linkABm1 = previousResB.linkDict[(spinSystemNumberBm1*hc+spinSystemNumberA)] 
                  newPeaks += (linkABm1.realPeaks)
                  newLinkScore += linkABm1.score
                  
              if  nextResA :
              
                spinSystemAp1 = nextResA.currentSpinSystemAssigned 
                
                if not spinSystemAp1.isJoker :
                  
                  spinSystemNumberAp1 = spinSystemAp1.spinSystemNumber
           
                  linkAAp1 = currentResidueA.linkDict[(spinSystemNumberA*hc+spinSystemNumberAp1)] 
                  oldPeaks += (linkAAp1.realPeaks)
                  oldLinkScore += linkAAp1.score
  
                  
              
          else :                                                                                                            # A and B are not sequential
          
            AisJoker = A.isJoker
            BisJoker = B.isJoker
          
            if not AisJoker and not BisJoker :
              
              spinSystemNumberA = A.spinSystemNumber
              spinSystemNumberB = B.spinSystemNumber
              
              if previousResA : 
              
                spinSystemAm1 = previousResA.currentSpinSystemAssigned
                
                if not spinSystemAm1.isJoker :
  
                  spinSystemNumberAm1 = spinSystemAm1.spinSystemNumber
                  
                  linkDict = previousResA.linkDict
              
                  linkAAm1 = linkDict[(spinSystemNumberAm1*hc+spinSystemNumberA)]
                  oldPeaks += (linkAAm1.realPeaks)
                  oldLinkScore += linkAAm1.score
                  
                  linkBAm1 = linkDict[(spinSystemNumberAm1*hc+spinSystemNumberB)] 
                  newPeaks += (linkBAm1.realPeaks)
                  newLinkScore += linkBAm1.score
                  
              if  nextResA :
              
                spinSystemAp1 = nextResA.currentSpinSystemAssigned 
                
                if not spinSystemAp1.isJoker :
                  
                  spinSystemNumberAp1 = spinSystemAp1.spinSystemNumber
                  
                  linkDict = currentResidueA.linkDict
           
                  linkAAp1 = linkDict[(spinSystemNumberA*hc+spinSystemNumberAp1)] 
                  oldPeaks += (linkAAp1.realPeaks)
                  oldLinkScore += linkAAp1.score
  
            
                  linkBAp1 = linkDict[(spinSystemNumberB*hc+spinSystemNumberAp1)] 
                  newPeaks += (linkBAp1.realPeaks)
                  newLinkScore += linkBAp1.score
                  
                  
              if previousResB : 
              
                spinSystemBm1 = previousResB.currentSpinSystemAssigned
                
                if not spinSystemBm1.isJoker :
                  
                  spinSystemNumberBm1 = spinSystemBm1.spinSystemNumber
                  
                  linkDict = previousResB.linkDict
                 
                  linkBBm1 = linkDict[(spinSystemNumberBm1*hc+spinSystemNumberB)]
                  oldPeaks += (linkBBm1.realPeaks)
                  oldLinkScore += linkBBm1.score
                  
                  linkABm1 = linkDict[(spinSystemNumberBm1*hc+spinSystemNumberA)] 
                  newPeaks += (linkABm1.realPeaks)
                  newLinkScore += linkABm1.score
                  
                  
              if  nextResB :
              
                spinSystemBp1 = nextResB.currentSpinSystemAssigned 
                
                if not spinSystemBp1.isJoker :
                  
                  spinSystemNumberBp1 = spinSystemBp1.spinSystemNumber
                  
                  linkDict = currentResidueB.linkDict
              
                  linkBBp1 = linkDict[(spinSystemNumberB*hc+spinSystemNumberBp1)] 
                  oldPeaks += (linkBBp1.realPeaks)
                  oldLinkScore += linkBBp1.score
  
                  linkABp1 = linkDict[(spinSystemNumberA*hc+spinSystemNumberBp1)] 
                  newPeaks += (linkABp1.realPeaks)
                  newLinkScore += linkABp1.score
                  
                  
  
                  
  
            
            elif AisJoker :
              
              spinSystemNumberB = B.spinSystemNumber
              
              if previousResA : 
              
                spinSystemAm1 = previousResA.currentSpinSystemAssigned
                
                if not spinSystemAm1.isJoker :
  
                  spinSystemNumberAm1 = spinSystemAm1.spinSystemNumber
                  
                  linkBAm1 = previousResA.linkDict[(spinSystemNumberAm1*hc+spinSystemNumberB)] 
                  newPeaks += (linkBAm1.realPeaks)
                  newLinkScore += linkBAm1.score
                  
              if  nextResA :
              
                spinSystemAp1 = nextResA.currentSpinSystemAssigned 
                
                if not spinSystemAp1.isJoker :
                  
                  spinSystemNumberAp1 = spinSystemAp1.spinSystemNumber
  
                  linkBAp1 = currentResidueA.linkDict[(spinSystemNumberB*hc+spinSystemNumberAp1)] 
                  newPeaks += (linkBAp1.realPeaks)
                  newLinkScore += linkBAp1.score
                  
                  
              if previousResB : 
              
                spinSystemBm1 = previousResB.currentSpinSystemAssigned
                
                if not spinSystemBm1.isJoker :
                  
                  spinSystemNumberBm1 = spinSystemBm1.spinSystemNumber
                 
                  linkBBm1 = previousResB.linkDict[(spinSystemNumberBm1*hc+spinSystemNumberB)]
                  oldPeaks += (linkBBm1.realPeaks)
                  oldLinkScore += linkBBm1.score
                  
                  
              if  nextResB :
              
                spinSystemBp1 = nextResB.currentSpinSystemAssigned 
                
                if not spinSystemBp1.isJoker :
                  
                  spinSystemNumberBp1 = spinSystemBp1.spinSystemNumber
              
                  linkBBp1 = currentResidueB.linkDict[(spinSystemNumberB*hc+spinSystemNumberBp1)] 
                  oldPeaks += (linkBBp1.realPeaks)
                  oldLinkScore += linkBBp1.score
                  
                  
                  
                  
                  
            elif BisJoker :
        
              spinSystemNumberA = A.spinSystemNumber
              
              if previousResA : 
              
                spinSystemAm1 = previousResA.currentSpinSystemAssigned
                
                if not spinSystemAm1.isJoker :
  
                  spinSystemNumberAm1 = spinSystemAm1.spinSystemNumber
              
                  linkAAm1 = previousResA.linkDict[(spinSystemNumberAm1*hc+spinSystemNumberA)]
                  oldPeaks += (linkAAm1.realPeaks)
                  oldLinkScore += linkAAm1.score
                  
              if  nextResA :
              
                spinSystemAp1 = nextResA.currentSpinSystemAssigned 
                
                if not spinSystemAp1.isJoker :
                  
                  spinSystemNumberAp1 = spinSystemAp1.spinSystemNumber
           
                  linkAAp1 = currentResidueA.linkDict[(spinSystemNumberA*hc+spinSystemNumberAp1)] 
                  oldPeaks += (linkAAp1.realPeaks)
                  oldLinkScore += linkAAp1.score
                  
                  
              if previousResB : 
              
                spinSystemBm1 = previousResB.currentSpinSystemAssigned
                
                if not spinSystemBm1.isJoker :
                  
                  spinSystemNumberBm1 = spinSystemBm1.spinSystemNumber
                  
                  linkABm1 = previousResB.linkDict[(spinSystemNumberBm1*hc+spinSystemNumberA)] 
                  newPeaks += (linkABm1.realPeaks)
                  newLinkScore += linkABm1.score
                  
                  
              if  nextResB :
              
                spinSystemBp1 = nextResB.currentSpinSystemAssigned 
                
                if not spinSystemBp1.isJoker :
                  
                  spinSystemNumberBp1 = spinSystemBp1.spinSystemNumber
  
                  linkABp1 = currentResidueB.linkDict[(spinSystemNumberA*hc+spinSystemNumberBp1)] 
                  newPeaks += (linkABp1.realPeaks)
                  newLinkScore += linkABp1.score
                  
                  
  
            
            
          
         
  
            
  
  
          peakSet = set(oldPeaks+newPeaks)
          DeltaScore = self.CcalcDeltaPeakScore(peakSet,oldPeaks,newPeaks) + newLinkScore - oldLinkScore
          
          
  
            
            
        
  
  
  
  
          if  DeltaScore == 0 :
          
            equals += 1
            
            B.currentResidueAssignment = currentResidueA
            A.currentResidueAssignment = currentResidueB
            
            currentResidueA.currentSpinSystemAssigned = B
            currentResidueB.currentSpinSystemAssigned = A
            
            score += DeltaScore
            
            
            
  
            
            for peak in peakSet :
              
              peak.degeneracy = peak.degeneracyTemp
              
  
  
  
          elif DeltaScore > 0 :
            
            improvements += 1
            
            B.currentResidueAssignment = currentResidueA
            A.currentResidueAssignment = currentResidueB
            
            currentResidueA.currentSpinSystemAssigned = B
            currentResidueB.currentSpinSystemAssigned = A
            
            score += DeltaScore
            
  
            
            for peak in peakSet :
              
              peak.degeneracy = peak.degeneracyTemp
              
  
            
  
          elif exp(AcceptanceConstant*DeltaScore) > cutoff() :     
            
            worse += 1
  
            B.currentResidueAssignment = currentResidueA
            A.currentResidueAssignment = currentResidueB
            
            currentResidueA.currentSpinSystemAssigned = B
            currentResidueB.currentSpinSystemAssigned = A
            
            score += DeltaScore
            
  
            
            for peak in peakSet :
              
              peak.degeneracy = peak.degeneracyTemp
              
  
            
            
  
          
          
          
  
      
      elif currentResidueA :                                                                      # spin system B is not assigned to any residue
        
        if currentResidueA.seqCode in B.allowedResidues :
          
          previousResA = currentResidueA.previousResidue
          nextResA = currentResidueA.nextResidue
  
          oldPeaks = []
          newPeaks = []
          
          oldLinkScore = 0
          newLinkScore = 0
          
          AisJoker = A.isJoker
          BisJoker = B.isJoker
          
          
          if not AisJoker and not BisJoker :
            
            spinSystemNumberA = A.spinSystemNumber
            spinSystemNumberB = B.spinSystemNumber
            
            if previousResA : 
            
              spinSystemAm1 = previousResA.currentSpinSystemAssigned
              
              if not spinSystemAm1.isJoker :
  
                spinSystemNumberAm1 = spinSystemAm1.spinSystemNumber
                
                linkDict =  previousResA.linkDict
            
                linkAAm1 = linkDict[(spinSystemNumberAm1*hc+spinSystemNumberA)]
                oldPeaks += (linkAAm1.realPeaks)
                oldLinkScore += linkAAm1.score
                
                linkBAm1 = linkDict[(spinSystemNumberAm1*hc+spinSystemNumberB)] 
                newPeaks += (linkBAm1.realPeaks)
                newLinkScore += linkBAm1.score
                
            if  nextResA :
            
              spinSystemAp1 = nextResA.currentSpinSystemAssigned 
              
              if not spinSystemAp1.isJoker :
                
                spinSystemNumberAp1 = spinSystemAp1.spinSystemNumber
                
                linkDict = currentResidueA.linkDict
         
                linkAAp1 = linkDict[(spinSystemNumberA*hc+spinSystemNumberAp1)] 
                oldPeaks += (linkAAp1.realPeaks)
                oldLinkScore += linkAAp1.score
  
          
                linkBAp1 = linkDict[(spinSystemNumberB*hc+spinSystemNumberAp1)] 
                newPeaks += (linkBAp1.realPeaks)
                newLinkScore += linkBAp1.score
          
          elif AisJoker :
            
            spinSystemNumberB = B.spinSystemNumber
            
            if previousResA : 
            
              spinSystemAm1 = previousResA.currentSpinSystemAssigned
              
              if not spinSystemAm1.isJoker :
  
                spinSystemNumberAm1 = spinSystemAm1.spinSystemNumber
                
                linkBAm1 = previousResA.linkDict[(spinSystemNumberAm1*hc+spinSystemNumberB)] 
                newPeaks += (linkBAm1.realPeaks)
                newLinkScore += linkBAm1.score
                
            if  nextResA :
            
              spinSystemAp1 = nextResA.currentSpinSystemAssigned 
              
              if not spinSystemAp1.isJoker :
                
                spinSystemNumberAp1 = spinSystemAp1.spinSystemNumber
         
                linkBAp1 = currentResidueA.linkDict[(spinSystemNumberB*hc+spinSystemNumberAp1)] 
                newPeaks += (linkBAp1.realPeaks)
                newLinkScore += linkBAp1.score
            
            
            
          elif BisJoker :
            
            spinSystemNumberA = A.spinSystemNumber
            
            if previousResA : 
            
              spinSystemAm1 = previousResA.currentSpinSystemAssigned
              
              if not spinSystemAm1.isJoker :
  
                spinSystemNumberAm1 = spinSystemAm1.spinSystemNumber
            
                linkAAm1 = previousResA.linkDict[(spinSystemNumberAm1*hc+spinSystemNumberA)]
                oldPeaks += (linkAAm1.realPeaks)
                oldLinkScore += linkAAm1.score
                
            if  nextResA :
            
              spinSystemAp1 = nextResA.currentSpinSystemAssigned 
              
              if not spinSystemAp1.isJoker :
                
                spinSystemNumberAp1 = spinSystemAp1.spinSystemNumber
         
                linkAAp1 = currentResidueA.linkDict[(spinSystemNumberA*hc+spinSystemNumberAp1)] 
                oldPeaks += (linkAAp1.realPeaks)
                oldLinkScore += linkAAp1.score
          
  
  
  
          
  
            
  
          peakSet = set(oldPeaks+newPeaks)
          DeltaScore = self.CcalcDeltaPeakScore(peakSet,oldPeaks,newPeaks) + newLinkScore - oldLinkScore
  
  
  
  
  
  
          if  DeltaScore == 0 :
            
            equals = equals + 1
            
            B.currentResidueAssignment = currentResidueA
            A.currentResidueAssignment = None
            
            currentResidueA.currentSpinSystemAssigned = B
            
            score += DeltaScore
            
  
            
            for peak in peakSet :
              
              peak.degeneracy = peak.degeneracyTemp
              
  
  
  
          elif DeltaScore > 0 :
            
            improvements = improvements + 1
            
            B.currentResidueAssignment = currentResidueA
            A.currentResidueAssignment = None
            
            currentResidueA.currentSpinSystemAssigned = B
            
            score += DeltaScore
            
  
            
            for peak in peakSet :
              
              peak.degeneracy = peak.degeneracyTemp
              
              
  
          elif exp(AcceptanceConstant*DeltaScore) > cutoff() :     
            
            worse = worse + 1
  
            B.currentResidueAssignment = currentResidueA
            A.currentResidueAssignment = None
            
            currentResidueA.currentSpinSystemAssigned = B
            
            score += DeltaScore
            
  
  
            
            for peak in peakSet :
              
              peak.degeneracy = peak.degeneracyTemp
              
            
            
  
          
  
      
              
      elif currentResidueB :                                                                      # spin system A is not assigned to any residue
        
        if currentResidueB.seqCode in A.allowedResidues :
  
          previousResB = currentResidueB.previousResidue
          nextResB = currentResidueB.nextResidue
  
          
          oldPeaks = []
          newPeaks = []
          
          oldLinkScore = 0
          newLinkScore = 0
          
          AisJoker = A.isJoker
          BisJoker = B.isJoker
          
          
          if not AisJoker and not BisJoker :
            
            spinSystemNumberA = A.spinSystemNumber
            spinSystemNumberB = B.spinSystemNumber
            
            if previousResB : 
            
              spinSystemBm1 = previousResB.currentSpinSystemAssigned
              
              if not spinSystemBm1.isJoker :
                
                spinSystemNumberBm1 = spinSystemBm1.spinSystemNumber
                
                linkDict = previousResB.linkDict
               
                linkBBm1 = linkDict[(spinSystemNumberBm1*hc+spinSystemNumberB)]
                oldPeaks += (linkBBm1.realPeaks)
                oldLinkScore += linkBBm1.score
                
                linkABm1 = linkDict[(spinSystemNumberBm1*hc+spinSystemNumberA)] 
                newPeaks += (linkABm1.realPeaks)
                newLinkScore += linkABm1.score
                
                
            if  nextResB :
            
              spinSystemBp1 = nextResB.currentSpinSystemAssigned 
              
              if not spinSystemBp1.isJoker :
                
                spinSystemNumberBp1 = spinSystemBp1.spinSystemNumber
                
                linkDict = currentResidueB.linkDict
            
                linkBBp1 = linkDict[(spinSystemNumberB*hc+spinSystemNumberBp1)] 
                oldPeaks += (linkBBp1.realPeaks)
                oldLinkScore += linkBBp1.score
  
                linkABp1 = linkDict[(spinSystemNumberA*hc+spinSystemNumberBp1)] 
                newPeaks += (linkABp1.realPeaks)
                newLinkScore += linkABp1.score
            
          elif AisJoker :
            
            spinSystemNumberB = B.spinSystemNumber
            
            if previousResB : 
            
              spinSystemBm1 = previousResB.currentSpinSystemAssigned
              
              if not spinSystemBm1.isJoker :
                
                spinSystemNumberBm1 = spinSystemBm1.spinSystemNumber
               
                linkBBm1 = previousResB.linkDict[(spinSystemNumberBm1*hc+spinSystemNumberB)]
                oldPeaks += (linkBBm1.realPeaks)
                oldLinkScore += linkBBm1.score
                
                
            if  nextResB :
            
              spinSystemBp1 = nextResB.currentSpinSystemAssigned 
              
              if not spinSystemBp1.isJoker :
                
                spinSystemNumberBp1 = spinSystemBp1.spinSystemNumber
            
                linkBBp1 = currentResidueB.linkDict[(spinSystemNumberB*hc+spinSystemNumberBp1)] 
                oldPeaks += (linkBBp1.realPeaks)
                oldLinkScore += linkBBp1.score
            
          elif BisJoker :
  
            spinSystemNumberA = A.spinSystemNumber
            
            if previousResB : 
            
              spinSystemBm1 = previousResB.currentSpinSystemAssigned
              
              if not spinSystemBm1.isJoker :
                
                spinSystemNumberBm1 = spinSystemBm1.spinSystemNumber
                
                linkABm1 = previousResB.linkDict[(spinSystemNumberBm1*hc+spinSystemNumberA)] 
                newPeaks += (linkABm1.realPeaks)
                newLinkScore += linkABm1.score
                
                
            if  nextResB :
            
              spinSystemBp1 = nextResB.currentSpinSystemAssigned 
              
              if not spinSystemBp1.isJoker :
                
                spinSystemNumberBp1 = spinSystemBp1.spinSystemNumber
  
                linkABp1 = currentResidueB.linkDict[(spinSystemNumberA*hc+spinSystemNumberBp1)] 
                newPeaks += (linkABp1.realPeaks)
                newLinkScore += linkABp1.score
            
  
  
          peakSet = set(oldPeaks+newPeaks)  
          DeltaScore = self.CcalcDeltaPeakScore(peakSet,oldPeaks,newPeaks) + newLinkScore - oldLinkScore
  
  
  
  
  
  
          if  DeltaScore == 0 :
            
            equals += 1
            
            B.currentResidueAssignment = None
            A.currentResidueAssignment = currentResidueB
            
            currentResidueB.currentSpinSystemAssigned = A
            
            score += DeltaScore
            
  
            
            for peak in peakSet :
              
              peak.degeneracy = peak.degeneracyTemp
              
  
  
          elif DeltaScore > 0 :
            
            improvements += 1
            
            B.currentResidueAssignment = None
            A.currentResidueAssignment = currentResidueB
            
            currentResidueB.currentSpinSystemAssigned = A
            
            score += DeltaScore
            
  
            
            for peak in peakSet :
              
              peak.degeneracy = peak.degeneracyTemp
              
  
              
  
          elif exp(AcceptanceConstant*DeltaScore) > cutoff() :     
            
            worse = worse + 1
  
            B.currentResidueAssignment = None
            A.currentResidueAssignment = currentResidueB
            
            currentResidueB.currentSpinSystemAssigned = A
            
            score += DeltaScore
            
  
            
            for peak in peakSet :
              
              peak.degeneracy = peak.degeneracyTemp
  
    
    self.score = score
    
    
          
  cdef double CcalcDeltaPeakScore(self,  set peakSet,list oldPeaks,list newPeaks):
      
    
    cdef double peakScoreNew = 0
    cdef double peakScoreOld = 0
    
    cdef aPeak peak 
  
    for peak in peakSet :
      
      peak.degeneracyTemp = peak.degeneracy
    
    for peak in oldPeaks :
      
      deg = peak.degeneracy
      peakScoreOld = peakScoreOld + 1.0/peak.degeneracy
      peak.degeneracyTemp -= 1
  
    for peak in newPeaks :
      
      peak.degeneracyTemp += 1
      
    for peak in newPeaks : 
      
      peakScoreNew += 1.0/peak.degeneracyTemp
      
    return peakScoreNew - peakScoreOld
  

  cdef void convertToPythonStyleDataModel(self) :
   
    cdef myDataModel DataModel
   
    cdef myChain chain
    
    cdef aResidue res
    
    cdef aSpectrum spec
    
    cdef aPeak peak
    
    cdef aDimension dim
    
    cdef list PeakList
    
    cdef list spinSystemList
    
    cdef simulatedPeak simulatedPeak
    
    cdef simulatedPeakContrib simulatedPeakContrib
    
    cdef mySpinSystem spinSystem
    
    cdef myResonance resonance
    
    cdef spinSystemLink spinSystemLink
   
    DataModel = self.DataModel
   
    chain = DataModel.myChain
  
    for res in chain.residues :
      
      res.createPythonStyleObject()
      
    chain.createPythonStyleObject()  
      

      
      
      
    for spec in DataModel.spectra :
        
      for peak in spec.peaks :
        
        for dim in peak.dimensions :
          
          dim.createPythonStyleObject()
          
        peak.createPythonStyleObject()
      
      spec.createPythonStyleObject()
      
      for peak in spec.peaks :
        
        peak.pyPeak.spectrum = spec.pySpectrum

      for peakList in spec.simulatedPeakMatrix :
        
        for simulatedPeak in peakList :
          
          for simulatedPeakContrib in simulatedPeak.simulatedPeakContribs :
          
            simulatedPeakContrib.createPythonStyleObject()
          
          simulatedPeak.createPythonStyleObject()
      
      
    for spinSystemList in DataModel.mySpinSystems.values() :
      
      for spinSystem in spinSystemList :
        
        for resonance in spinSystem.resonanceDict.values():
          
          resonance.createPythonStyleObject()
        
        spinSystem.createPythonStyleObject()
        
        
        
    for res in chain.residues :
            
      for spinSystem in res.solutions :
        
        res.pyResidue.solutions.append(spinSystem.pySpinSystem)

        
      for key, spinSystemLink in res.linkDict.items() :
        
        spinSystemLink.createPythonStyleObject()
        
        res.pyResidue.linkDict[key] = spinSystemLink.pySpinSystemLink
        
    DataModel.createPythonStyleObject()
    DataModel.pyDataModel.amountOfRepeats = self.amountOfRepeats
        
        
      
  cdef void scoreInitialAssignment(self) :
    
    cdef list residues
    
    cdef aResidue res
    
    cdef aResidue nextRes
    
    cdef double score
    
    cdef dict linkDict
    
    cdef int keyA
    
    cdef int keyB
    
    cdef int key
    
    cdef spinSystemLink link
    
    cdef list peaks
    
    cdef aPeak peak
    
    cdef double peakScore
    
    
    
    residues = self.DataModel.myChain.residues
    
    score = 0.0
    
    for res in residues :
      
      nextRes = res.nextResidue
      
      if nextRes :
        
        keyA = res.currentSpinSystemAssigned.spinSystemNumber
        keyB = nextRes.currentSpinSystemAssigned.spinSystemNumber
        
        key = keyA*self.hc + keyB
        
        linkDict = res.linkDict
          
        if key in linkDict :
          
          link = linkDict[key]
          
          peaks = link.realPeaks
          
          peakScore = 0.0
          
          for peak in peaks :
            
            peakScore = peakScore + 1.0/peak.degeneracy
            
          score += peakScore + link.score
          
          
    self.score = score
  
  

cdef class myDataModel :

  cdef list spectra

  cdef myChain myChain

  cdef dict mySpinSystems

  cdef dict previouslyAssignedSpinSystems

  cdef dict justTypedSpinSystems

  cdef dict tentativeSpinSystems

  cdef dict untypedSpinSystems

  cdef dict allSpinSystemsWithoutAssigned

  cdef dict jokerSpinSystems
  
  cdef list Hresonances 

  cdef list Nresonances
    
  cdef list CAresonances
  
  cdef list CBresonances
    
  cdef list COresonances 
    
  cdef list CXresonances 
  
  cdef list HXresonances
  
  cdef list NXresonances
    
  cdef dict residueTypeFrequencyDict 

  cdef autoAssign auto
  
  cdef object pyDataModel


  def __init__(self, autoAssign auto):
    
    self.auto = auto

    self.spectra = []

    self.myChain = myChain()
        
    self.mySpinSystems = {}
      
    self.previouslyAssignedSpinSystems = {}
        
    self.justTypedSpinSystems = {}
    
    self.tentativeSpinSystems = {}
    
    self.untypedSpinSystems = {}
    
    self.allSpinSystemsWithoutAssigned = {}
    
    self.jokerSpinSystems = {}
        
    self.Nresonances = []
    
    self.CAresonances = []
    
    self.CBresonances = []
    
    self.COresonances = []
    
    self.CXresonances = []
    
    self.HXresonances = []
    
    self.NXresonances = []
    
    self.Hresonances = []
    
    self.residueTypeFrequencyDict = {}


  def setupChain(self):

    #self.myChain.passData(self.auto, self)
    self.myChain.setupCcpnChain(self.auto.chain)
    self.myChain.setupResidues()
    self.countResidueTypeFrequency()
    self.linkResiduesTogether()

    
  cdef void countResidueTypeFrequency(self):
    
    cdef aResidue res
    
    for res in self.myChain.residues :
      
      if res.ccpCode in self.residueTypeFrequencyDict:
        
        self.residueTypeFrequencyDict[res.ccpCode] = (self.residueTypeFrequencyDict[res.ccpCode] + 1)
        
      else :
        
        self.residueTypeFrequencyDict[res.ccpCode] = 1
        
  cdef void linkResiduesTogether(self):
    
    cdef list residues
    
    cdef int i
    
    cdef aResidue res
    
    cdef aResidue nextResidue
    
    residues = self.myChain.residues
    
    
    
    for i,  res in enumerate(residues[:-1])  :
      
      nextResidue = residues[i+1]
      
      res.nextResidue = nextResidue
      
      nextResidue.previousResidue = res
    
        

  def setupSpectra(self):

    for spectrum in self.auto.selectedSpectra :

      newspectrum = aSpectrum()
      
      #newspectrum.passData(self)

      newspectrum.setupCcpnSpectrum(spectrum.ccpnSpectrum)
      
      newspectrum.setupLabellingScheme(spectrum.labellingScheme)
      
      newspectrum.setupPeakList(spectrum.peakList)

      self.spectra.append(newspectrum)
      
      
    








  cdef void createPythonStyleObject(self):
    
    cdef aSpectrum spec
    
    cdef mySpinSystem spinSystem 
    
    cdef list listWithSpinSystems
    
    cdef list listWithPySpinSystems
    
    cdef str key
    
    self.pyDataModel = pyDataModel()
    
    self.pyDataModel.chain = self.myChain.pyChain
    
    for spec in self.spectra :
    
      self.pyDataModel.spectra.append(spec.pySpectrum)
      
    for key, listWithSpinSystems in self.mySpinSystems.items():
    
      listWithPySpinSystems = []
    
      for spinSystem in listWithSpinSystems :
        
        listWithPySpinSystems.append(spinSystem.pySpinSystem)
        
      self.pyDataModel.spinSystems[key] = listWithPySpinSystems


cdef class myChain :

  cdef list residues 
    
  cdef dict residuesByCcpCode 

  cdef object ccpnChain 

  #cdef object myDataModel
  #
  #cdef autoAssign auto
  
  cdef object pyChain
  
  def __init__(self):

    self.residues = []
    
    self.residuesByCcpCode = {}

    self.ccpnChain = None
    
  def __getstate__(self):
    state = dict(self.__dict__)
    if 'ccpnChain' in state :
      del state['ccpnChain']
    if 'popup' in state :
      del state['popup']
    return state


    
    
  #def passData(self,autoAssign auto, myDataModel):
  #
  #  self.auto = auto
  #  self.myDataModel = myDataModel


  def setupCcpnChain(self, ccpnChain) :

    self.ccpnChain = ccpnChain

  def setupResidues(self):
    

    for res in self.ccpnChain.sortedResidues() :
      

      newresidue = aResidue()
      
      #newresidue.passData(self.auto, self.myDataModel)

      newresidue.setupCcpnResidue(res)
      
      newresidue.seqCode = res.seqCode
      
      newresidue.ccpCode = res.ccpCode

      newresidue.setupChainLink(self)
      
      newresidue.setupAtoms()

      self.residues.append(newresidue)
      
      if res.ccpCode in self.residuesByCcpCode :
          
        self.residuesByCcpCode[res.ccpCode].append(newresidue)
        
      else :
        
        self.residuesByCcpCode[res.ccpCode] = [newresidue]
        

    
    



  cdef void createPythonStyleObject(self):
    
    cdef aResidue res
    
    self.pyChain = pyChain()
    
    for res in self.residues :
      
      self.pyChain.residues.append(res.pyResidue)


cdef class aResidue :

  cdef aResidue previousResidue
  
  cdef aResidue nextResidue
  
  cdef object ccpnResidue
  
  cdef str ccpCode
  
  cdef myChain chain
  
  cdef list atoms
  
  cdef dict atomsByName
  
  cdef dict atomsByCcpnChemAtom
  
  cdef list solutions
  
  cdef mySpinSystem currentSpinSystemAssigned
  
  cdef int seqCode
  
  cdef dict linkDict

  #cdef object myDataModel
  #
  #cdef autoAssign auto
  
  cdef object pyResidue
  
  #cdef object userDefinedSolution
  
  

  def __init__(self):

    self.ccpnResidue = None

    self.ccpCode = None

    self.chain = None    # Parent link

    self.atoms = []
    
    self.solutions = []
    
    #self.userDefinedSolution = None
    
    self.linkDict = {}
    
    self.atomsByName = {}
    
    self.atomsByCcpnChemAtom = {}
    
    
    
  def __getstate__(self):
    state = dict(self.__dict__)
    if 'ccpnResidue' in state :
      del state['ccpnResidue']
    if 'popup' in state :
      del state['popup']
    return state
    
  #def passData(self, autoAssign auto, myDataModel):
  #
  #  self.auto = auto
  #  self.myDataModel = myDataModel

  def setupCcpnResidue(self, ccpnResidue):

    self.ccpnResidue = ccpnResidue

  def setupCcpCode(self, ccpCode):

    self.ccpCode = ccpCode

  def setupChainLink(self, chain):

    self.chain = chain

  def setupAtoms(self):
      
    for atom in self.ccpnResidue.atoms :
      
      newatom = anAtom()
      #newatom.passData(self.auto, self.myDataModel)
      newatom.setupCcpnAtom(atom)
      newatom.setupAtomName(atom.chemAtom.name)
      newatom.setupResidueLink(self)
      self.atoms.append(newatom)
      self.atomsByName[atom.chemAtom.name] = newatom
      self.atomsByCcpnChemAtom[atom.chemAtom] = newatom
      
      





  cdef void createPythonStyleObject(self):
    
      self.pyResidue = pyResidue()
      
      self.pyResidue.seqCode = self.seqCode
      
      self.pyResidue.ccpCode = self.ccpCode


cdef class anAtom :

  #cdef object myDataModel
  #
  #cdef autoAssign auto

  cdef object ccpnAtom

  cdef str atomName

  cdef aResidue residue

  cdef list assignmentPossibilityDimensions
  
  cdef dict labelInfoTemp 

  def __init__(self):

    self.assignmentPossibilityDimensions = []
    
    self.labelInfoTemp = {}
    

  #def passData(self,autoAssign auto, myDataModel):
  #
  #  self.auto = auto
  #  self.myDataModel = myDataModel

  def setupCcpnAtom(self, atom):

    self.ccpnAtom = atom

  def setupAtomName(self, atomName):

    self.atomName = atomName

  def setupResidueLink(self, residue):
  
    self.residue = residue


cdef class aSpectrum :
  
  
  cdef str name
  
  cdef object ccpnSpectrum
 
  cdef list simulatedPeakMatrix
 
  cdef list peaks
 
  cdef object thisSpectrumIsUsed

  cdef object intraresidualPeaksInThisSpectrum
  
  cdef object sequentialPeaksInThisSpectrum
  
  cdef object longRangePeaksInThisSpectrum
  
  #cdef myDataModel myDataModel
  
  cdef object peakList
 
  cdef object labellingScheme
  
  cdef object pySpectrum
 
   

  def __init__(self):
    
    self.name = None

    self.ccpnSpectrum = None
    
    self.simulatedPeakMatrix = []

    self.peaks = []

    self.labellingScheme = True
    
    self.thisSpectrumIsUsed = False
    
    self.intraresidualPeaksInThisSpectrum = False
    
    self.sequentialPeaksInThisSpectrum = False
    
    self.longRangePeaksInThisSpectrum = False
    
    self.labellingScheme = None
    
    
  def __getstate__(self):
    state = dict(self.__dict__)
    if 'ccpnSpectrum' in state :
      del state['ccpnSpectrum']
    if 'labellingScheme' in state :
      del state['labellingScheme']      
    if 'popup' in state :
      del state['popup']
    return state
    
  #def passData(self,myDataModel):
  #
  #  self.myDataModel = myDataModel


  def setupCcpnSpectrum(self, ccpnSpectrum):
    
    self.ccpnSpectrum = ccpnSpectrum
    self.name = self.ccpnSpectrum.name
    
    

    
    
  def changeSpectrumUse(self, TrueOrFalse):

    self.thisSpectrumIsUsed = TrueOrFalse
    
  def changeIntra(self, TrueOrFalse):
    
    self.intraresidualPeaksInThisSpectrum = TrueOrFalse
    
  def changeSequential(self, TrueOrFalse):
    
    self.sequentialPeaksInThisSpectrum = TrueOrFalse
    
  def changeLongRange(self, TrueOrFalse):
    
    self.longRangePeaksInThisSpectrum = TrueOrFalse



  def setupLabellingScheme(self, labellingScheme):

    self.labellingScheme = labellingScheme
    
  def setupPeakList(self,  peakList):
    
    self.peakList = peakList

  def setupPeaks(self):

    peaks = self.peakList.peaks

    for peak in peaks:
      
      newpeak = aPeak()
      #newpeak.passData(self.myDataModel)
      newpeak.setupApeak(peak)
      newpeak.setLinktoSpectrum(self)
      self.peaks.append(newpeak)


  #def setupSimulatedPeakMatrix(self):
  #  
  #  
  #
  #  
  #  i = 0
  #  
  #  residues = self.myDataModel.myChain.residues
  #
  #
  #  for residue in residues[:-1] :
  #    
  #    self.simulatedPeakMatrix.append([])

      


  cdef void createPythonStyleObject(self):
    
    cdef aPeak peak

    self.pySpectrum= pySpectrum()
    
    self.pySpectrum.name = self.name
    
    for peak in self.peaks :
      
      self.pySpectrum.peaks.append(peak.pyPeak)


cdef class aPeak :

  cdef int degeneracy
  
  cdef int degeneracyTemp
  
  cdef object ccpnPeak
  
  cdef list dimensions
  
  cdef aSpectrum spectrum
  
  cdef int serial
  
  cdef int peakListSerial
  
  #cdef autoAssign auto
  #
  #cdef myDataModel myDataModel

  cdef object intraResidual
  
  cdef object pyPeak
  
  
  def __init__(self):
    
    self.ccpnPeak = None
    
    self.dimensions = []
    
    self.degeneracy = 0
    
    self.degeneracyTemp = 0
    
    self.intraResidual = False
    

    
  #def passData(self,myDataModel):
  #
  #  self.myDataModel = myDataModel
    
  
  def setupApeak(self, ccpnPeak):
    
    self.ccpnPeak = ccpnPeak
    self.setupDimensions()
    self.checkForIntraResidualAssignment()
    self.serial = ccpnPeak.serial
    self.peakListSerial = ccpnPeak.peakList.serial

    
    
  def setLinktoSpectrum(self,  spectrum):
    
    self.spectrum = spectrum
    
  def setupDimensions(self):
    
    for dim in self.ccpnPeak.peakDims :

      dimension = aDimension()
      #dimension.passData(self.myDataModel)
      dimension.setupCcpnDim(dim)
      dimension.setupPeakLink(self)
      self.dimensions.append(dimension)
      
  def checkForIntraResidualAssignment(self):
    
    lists = []
    
    for dim in self.ccpnPeak.peakDims :
      
      spinSystems = []  
        
      for contrib in dim.peakDimContribs :
        
        if contrib.resonance.resonanceGroup :
          
          spinSystems.append(contrib.resonance.resonanceGroup)
          
      lists.append(spinSystems)
      
    intra = False
    for element in lists[0] :
     
      inAllLists = True
      
      for list in lists[1:] :
      
        if element not in list :
         
          inAllLists = False 
          
      if inAllLists :

        intra = True
        
    self.intraResidual = intra





  cdef void createPythonStyleObject(self):
    
    cdef aDimension dim
    
    self.pyPeak =pyPeak()
    
    self.pyPeak.serial = self.serial
        
    self.pyPeak.peakListSerial = self.peakListSerial
    
    for dim in self.dimensions :
      
      self.pyPeak.dimensions.append(dim.pyDimension)

 
cdef class aDimension :
  
  cdef object ccpnDim
  
  cdef double ppmValue
  
  cdef int dimNumber
  
  cdef list possibleContributions
  
  cdef list nonLabelledResonances
  
  cdef aPeak peak 
  
  #cdef myDataModel myDataModel
  #
  #cdef autoAssign auto
  
  cdef pyDimension
  
  
  
  def __init__(self):

    self.ccpnDim = None
      
    self.possibleContributions = []                         # All resonances in the resonanceList that could potentially contribute to this dimension of the peak 
    
    self.nonLabelledResonances = []                     # Here all resonances are gathered that can not contribute to the peak because of the labelling scheme. They are collected anywya to search for peaks that explicitely should NOT be there.

    self.peak = None
    
  def __getstate__(self):
    state = dict(self.__dict__)
    if 'ccpnDim' in state :
      del state['ccpnDim']
    if 'popup' in state :
      del state['popup']
    if 'possibleContributions' in state :
      del state['possibleContributions']
    if 'nonLabelledResonances' in state :
      del state['nonLabelledResonances']
    if 'peak' in state :
      del state['peak']
    
    return state
    
  #def passData(self,myDataModel):
  #
  #  self.myDataModel = myDataModel


  def setupCcpnDim(self, dim):

    self.ccpnDim = dim
    self.ppmValue = self.ccpnDim.value
    #self.dimNumber = self.ccpnDim.dim
    self.dimNumber = self.ccpnDim.dataDim.expDim.refExpDim.dim    



  def setupPeakLink(self, peak):

    self.peak = peak
    








  
  cdef void createPythonStyleObject(self) :
    
    self.pyDimension = pyDimension()
          
    self.pyDimension.ppmValue = self.ppmValue
    
    self.pyDimension.dimNumber = self.dimNumber

 
cdef class spinSystemLink :

  cdef public int score
  
  cdef public list realPeaks
  
  cdef mySpinSystem spinSystem1
  
  cdef mySpinSystem spinSystem2
  
  cdef list simulatedPeaks
  
  cdef list notFoundSimulatedPeaks
  
  cdef list simulatedPeaksThatShouldNotHaveBeenThere
  
  cdef list peaksThatShouldNotHaveBeenThere
  
  #cdef autoAssign auto
  #
  #cdef myDataModel myDataModel
  
  cdef object pySpinSystemLink
  
  def __init__(self):
    

    self.score = 0
    
    self.spinSystem1 = None
    self.spinSystem2 = None
    
    self.simulatedPeaks = []
    self.realPeaks = []
    self.notFoundSimulatedPeaks = []
    
    self.simulatedPeaksThatShouldNotHaveBeenThere = []
    self.peaksThatShouldNotHaveBeenThere = []
    
  
  cdef void createPythonStyleObject(self) :
    
    cdef aPeak peak
    
    cdef simulatedPeak simulatedPeak
    
    self.pySpinSystemLink = pySpinSystemLink()
    
    self.pySpinSystemLink.spinSystem1 = self.spinSystem1.pySpinSystem
    
    self.pySpinSystemLink.spinSystem2 = self.spinSystem2.pySpinSystem
    
    for peak in self.realPeaks :
      
      self.pySpinSystemLink.realPeaks.append(peak.pyPeak)
      
    for simulatedPeak in self.simulatedPeaks :
      
      self.pySpinSystemLink.simulatedPeaks.append(simulatedPeak.pySimulatedPeak)

 
cdef class simulatedPeak :
  
  cdef double colabelling
  
  cdef list simulatedPeakContribs
  
  cdef aSpectrum spectrum
  
  cdef object pySimulatedPeak
  
  
    
  def __init__(self):

    self.colabelling = 0.0
    
    self.simulatedPeakContribs = []
    
    self.spectrum = None
    
  cdef void createPythonStyleObject(self) :
    
    cdef simulatedPeakContrib contrib
    
    self.pySimulatedPeak = pySimulatedPeak()
    
    self.pySimulatedPeak.colabelling = self.colabelling
    
    for contrib in self.simulatedPeakContribs :
      
      self.pySimulatedPeak.simulatedPeakContribs.append(contrib.pySimulatedPeakContrib)
 
  
cdef class simulatedPeakContrib:
  
  
  cdef str ccpCode
  
  cdef str atomName
  
  cdef int dimNumber
  
  cdef int firstOrSecondResidue
  
  cdef object pySimulatedPeakContrib
    
  def __init__(self):

    self.ccpCode = None
    
    self.atomName = None
    
  cdef void createPythonStyleObject(self) :
    
    self.pySimulatedPeakContrib = pySimulatedPeakContrib()
    
    self.pySimulatedPeakContrib.ccpCode = self.ccpCode
    
    self.pySimulatedPeakContrib.atomName = self.atomName
    
    self.pySimulatedPeakContrib.dimNumber =  self.dimNumber
    
    self.pySimulatedPeakContrib.firstOrSecondResidue = self.firstOrSecondResidue


cdef class mySpinSystem :

  cdef public aResidue currentResidueAssignment
  
  cdef public int spinSystemNumber
  
  cdef public list exchangeSpinSystems
  
  cdef str ccpCode
  
  cdef dict resonanceDict
  
  cdef object ccpnResonanceGroup
  
  cdef int ccpnSeqCode
  
  cdef object isJoker
  
  cdef list solutions
  
  cdef list userDefinedSolutions
  
  cdef list tentativeCcpCodes
  
  cdef list tentativeSeqCodes
  
  cdef set allowedResidues
  
  cdef dict aminoAcidProbs
  
  #cdef autoAssign auto
  #
  #cdef myDataModel myDataModel
  
  cdef object pySpinSystem
  
  
  
  
  
  
  
  
  
  def __init__(self):

    self.ccpCode = None
    
    self.resonanceDict = {}
    
    self. ccpnResonanceGroup = None
    
    self.currentResidueAssignment = None   
    
    self.isJoker = False
    
    self.solutions = []
    
    self.userDefinedSolutions = []
    
    self.tentativeCcpCodes = []
    
    self.tentativeSeqCodes = []
    
    self.exchangeSpinSystems = []
    
    self.allowedResidues = set()                                      # This is later on going to be a Frozen Set for fast membership testing during the annealing run. If the set is empty that means everything residue is allowed, if has members, only these residues are allowed.
    
 
    
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


cdef class myResonance :
  
  cdef mySpinSystem mySpinSystem
  
  cdef double CS
  
  cdef object isotope               # String?
  
  cdef str atomName 
  
  cdef dict peakDimsLib
  
  cdef dict peakDimsLibUnlabelled
  
  cdef object ccpnResonance
  
  cdef object pyResonance
  
  
  
  def __init__(self):

    self.mySpinSystem = None
    
    self.isotope = None
    
    self.atomName = None
    
    self.peakDimsLib = {}
    
    self.peakDimsLibUnlabelled = {}                     # This is a dictionary with all dimensions of peaks that are in spectra where this resonance is not labelled. Helps to find peaks that should explicitely NOT be there.
    
    self.ccpnResonance = None
    











  cdef void createPythonStyleObject(self) :
    
    self.pyResonance = pyResonance()
    
    self.pyResonance.CS = self.CS
    
    self.pyResonance.atomName = self.atomName
    
    