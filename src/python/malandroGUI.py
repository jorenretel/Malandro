

def open_reference(argServer):
  """Descrn: Opens the macro.
     Inputs: ArgumentServer
     Output: None
  """
  print 'version...'

  program = connector(argServer.parent)
  
def lockUntillResults(function) :
  '''Decorator to lock certain functions untill results are present.
  '''

  def wrapperFunction(self,*args, **kwargs) :
    
    if self.connector.results :
  
      function(self,*args, **kwargs)
      
    else :
      
      string = 'There are no results since the annealing did not run yet.'
      showWarning('Run Annealing first', string,  parent=self)
    
  return wrapperFunction  

import Tkinter
import re
import os
import random
import math
import time

import cProfile

import cPickle
#import gzip

from sys import setrecursionlimit
from memops.universal.Io import joinPath, splitPath

from memops.general import Implementation

#from memops.general import Implementation

from memops.gui.Button import Button
from memops.gui.ButtonList import ButtonList, UtilityButtonList
from memops.gui.Canvas import Canvas
from memops.gui.CheckButton import CheckButton
from memops.gui.Color import black, hsbToRgb, inverseGrey, standardColors, hexRepr
from memops.gui.FloatEntry import FloatEntry
from memops.gui.IntEntry import IntEntry
from memops.gui.Frame import Frame
from memops.gui.LabelFrame import LabelFrame
from memops.gui.Label import Label
from memops.gui.PartitionedSelector import PartitionedSelector
from memops.gui.PulldownList import PulldownList
from memops.gui.PulldownMenu import PulldownMenu
from memops.gui.RadioButtons import RadioButtons
from memops.gui.ScrolledCanvas import ScrolledCanvas
from memops.gui.Spacer import Spacer
from memops.gui.TabbedFrame import TabbedFrame
from memops.gui.ToggleLabel import ToggleLabel
from memops.gui.FileSelect import FileSelect,  FileType
from memops.gui.MessageReporter import showWarning,  showYesNo, showMulti
from memops.gui.ScrolledGraph import ScrolledGraph
from memops.gui.Entry import Entry

from ccpnmr.analysis.core.AssignmentBasic import isAtomAssigned, getShiftLists
from ccpnmr.analysis.core.MoleculeBasic import STANDARD_ISOTOPES, unicodeGreek

from ccpnmr.analysis.core.ChemicalShiftBasic import getShiftsChainProbabilities

from ccpnmr.analysis.popups.BasePopup import BasePopup
from ccpnmr.analysis.core.MoleculeBasic import getResidueCode
from ccpnmr.analysis.core.Util import stringFromExperimentSpectrum, getSpectrumPosContourColor

from ccp.util.LabeledMolecule import getIsotopomerSingleAtomFractions, getIsotopomerAtomPairFractions,  atomPairFrac,  atomPairFractions

from memops.gui.ScrolledGraph import ScrolledGraph
from memops.gui.ScrolledDensityMatrix import ScrolledDensityMatrix
from memops.gui.ScrolledMatrix  import ScrolledMatrix

from ccp.general.Constants import standardResidueCcpCodes

from ccpnmr.analysis.macros.ArgumentServer import ArgumentServer
from ccpnmr.analysis.core.MoleculeBasic import greekSortAtomNames

from ccpnmr.analysis.core.Util import getAnalysisDataDim

from ccpnmr.analysis.core.MarkBasic import createPeakMark, createNonPeakMark

from ccpnmr.analysis.core.AssignmentBasic import assignResToDim #, assignSpinSystemResidue


import ccpnmr.analysis.core.WindowBasic as WindowBasic
from ccpnmr.analysis.core.WindowBasic import getSpectrumWindowView

import modifiedAssignmentBasic
reload(modifiedAssignmentBasic)
from modifiedAssignmentBasic import assignSpinSystemstoResidues

from src.cython.malandro import Malandro, myDataModel

import pythonStyleClasses

reload(pythonStyleClasses)

from pythonStyleClasses import *


AMINO_ACIDS = standardResidueCcpCodes['protein']

import assignMentTransferTab
reload(assignMentTransferTab)
from assignMentTransferTab import AssignMentTransferTab
import spectrumSelectionTab
reload(spectrumSelectionTab)
from spectrumSelectionTab import SpectrumSelectionTab
import saveLoadTab
reload(saveLoadTab)
from saveLoadTab import SaveLoadTab
import annealingSettingsTab
reload(annealingSettingsTab)
from annealingSettingsTab import AnnealingSettingsTab

class connector(object):
  '''This is just a little class that contains all settings the user configures in the GUI'''
  
  def __init__(self, guiParent):
    
    self.chain = None
    
    self.useAssignments = True
    
    self.useTentative = True
    
    self.selectedSpectra = None
    
    self.amountOfSteps = 0
    
    self.amountOfRepeats = 0
    
    self.shiftList = None
    
    self.nmrProject = None
    
    self.project = None
    
    self.minIsoFrac = None
    
    self.sourceName = None
    
    self.ranAnnealling = False
    
    self.ranPreCalculations = False
    
    self.results = None
    
    self.GUI = ViewAssignmentPopup(guiParent, self)
    
    self.auto = Malandro()
    
    self.addEnergyPoint = self.GUI.annealingSettingsTab.addEnergyPoint
        
  def update(self):
    
    '''
    Pulls all settings out of the GUI and gives them to the assignment algorithm written in Cython
    '''
    
    self.chain = self.GUI.annealingSettingsTab.chain
    
    self.useAssignments =self.GUI.annealingSettingsTab.useAssignmentsCheck.get()
    
    self.useTentative = self.GUI.annealingSettingsTab.useTentativeCheck.get()
    
    self.typeSpinSystems = self.GUI.annealingSettingsTab.typeSpinSystemsCheck.get()
    
    self.useDimenionalAssignments = self.GUI.annealingSettingsTab.useDimenionalAssignmentsCheck.get()
    
    self.amountOfRepeats = self.GUI.annealingSettingsTab.amountOfRepeats
    
    self.amountOfSteps = self.GUI.annealingSettingsTab.amountOfSteps
    
    self.shiftList = self.GUI.annealingSettingsTab.shiftList
    
    self.nmrProject = self.GUI.nmrProject
    
    self.project = self.GUI.project
    
    self.minIsoFrac = self.GUI.annealingSettingsTab.minIsoFrac
    
    #self.sourceName = self.GUI.sourceName
    
    self.selectedSpectra = []
    
    for spec in self.GUI.spectrumSelectionTab.specConfigList :
      
      if spec.used :
        
        self.selectedSpectra.append(spec)
            
    self.GUI.annealingSettingsTab.updateAcceptanceConstantList()
    
    self.acceptanceConstantList = self.GUI.annealingSettingsTab.acceptanceConstantList
    
    self.auto.updateSettings(self)
    
  def runAllCalculations(self):
    
    self.update()
    
    
    self.preCalculateDataModel()
    self.startAnnealing()
    
  def preCalculateDataModel(self):
    
    print 'updating'
    
    self.update()
    
    if self.checkPoint1() :
    
      #print 'setting up data model'
      
      #self.auto.setupMyDataModel()
      
      print 'precalculate data model'
      
      self.auto.preCalculateDataModel()
      
      self.ranPreCalculations = True
    
  def startAnnealing(self):
    
    if self.checkPoint2() :

      self.update()
      
      self.auto.startMonteCarlo()
    
      self.ranAnnealling = True
      
      self.getResults()
        
      self.GUI.updateResultsTab()
      self.GUI.assignTab.update()
    
  def checkPoint1(self):
    
    if self.checkSpectraAndLabelling() and self.checkChain():
      
      return True
      
    else :
    
      return False  
    
  def checkPoint2(self) :
    
    if self.ranPreCalculations :
      
      if self.GUI.annealingSettingsTab.updateAcceptanceConstantList():
        
        return True
    
    else :
      
      string = 'Without running the pre-calculations, the annealing can not run.'
      
      showWarning('Run pre-calculations first', string,  parent=self.GUI)
      
    return False  
    
  def checkSpectraAndLabelling(self) :
    
    '''
    Check whether some spectra are selected by the user and whether a reference experiment has been set
    for this experiment (needed to simulate the spectra) and if the user said the labelling scheme should
    come automatically from the labelled sample, the labelled mixture should be set.
    '''
    
    if not self.selectedSpectra :
      
      string = 'Select one or more spectra.'
        
      showWarning('No Spectra Selected', string,  parent=self.GUI)
       
      return False
      
      
  
    for spec in self.selectedSpectra :
     
      name = spec.ccpnSpectrum.experiment.name
     
      if not spec.ccpnSpectrum.experiment.refExperiment :
       
        string = name +' has no reference experiment, go to Experiment --> Experiments --> Experiment Types and select a experiment type.'
        
        showWarning('No Reference Experiment', string,  parent=self.GUI)
       
        return False
       

      mixtures = spec.ccpnSpectrum.experiment.labeledMixtures
     
      if spec.labellingScheme == True and not mixtures :
       
        string = name + ' is not connected to any labelled sample. In order to determine the labelling scheme automatically from the experiment, this has to be set. Go to Menu-Modecule-Isotope Labelling to do so.'
        
        showWarning('No Labelled Sample', string,  parent=self.GUI)
       
        return False
       
       
      elif len(mixtures) > 1 :
       
        string = name +' is connected to more than one labelled sample. This is physically impossible. Go to Menu-Modecule-Isotope Labelling if you want to change this.'

        showWarning('More than one Labelled Sample for Spectrum', string,  parent=self.GUI)
       
        return False
       
    return True
    
  def checkChain(self):
    '''
    Checks whether a chain is configured and whether this chain has residues.
    '''
   
     
    if not self.chain :
       
      string = 'Setup a molecular chain in Molecule --> Molecules'

      showWarning('No chain selected', string,  parent=self.GUI)
         
      return False
      
    if not self.chain.residues :  
      
      string = 'The selected chain does not have residues, set up the residues in Molecule --> Molecules --> Sequences'

      showWarning('No Residues', string,  parent=self.GUI)
      
      return False
    
    #stereoCcp = set()
    #for res in self.chain.residues :
    #  
    #  print res.ccpCode
    #  
    #  if res.seqCode == 82 :
    #    
    #    for atom in res.atoms :
    #      print 'Asp 82 hereeeee'
    #      print atom.name
    #      print atom.atomSet.name
    #      print atom.atomSet.atoms
    #  
    #  for atom in res.atoms :
    #    
    #    if len(atom.atomSet.atoms) > 1 :
    #      
    #      stereoCcp.add(res)
    #    
    #print 'stereo:'
    #for res in stereoCcp :
    #  
    #  print res.ccpCode
    #  
    #  for atom in res.atoms :
    #    
    #    if len(atom.atomSet.atoms) > 1 :
    #    
    #      print '---'
    #      print atom.name
    #      print atom.atomSet.name
    #      print atom.atomSet.atoms
      
         
    return True
        
  def getResults(self):
    '''
    Gets the results from the assignment algorithm.
    '''
    
    if self.ranAnnealling :

      self.results = self.auto.getResults()
      
  def updateInfoText(self,string):
    
    print string
    
    self.GUI.updateInfoText(string)
    
  def saveDataToPyc(self, fileName):
    
    if self.results :
      
      if not fileName[-4:] == '.pyc' :
        
        fileName = fileName + '.pyc'
      
      if os.path.exists(fileName) and not showYesNo('File exists', 'File "%s" exists, overwrite?' % fileName, parent=self.GUI) :
          
        return
      
      self.updateInfoText('Saving results...')
        
      writeFile = open(fileName,  'wb')
      
      #cPickle.dump(self.results,  writeFile)
      
      cPickle.dump(self.results, writeFile, protocol=2)
      
      writeFile.close()
      
      self.updateInfoText('Results succesfully saved.')
     
       
    
    else :
      
      string = 'There are no results to save at the moment.'
    
      showWarning('No results to save', string,  parent=self.GUI)
      
  def loadDataFromPyc(self, fileName):
    
    self.updateInfoText('Loading file...')
    
    self.project = self.GUI.project                                                       #TODO : fix this, getting these values from the GUI is ridiculous
    self.nmrProject = self.GUI.nmrProject
    
    results = None
    
    fileExists = False
        
    if os.path.exists(fileName) :
      
      fileExists = True
      
    elif  os.path.exists(fileName + '.pyc') :
      
      fileName = fileName + '.pyc'
      
      fileExists = True
    
    if fileExists :
      
      try :
   
        readFile = open(fileName,  'rb') 
    
        results = cPickle.load(readFile)
        
        readFile.close()
    
      except :
        
        string = 'Probably the file you try to open is not of the correct type.'
      
        showWarning('Can not open file', string,  parent=self.GUI)
        
      if results :
                
        if isinstance(results, myDataModel) :
        
          self.results = results
          self.updateInfoText('File loaded succesfully.')
          self.updateInfoText('Re-connecting to object in the ccpn analysis project...')
          self.results.connectToProject(self.project,self.nmrProject)
          self.updateInfoText('Done')
          return
          
        else :
          
          string = 'You are trying to open a file that was not created using this macro.'
      
          showWarning('Can not open file', string,  parent=self.GUI)
        
            
    else :
      
      string = 'The file you selected does not exist.'
      
      showWarning('No Such File', string,  parent=self.GUI)
      
    self.updateInfoText(' ')  


class spectrumSettings(object):
  
  '''
  A small class to temporarely store some info about spectra the user configures in the GUI.
  '''

  def __init__(self):
    
    self.ccpnSpectrum = None
    
    self.labellingScheme = True
    
    self.peakList = None
    
    self.used = False
   
  def setupLabellingScheme(self, labellingScheme):

    self.labellingScheme = labellingScheme
    
  def setupPeakList(self,  peakList):
  
    self.peakList = peakList  
    
  def changeSpectrumUse(self,  used) :
    
    self.used = used
    
    if not self.ccpnSpectrum.experiment :
      
      print 'This spectrum has not been connected to a experiment type. If you want to use this spectrum, go to experiments --> experiment and configure this.'
 
 
class ViewAssignmentPopup(BasePopup):
  '''
  The main popup that is shown when the macro is loaded.
  '''

  def __init__(self, parent, controler, *args, **kw):
    
  
    self.font = 'Helvetica 10'
    self.sFont = 'Helvetica %d'
    self.project   = parent.project
    self.guiParent = parent

    self.sourceName = 'RefDB'
    
    self.connector = controler
    
    self.selectedLinkA = None
    self.selectedLinkB = None
    self.selectedResidueA = None
    self.selectedResidueB = None
    
    self.energyDataSets = [[]]
    
    self.waiting   = False
    
    BasePopup.__init__(self, parent, title= "Assignment Suggestions", **kw)
    
  def open(self):
  
    self.updateAfter()
    BasePopup.open(self)

  def body(self, guiFrame):

    self.geometry('900x700')

    self.minIsoFrac = 0.1
    self.chain         = None
    self.spectra       = []
    self.selectedPeak = None
    self.amountOfRepeats = 10
    self.amountOfSteps = 10000
    self.acceptanceConstantList = [0, 0.1, 0.2,0.4, 0.8,1.0, 1.1, 1.2, 1.4, 1.6,2.0, 2.4,2.8, 3.2,4.0,4.5,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0]
    

    guiFrame.grid_columnconfigure(0, weight=1)
    guiFrame.grid_rowconfigure(0, weight=1)

    tabbedFrame = TabbedFrame(guiFrame,options=['Spectra Properties',  'Network Anchoring',  'Results','Bulk Transfer Assignments To Project', 'Save and Load'],
                              callback=self.toggleTab, grid=(0,0))
    self.tabbedFrame = tabbedFrame

    autoFrame,  NAFrame,  resultsFrame, assignFrame, saveFrame = tabbedFrame.frames
    
    self.spectrumSelectionTab = SpectrumSelectionTab(self,autoFrame)
    self.assignTab = AssignMentTransferTab(self, assignFrame)
    self.saveLoadTab = SaveLoadTab(self,saveFrame)
    self.annealingSettingsTab = AnnealingSettingsTab(self,NAFrame)
    
    
    self.resultsResidueNumber = 1
  
    self.waiting = False
    self.updateAfter()
    
    self.bottomButtons = UtilityButtonList(tabbedFrame.sideFrame, helpUrl='www.google.com')
    self.bottomButtons.grid(row=0, column=0, sticky = 'e')

    self.administerNotifiers(self.registerNotify)
  
  
    resultsFrame.expandGrid(5,0)

    resultTopFrame = LabelFrame(resultsFrame, text='Which results to show')
    resultTopFrame.grid(row=0, column=0, sticky='ew')
    
    self.resultsResidueNumber = 3

    texts    = [' < ']
    commands = [self.resultsPrevResidue]
    self.resultsPreviousButton = ButtonList(resultTopFrame,commands=commands, texts=texts)
    self.resultsPreviousButton.grid(row=0, column=1, sticky='nsew')    
    
    tipText = 'The Number of the residue in the sequence to display results for'
    self.resultsResidueNumberEntry = IntEntry(resultTopFrame, grid=(0,2), width=7, text=3,
                                      returnCallback=self.resultsUpdateAfterEntry,
                                      tipText=tipText)
    self.resultsResidueNumberEntry.bind('<Leave>', self.resultsUpdateAfterEntry, '+')
    
    
    texts    = [' > ']
    commands = [self.resultsNextResidue]
    self.resultsNextButton = ButtonList(resultTopFrame,commands=commands, texts=texts)
    self.resultsNextButton.grid(row=0, column=3, sticky='nsew')   
    
    
    selectCcpCodes = ['residue'] + AMINO_ACIDS
    self.resultsSelectedCcpCode = 'residue'
    
    
    tipText = 'Instead of going through the sequence residue by residue, jump directly to next amino acid of a specific type.'
    resultsSelectCcpCodeLabel     = Label(resultTopFrame, text = 'Directly jump to previous/next:', grid=(0,4))
    self.resultsSelectCcpCodePulldown  = PulldownList(resultTopFrame, callback=self.resultsChangeSelectedCcpCode, texts=selectCcpCodes,
                                         index=selectCcpCodes.index(self.resultsSelectedCcpCode), 
                                         grid=(0,5), tipText=tipText)
                                         
    self.selectedSolution = 1                               
                                         
    texts    = [' < ']
    commands = [self.resultsPrevSolution]
    self.resultsPreviousSolutionButton = ButtonList(resultTopFrame,commands=commands, texts=texts)
    self.resultsPreviousSolutionButton.grid(row=0, column=6, sticky='nsew')    
    
    tipText = 'If you rean the algorithm more than once, you can select the solution given by the different runs.'
    self.resultsSolutionNumberEntry = IntEntry(resultTopFrame, grid=(0,7), width=7, text=1,
                                      returnCallback=self.solutionUpdateAfterEntry,
                                      tipText=tipText)
    self.resultsSolutionNumberEntry.bind('<Leave>', self.solutionUpdateAfterEntry, '+')
    
    
    texts    = [' > ']
    commands = [self.resultsNextSolution]
    self.resultsNextSolutionButton = ButtonList(resultTopFrame,commands=commands, texts=texts)
    self.resultsNextSolutionButton.grid(row=0, column=8, sticky='nsew')   
    
    
    
    texts    = ['adopt this solution']
    commands = [self.adoptSolution]
    self.adoptButton = ButtonList(resultTopFrame,commands=commands, texts=texts)
    self.adoptButton.grid(row=0, column=9, sticky='nsew')  
    

    
    
    
    resultsSecondFrame = Frame(resultsFrame)#LabelFrame(resultsFrame, text='Spin Systems')
    resultsSecondFrame.grid(row=2, column=0, sticky='nsew')
    
    resultsSecondFrame.grid_columnconfigure(0,  weight=1) 
    resultsSecondFrame.grid_columnconfigure(1,  weight=1)        
    resultsSecondFrame.grid_columnconfigure(2,  weight=1)    
    resultsSecondFrame.grid_columnconfigure(3,  weight=1)    
    resultsSecondFrame.grid_columnconfigure(4,  weight=1) 

    headingList = [ '#',  '%']

    tipTexts = [ 'Spinsystem number {} indicates serial of the spinsystem. If the spinsystem was already assigned to a residue, the residue number is shown aswell',  'percentage of the solutions that connected this spinsystem to this residue']

    editWidgets = [None, None]

    editSetCallbacks = [None, None]
    
    self.displayResultsTables = []
    self.residueLabels = []
    
    # Sorry for this one, but it saves a lot of lines of code with arbitrary function names. createCallBackFunction takes one argument 'tableNumber'
    # and returns a function that also takes one archument 'spinSystem'. This function is called in the end when the user
    # clicks on the table to select a spinsystem. What this function does is calling self.selectSpinSystem(tableNumber,spinSystem).
    # In this way it is possible to create five tables in a loop but still tell self.selectSpinSystem from who the call came. 
    createCallbackFunction = lambda tableNumber : (lambda spinSystem :  self.selectSpinSystem(tableNumber, spinSystem) )
    
    for i in range(5) :
      
      label = Label(resultsSecondFrame, text='residue')
      label.grid(row=0, column=i)
      
      editGetCallbacks = [createCallbackFunction(i)]*3
      
      displayResultsTable = ScrolledMatrix(resultsSecondFrame, headingList=headingList,
                                       editWidgets=editWidgets, multiSelect=False,
                                       editGetCallbacks=editGetCallbacks,
                                       editSetCallbacks=editSetCallbacks,
                                       tipTexts=tipTexts)
      displayResultsTable.grid(row=2, column=i, sticky='nsew')
      displayResultsTable.sortDown = False
      
      self.residueLabels.append(label)
      self.displayResultsTables.append(displayResultsTable)
      
    resultsFirstFrame = Frame(resultsSecondFrame) #LabelFrame(resultsFrame, text='Sequence Fragment')
    resultsFirstFrame.grid(row=1, column=0, sticky='ew', columnspan=5)
    
    resultsFirstFrame.grid_rowconfigure(0, weight=1)
    resultsFirstFrame.grid_rowconfigure(1, weight=1)
    resultsFirstFrame.grid_columnconfigure(0,  weight=1)
    
    texts    = [' res 1 ',  ' link ',  ' res 2 ',  ' link ', ' res 3 ', ' link ', ' res 4 ', ' link ', ' res 5 ']
    commands = [lambda:self.selectRelativeResidue(1,True), lambda: self.selectLink(1,True), lambda:self.selectRelativeResidue(2,True), lambda: self.selectLink(2,True), lambda:self.selectRelativeResidue(3,True), lambda: self.selectLink(3,True), lambda:self.selectRelativeResidue(4,True), lambda: self.selectLink(4,True), lambda:self.selectRelativeResidue(5,True)]
    self.sequenceButtons = ButtonList(resultsFirstFrame,commands=commands, texts=texts)
    self.sequenceButtons.grid(row=0, column=0, sticky='nsew')   
    

    
    for n, button in enumerate(self.sequenceButtons.buttons) :
      
      if n%2 :
      
        button.grid(column=n, sticky='ns')
        
        self.sequenceButtons.grid_columnconfigure(n, weight=0)
        
      else :
      
        self.sequenceButtons.grid_columnconfigure(n, uniform=2)  
      


    spacer = Spacer(resultsFirstFrame)
    spacer.grid(row=1, column=0, sticky='nsew') 
    
    texts    = [' res 1 ',  ' link ',  ' res 2 ',  ' link ', ' res 3 ', ' link ', ' res 4 ', ' link ', ' res 5 ']
    commands = commands = [lambda:self.selectRelativeResidue(1,False), lambda: self.selectLink(1,False), lambda:self.selectRelativeResidue(2,False), lambda: self.selectLink(2,False), lambda:self.selectRelativeResidue(3,False), lambda: self.selectLink(3,False), lambda:self.selectRelativeResidue(4,False), lambda: self.selectLink(4,False), lambda:self.selectRelativeResidue(5,False)]
    self.sequenceButtonsB = ButtonList(resultsFirstFrame,commands=commands, texts=texts)
    self.sequenceButtonsB.grid(row=2, column=0, sticky='nsew')    
 

    for n, button in enumerate(self.sequenceButtonsB.buttons) :
      
      if n%2 :
      
        button.grid(column=n, sticky='ns')
        
        self.sequenceButtonsB.grid_columnconfigure(n, weight=0)
        
      else :
      
        self.sequenceButtonsB.grid_columnconfigure(n, uniform=2)  
    
    resultsFrame.grid_rowconfigure(3, weight=2)
    
    resultsThirdFrame = Frame(resultsFrame)
    resultsThirdFrame.grid(row=3, column=0, sticky='nsew')
    
    
    resultsThirdFrame.grid_rowconfigure(0, weight=1)
    resultsThirdFrame.grid_columnconfigure(0,  weight=1) 
   
    
    tabbedFrameB = TabbedFrame(resultsThirdFrame,options=['Peaks', 'Spin System'],
                              callback=self.toggleTab, grid=(0,0))
    self.tabbedFrameB = tabbedFrame

    PeakFrame, SpinSystemFrame = tabbedFrameB.frames
    

    SpinSystemFrame.grid_rowconfigure(0, weight=1)
    PeakFrame.grid_rowconfigure(1, weight=1)

    SpinSystemFrame.grid_columnconfigure(0,  weight=1)    
    PeakFrame.grid_columnconfigure(0,  weight=1)    
        
    headingList = [ 'residue','assigned to in project','user defined sequence', 'selected annealing result','%']
    
    tipTexts = [None, None, None, None, None]
    
    editWidgets =[None, None, None, None, None]

    editGetCallbacks = [None, None, None, None, None]

    editSetCallbacks =[None, None, None, None, None]
 

    self.spinSysTable = ScrolledMatrix(SpinSystemFrame,headingList=headingList,
                                       editWidgets=editWidgets, multiSelect=False,
                                       editGetCallbacks=editGetCallbacks,
                                       editSetCallbacks=editSetCallbacks,
                                       tipTexts=tipTexts)
    self.spinSysTable.grid(row=0, column=0, sticky='nsew')
    
    buttonFrameinPeakFrame = Frame(PeakFrame)
    buttonFrameinPeakFrame.grid(sticky='ew')

    self.findButton  = Button(buttonFrameinPeakFrame, text=' Go to Peak ',
                              borderwidth=1, padx=2, pady=1,command=self.findPeak,
                              tipText='Locate the currently selected peak in the specified window')
                              
    self.findButton.grid(row=0, column=0, sticky='e')
    
    label = Label(buttonFrameinPeakFrame,text='in window:')
    
    label.grid(row=0, column=1, sticky='w')

    self.windowPulldown = PulldownList(buttonFrameinPeakFrame, callback=self.selectWindowPane,
                                      tipText='Choose the spectrum window for locating peaks or strips')
                                    
    self.windowPulldown.grid(row=0, column=2, sticky='w')
    
    self.assignSelectedPeaksButton = Button(buttonFrameinPeakFrame, text='Assign Resonances to Peak(s)',
                              borderwidth=1, padx=2, pady=1,command=self.assignSelectedPeaks,
                              tipText='Assign resonances to peak dimensions, this of course only works when the peak is found in the spectrum.')

    self.assignSelectedPeaksButton.grid(row=0, column=3, sticky='ew')
    
    self.assignSelectedSpinSystemsToResiduesButton = Button(buttonFrameinPeakFrame, text='Assign Spinsystems to Residues',
                              borderwidth=1, padx=2, pady=1,command=self.assignSelectedSpinSystemsToResidues,
                              tipText='Assign spinsystems to residues')

    self.assignSelectedSpinSystemsToResiduesButton.grid(row=0, column=4, sticky='ew')

    headingList = [ '#','spectrum','Dim1', 'Dim2','Dim3', 'c.s. dim1','c.s. dim2','c.s. dim3', 'colabelling']
    
    tipTexts = ['Peak number, only present when the peak was actually found in the spectrum.',
                'Name of the spectrum',
                'Name of atomSet measured in this dimension. Dimension number corresponds to Ref Exp Dim as indicated by going in the main menu to Experiment-->Experiments-->Experiment Type',
                'Name of atomSet measured in this dimension. Dimension number corresponds to Ref Exp Dim as indicated by going in the main menu to Experiment-->Experiments-->Experiment Type',
                'Name of atomSet measured in this dimension. Dimension number corresponds to Ref Exp Dim as indicated by going in the main menu to Experiment-->Experiments-->Experiment Type',
                'Chemical Shift', 'Chemical Shift', 'Chemical Shift',
                'Colabbeling fraction over all nuclei that are on the magnetization transfer pathway during the experiment that gave rise to the peak, including visited nuclei that were not measured in any of the peak dimensions']
    
    #editWidgets = [None, None, None, None, None, None, None, None, None]
    
    editGetCallbacks = [None, None, None, None, None, None, None, None, None]
    
    #editGetCallbacks = [self.selectPeak, self.selectPeak, self.selectPeak, self.selectPeak, self.selectPeak, self.selectPeak, self.selectPeak, self.selectPeak, self.selectPeak]

    editSetCallbacks = [None, None, None, None, None, None, None, None, None]
    
    self.displayPeakTable = ScrolledMatrix(PeakFrame,headingList=headingList,
                                       editWidgets=editWidgets, multiSelect=True,
                                       editGetCallbacks=editGetCallbacks,
                                       editSetCallbacks=editSetCallbacks,
                                       tipTexts=tipTexts)
    self.displayPeakTable.grid(row=1, column=0, sticky='nsew')
    
    self.infoLabel = Label(guiFrame,text=' ')
    
    self.infoLabel.grid(row=2, column=0, sticky='w')
    
    

    self.windowPane = 'None'
    self.updateWindows()

  @lockUntillResults              
  def showResults(self):
      
    self.updateResultsTable()
   
  @lockUntillResults  
  def selectLink(self,number, topRow) :
    
    if topRow :
      
      self.selectedResidueA = None
      self.selectedResidueB = None
      self.selectedLinkA = number
      self.selectedLinkB = None
      
    else :
      
      self.selectedResidueA = None
      self.selectedResidueB = None
      self.selectedLinkA = None
      self.selectedLinkB = number
      
    self.updateButtonHighLights()
    self.updateLink()
    
  @lockUntillResults  
  def selectRelativeResidue(self,number,topRow) :
    
    if topRow :
      
      self.selectedResidueA = number
      self.selectedResidueB = None
      self.selectedLinkA = None
      self.selectedLinkB = None
      
    else :
      
      self.selectedResidueA = None
      self.selectedResidueB = number
      self.selectedLinkA = None
      self.selectedLinkB = None
      
    self.updateButtonHighLights()
    self.updateLink()

  def updateLink(self) :
    '''
    Checks for any selected link (self.selectedLinkA or self.selectedLinkB) and calls
    updatePeakTable with the correct residue Object and spinsystem Objects.
    '''
    
    number = self.selectedLinkA or self.selectedLinkB or self.selectedResidueA or self.selectedResidueB
    
    if not number:
      
      self.emptyPeakTable()
      return
    
    DataModel = self.connector.results
    resNumber = self.resultsResidueNumber
    chain = DataModel.myChain
    residues = chain.residues    
    solutionNumber = self.selectedSolution-1
    
    
    
    if self.selectedResidueA:
      
      res = residues[resNumber-4 + number]
      spinSystem = res.solutions[solutionNumber]
      
      self.selectedLink = None
      if res and spinSystem :
        self.selectedLink = res.getIntraLink(spinSystem)

      #self.updatePeakTableIntra(res, spinSystem)
      self.updateSpinSystemTable(spinSystem)
      
    elif self.selectedResidueB :
      
      res = residues[resNumber-4 + number]
      spinSystem = res.userDefinedSolution
      
      self.selectedLink = None
      if res and spinSystem :
        self.selectedLink = res.getIntraLink(spinSystem)
      
      #self.updatePeakTableIntra(res, spinSystem)
      self.updateSpinSystemTable(spinSystem)
       
    elif self.selectedLinkA:
      
      resA = residues[resNumber-4 + number]
      resB = residues[resNumber-3 + number]
      
      spinSystemA = resA.solutions[solutionNumber]
      spinSystemB = resB.solutions[solutionNumber]
      
      self.selectedLink = None
      if resA and spinSystemA and spinSystemB :
        self.selectedLink = resA.getLink(spinSystemA,spinSystemB)
        
      #self.updatePeakTable(resA, spinSystemA, spinSystemB)
      
    elif self.selectedLinkB : #and resA.userDefinedSolution and resB.userDefinedSolution:  
        
      resA = residues[resNumber-4 + number]
      resB = residues[resNumber-3 + number]
      
      spinSystemA = resA.userDefinedSolution
      spinSystemB = resB.userDefinedSolution
      
      self.selectedLink = None
      if resA and spinSystemA and spinSystemB :
        self.selectedLink = resA.getLink(spinSystemA,spinSystemB)
      
      #self.updatePeakTable(resA, spinSystemA, spinSystemB)
      
    self.updatePeakTable()
    
  def emptyPeakTable(self) :
    
    self.displayPeakTable.update(objectList=[],textMatrix=[], colorMatrix=[])

  def updatePeakTable(self):
    '''
    Updates the peak table to show the peaks that are found for a sequencial pair of
    spinsystems A and B. If there is not a linkobject found for spinsystems A and B the
    table is emptied. Also sets the selected peak to None.
    '''
    
    link = self.selectedLink
    
    if not link :
      
      self.emptyPeakTable()
      
    else :
      
      resA, resB = link.getResidues()
      spinSystemA, spinSystemB = link.getSpinSystems()
      
      data = []
      
      objectList = []
      
      #peakLinks = link.getPeakLinks()

      peakLinks = link.getAllPeakLinks()

      for peakLink in peakLinks : #link.getPeakLinks() :

        realPeak = peakLink.getPeak()

        simPeak = peakLink.getSimulatedPeak()
        
        oneRow = [None, None, None, None, None, None, None, None, None]
        
        oneRow[1] = simPeak.getSpectrum().name

        oneRow[8] = simPeak.colabelling

        for simulatedPeakContrib in simPeak.getContribs() :

          atomName = simulatedPeakContrib.getAtomName()

          ccpCode = simulatedPeakContrib.getCcpCode()

          dimNumber = simulatedPeakContrib.getDimNumber()

          if resA is simulatedPeakContrib.getResidue() :

            spinSystemDescription = spinSystemA.getDescription(noSerialWhenSeqCodeIsPresent=True)

          else :
            
            spinSystemDescription = spinSystemB.getDescription(noSerialWhenSeqCodeIsPresent=True)
            
          oneRow[dimNumber+1] = '%s %s' %(spinSystemDescription, atomName) #ccpCode + '{' + str(spinSystemNumber) +'} ' + atomName
          
        if realPeak :

          oneRow[0] = realPeak.getSerial()

          for dim in realPeak.getDimensions() :
              
            oneRow[dim.getDimNumber()+4] =  dim.getChemicalShift()
            
        else :
          
          for resonance, simulatedPeakContrib in zip(peakLink.getResonances(), simPeak.getContribs()) :
            
            if resonance :
            
              oneRow[simulatedPeakContrib.getDimNumber()+4] = resonance.getChemicalShift()
              
            else :
              
              oneRow[simulatedPeakContrib.getDimNumber()+4] = '?'
          
        data.append(oneRow)
        objectList.append(peakLink)
      
      self.displayPeakTable.update(objectList=objectList,textMatrix=data)

      self.selectedPeak = None  

  def findPeak(self):
    
    if not self.windowPane :
      
      return
    
    selectedPeakLinks = self.displayPeakTable.currentObjects
    
    if not selectedPeakLinks :
      
      self.updateInfoText('Please select a peak first.')
      return
    
    if len(selectedPeakLinks) > 1 :
      
      self.updateInfoText('Can only go to one peak at a time.')
      return
    
    selectedPeakLink = selectedPeakLinks[0]
    selectedPeak = selectedPeakLink.getPeak()
    
    if selectedPeak :

      ccpnPeak = selectedPeak.getCcpnPeak()
      createPeakMark(ccpnPeak, lineWidth=2.0)
      
      windowFrame = self.windowPane.getWindowFrame()
      windowFrame.gotoPeak(ccpnPeak)
      
    else :
      
      simPeak = selectedPeakLink.getSimulatedPeak()
      spectrum = simPeak.getSpectrum()
      ccpnSpectrum = spectrum.getCcpnSpectrum()
      view = getSpectrumWindowView(self.windowPane, ccpnSpectrum)
      
      if not view :
        
        self.updateInfoText('This peak cannot be displayed in the window you chose.')
      
      axisMappingByRefExpDimNumber = {}
      
      for axisMapping in view.axisMappings :
        
        refExpDimNumber = axisMapping.analysisDataDim.dataDim.expDim.refExpDim.dim
      
        axisMappingByRefExpDimNumber[refExpDimNumber] = axisMapping
        
        
      positionToGoTo = {}
      markPosition = []
      axisTypes = []
      
      for resonance, contrib in zip(selectedPeakLink.getResonances(), simPeak.getContribs()) :
        
        dimNumber = contrib.getDimNumber()
        
        axisMapping = axisMappingByRefExpDimNumber.get(dimNumber)
        
        label = axisMapping.label
          
        if resonance :
          
          axisType = axisMapping.axisPanel.axisType
          chemicalShift = resonance.getChemicalShift()
          
          positionToGoTo[label] = chemicalShift
          markPosition.append(chemicalShift)
          axisTypes.append(axisType)
          
        else :                                                    # Not drawing a mark at this chemical shift, just hoovering to the good region in the spectrum 
          
          ccpCode = contrib.getResidue().getCcpCode()
          atomName = contrib.getAtomName()
          
          medianChemicalShift = self.getMedianChemicalShift(ccpCode, atomName)
          
          if medianChemicalShift :
          
            positionToGoTo[label] = medianChemicalShift
        
      
      if positionToGoTo :
        
        windowFrame = self.windowPane.getWindowFrame()
        windowFrame.gotoPosition(positionToGoTo)
        
      if markPosition :  
        
        createNonPeakMark(markPosition, axisTypes)
        
  def assignSelectedPeaks(self) :

    selectedPeakLinks = self.displayPeakTable.currentObjects
    
    for pl in selectedPeakLinks :
      
      peak = pl.getPeak()
      
      if peak :
        
        for resonance, dimension in zip(pl.getResonances(), peak.getDimensions()) :
          
          ccpnResonance = resonance.getCcpnResonance()
          ccpnDimension = dimension.getCcpnDimension()
          
          assignResToDim(ccpnDimension, ccpnResonance)
          
  def assignSelectedSpinSystemsToResidues(self) :
    
    link = self.selectedLink
    
    residues = link.getResidues()
    spinSystems = link.getSpinSystems()
    
    ccpnSpinSystems = []
    ccpnResidues = []
    
    for spinSys,res in zip(spinSystems, residues) :
      
      if spinSys and res:
        
        ccpnSpinSystems.append(spinSys.getCcpnResonanceGroup())
        ccpnResidues.append(res.getCcpnResidue())
    
    assignSpinSystemstoResidues(ccpnSpinSystems, ccpnResidues, guiParent=self)
    
    self.updateResultsTopRowButtons()    
    self.updateResultsBottomRowButtons()
    self.updateResultsTable()
    self.updatePeakTable() 

  def getMedianChemicalShift(self, ccpCode, atomName) :
    
    nmrRefStore = self.project.findFirstNmrReferenceStore(molType='protein',ccpCode=ccpCode)

    chemCompNmrRef = nmrRefStore.findFirstChemCompNmrRef(sourceName='RefDB')
    
    chemCompVarNmrRef = chemCompNmrRef.findFirstChemCompVarNmrRef(linking='any',descriptor='any')
    
    if chemCompVarNmrRef :
      
      chemAtomNmrRef = chemCompVarNmrRef.findFirstChemAtomNmrRef(name=atomName)
      
      if chemAtomNmrRef :
        
        distribution   = chemAtomNmrRef.distribution
        
        maxIndex = max([(value,index) for index, value in enumerate(distribution)])[1]
    
        return chemAtomNmrRef.refValue + chemAtomNmrRef.valuePerPoint*(maxIndex-chemAtomNmrRef.refPoint)
        
    return None
  
  def selectWindowPane(self, windowPane):
  
    if windowPane is not self.windowPane:
      self.windowPane = windowPane

  def updateWindows(self):
  
    index = 0
    windowPane = None
    windowPanes = []
    names = []
    peakList = None
    tryWindows = WindowBasic.getActiveWindows(self.project)
    
    windowData = []
    getName = WindowBasic.getWindowPaneName
    
    for window in tryWindows:
      for windowPane0 in window.spectrumWindowPanes:
        #if WindowBasic.isSpectrumInWindowPane(windowPane0, spectrum):
        windowData.append( (getName(windowPane0), windowPane0) )
      
      windowData.sort()
      names = [x[0] for x in windowData]
      windowPanes = [x[1] for x in windowData]
    
    if windowPanes:
      if windowPane not in windowPanes:
        windowPane = windowPanes[0]
        
      index = windowPanes.index(windowPane)
        
    else:
      windowPane = None
    
    self.selectWindowPane(windowPane)
    
    self.windowPulldown.setup(names, windowPanes, index)
     
  def selectSpinSystem(self, number, spinSystem):
    
    res = self.connector.results.myChain.residues[self.resultsResidueNumber-3 + number]
    
    oldSpinSystemForResidue = res.userDefinedSolution
    
    if oldSpinSystemForResidue and res.getSeqCode() in oldSpinSystemForResidue.userDefinedSolutions :
    
      oldSpinSystemForResidue.userDefinedSolutions.remove(res.getSeqCode())
    
    res.userDefinedSolution = spinSystem
    
    spinSystem.userDefinedSolutions.append(res.getSeqCode())
    
    #self.updateSpinSystemTable(spinSystem)
    self.updateLink()
    self.updateResultsBottomRowButtons()
    
  def updateSpinSystemTable(self, spinSystem):
    
    if not spinSystem :
      
      self.emptySpinSystemTable()
      return
    
    DataModel = self.connector.results
    
    residues = DataModel.myChain.residues
    
    data = []
    colorMatrix = []
    
    for residueNumber in spinSystem.allowedResidues :
      
      residue = residues[residueNumber-1]
        
      oneRow = []
      
      oneRowColor = []
      
      string = str(residue.getSeqCode()) + ' ' + residue.getCcpCode()

      oneRow.append(string)
      
      ccpnSeqCode = spinSystem.getCcpnSeqCode()
          
      if ccpnSeqCode and ccpnSeqCode == residue.getSeqCode() :                          # Assigned in the project to this residue
       
        oneRow.append('x') 
        
      else :
        
        oneRow.append(None)
        
      if residue.getSeqCode() in spinSystem.userDefinedSolutions :                                # The user selected this res for this spinsystem (could be more than one res for which this happens)
        
        oneRow.append('x')
        
      else :
        
        oneRow.append(None)
        
      if residue.solutions[self.selectedSolution-1] == spinSystem :
        
        oneRow.append('x')
        
      else :
        
        oneRow.append(None)
      
      if spinSystem.solutions :
        
        percentage = spinSystem.solutions.count(residueNumber)/float(len(spinSystem.solutions))*100.0
        
      else :
        
        percentage = 0
        
      oneRow.append(int(percentage))
      color = self.pickColorByPercentage(percentage)
      oneRowColor = [color]*5
         
      data.append(oneRow)    
      colorMatrix.append(oneRowColor)
          
    self.spinSysTable.update(objectList=data,textMatrix=data, colorMatrix=colorMatrix)

    
    self.spinSysTable.sortDown = False
    self.spinSysTable.sortLine(-1,  noUpdate=True)
    
  def emptySpinSystemTable(self) :
    
    self.spinSysTable.update(objectList=[],textMatrix=[], colorMatrix=[])
  
  @lockUntillResults   
  def adoptSolution(self):
    
    DataModel = self.connector.results
    
    selectedSolution = self.selectedSolution
    
    for res in DataModel.myChain.residues :
      
      res.userDefinedSolution = res.solutions[selectedSolution]
     
    self.updateLink()
    self.updateResultsBottomRowButtons()
  
  @lockUntillResults    
  def resultsPrevSolution(self):

    if self.selectedSolution != 1:
      self.selectedSolution = self.selectedSolution - 1
      self.resultsSolutionNumberEntry.set(self.selectedSolution)
      
      self.updateLink()
      self.updateResultsTopRowButtons()
 
  @lockUntillResults     
  def resultsNextSolution(self):
    
    amountOfRepeats = len(self.connector.results.myChain.residues[0].solutions)
    
    if self.selectedSolution < amountOfRepeats:
      self.selectedSolution = self.selectedSolution + 1
      self.resultsSolutionNumberEntry.set(self.selectedSolution)
      
      self.updateLink()
      self.updateResultsTopRowButtons()
  
  @lockUntillResults
  def resultsPrevResidue(self):
    
    residues = self.connector.results.myChain.residues
    #chainLength = len(residues)
    
    new_value = self.resultsResidueNumber
    if self.resultsSelectedCcpCode == 'residue' :
      if self.resultsResidueNumber != 3 :
        new_value = self.resultsResidueNumber - 1
    else :
      for res in residues :
        if res.getSeqCode == self.resultsResidueNumber :
          break
        elif res.getCcpCode() == self.resultsSelectedCcpCode :
          new_value = res.seqCode
          if new_value < 3 :
            new_value = 3
    if self.resultsResidueNumber != new_value :
      self.resultsResidueNumber = new_value
      
      self.resultsResidueNumberEntry.set(self.resultsResidueNumber)
      
      self.updateLink()
      self.updateResultsTopRowButtons()
      self.updateResultsBottomRowButtons()
      self.updateResultsTable()
      self.updateResidueLabels()
      
  @lockUntillResults
  def resultsNextResidue(self):
    
    residues = self.connector.results.myChain.residues
    chainLength = len(residues)
    
    new_value = self.resultsResidueNumber
    if self.resultsSelectedCcpCode == 'residue' :
      if self.resultsResidueNumber != chainLength-2:
        new_value = self.resultsResidueNumber + 1
    else :
      for res in residues[(self.resultsResidueNumber):] :
        if res.getCcpCode() == self.resultsSelectedCcpCode :
          new_value = res.seqCode
          if new_value > chainLength-2 :
            new_value = chainLength-2
          break
    if self.resultsResidueNumber != new_value :
      self.resultsResidueNumber = new_value
      self.resultsResidueNumberEntry.set(self.resultsResidueNumber)
      
      self.updateLink()
      self.updateResultsTopRowButtons()
      self.updateResultsBottomRowButtons()
      self.updateResultsTable()
      self.updateResidueLabels()
  
  def resultsChangeSelectedCcpCode(self, ccpCode):
  
    self.resultsSelectedCcpCode = ccpCode
  
  @lockUntillResults
  def resultsUpdateAfterEntry(self, event=None):
    '''
    Update for entry of residue number in strip plots
    '''

    residues = self.connector.results.myChain.residues

    value = self.resultsResidueNumberEntry.get()
    if value == self.resultsResidueNumber:
      return 
    else :
      self.resultsResidueNumber = value
    if value < 3:
      self.resultsResidueNumberEntry.set(3)
      self.resultsResidueNumber = 3

    elif value > len(residues)-2:
      self.resultsResidueNumber = len(residues)-2
      self.resultsResidueNumberEntry.set(self.resultsResidueNumber)

    else :
      self.resultsResidueNumberEntry.set(self.resultsResidueNumber)
      
    self.updateLink()
    self.updateResultsTopRowButtons()
    self.updateResultsBottomRowButtons()
    self.updateResultsTable()
    self.updateResidueLabels()

  @lockUntillResults
  def solutionUpdateAfterEntry(self, event=None):
    '''
    Update for entry of residue number in strip plots
    '''
  
    Nsolutions = len(self.connector.results.myChain.residues[0].solutions)

    value = self.resultsSolutionNumberEntry.get()
    if value == self.selectedSolution:
      return 
    else :
      self.selectedSolution = value
    if value < 1:
      self.resultsSolutionNumberEntry.set(1)
      self.selectedSolution = 1
    elif value > Nsolutions:
      self.selectedSolution = Nsolutions
      self.resultsSolutionNumberEntry.set(self.selectedSolution)
    else :
      self.resultsSolutionNumberEntry.set(self.selectedSolution)
      
    self.updateLink()
    self.updateResultsTopRowButtons()      
    
  def updateResultsTab(self) :
    
    self.updateLink()
    self.updateResidueLabels()
    self.updateResultsTable()
    self.updateResultsTopRowButtons()
    self.updateResultsBottomRowButtons()
    
  def updateResultsTable(self):
    
    resNumber = self.resultsResidueNumber    
    
    DataModel = self.connector.results
    
    chain = DataModel.myChain
    
    residues = chain.residues
  
    resA = residues[resNumber -3]
    resB = residues[resNumber -2]
    resC = residues[resNumber - 1]
    resD = residues[resNumber]
    resE = residues[resNumber +1]
    
    resList = [resA, resB, resC, resD, resE]
    tableList = self.displayResultsTables
    
    
    for res,  table in zip(resList, tableList) :
      
      ccpCode = res.ccpCode
      
      spinSystemsWithThisCcpCode = DataModel.getSpinSystems()[ccpCode]
      
      data = []
      colorMatrix = []
      objectList = []
      
      jokers = []
      realSpinSystems = []
      
      for spinSys in spinSystemsWithThisCcpCode :
        
        if spinSys.getIsJoker() :
          
          jokers.append(spinSys)
          
        else :
        
          realSpinSystems.append(spinSys)  
      
      for spinsys in realSpinSystems :
        
        oneRow = []
        oneRowColor = []
        
        spinSystemInfo = spinsys.getDescription() #self.getStringDescriptionOfSpinSystem(spinsys)
          
        oneRow.append(spinSystemInfo)
          
        assignmentPercentage = int(float(res.solutions.count(spinsys)) / len(res.solutions) * 100.0)
        
        oneRow.append(assignmentPercentage)
        
        objectList.append(spinsys)
        
        color = self.pickColorByPercentage(assignmentPercentage)
        
        oneRowColor = [color, color]
        
        data.append(oneRow)
        colorMatrix.append(oneRowColor)
        
        
      if jokers :
        
        oneRow = ['Joker']
        
        NumberOfAssignmentsToJoker = 0
        
        for spinSys in jokers :
        
          NumberOfAssignmentsToJoker += res.solutions.count(spinSys)
        
        assignmentPercentage = int(float(NumberOfAssignmentsToJoker) / len(res.solutions) * 100.0)
        
        oneRow.append(assignmentPercentage)
        
        color = self.pickColorByPercentage(assignmentPercentage)
        
        oneRowColor = [color, color]
        
        data.append(oneRow)
        colorMatrix.append(oneRowColor)
        objectList.append(jokers[0])
      
      
      percentages = [datapoint[1] for datapoint in data]
      
      tableData = sorted(zip(percentages,data,objectList,colorMatrix), reverse=True)
      
      percentage, data, objectList, colorMatrix = zip(*tableData)
        
      table.update(objectList=objectList,textMatrix=data,colorMatrix=colorMatrix)

  def updateResidueLabels(self) :
        
    resList = self.getCurrentlyDisplayedResidues()
    labels = self.residueLabels
    
    for residue, label in zip(resList, labels) :
      
      text = str(residue.getSeqCode()) + ' ' + residue.getCcpCode()
      
      label.set(text)

  def updateResultsTopRowButtons(self) :
    
    resList = self.getCurrentlyDisplayedResidues()
    
    buttons = self.sequenceButtons.buttons[::2]
      
    for button,  res in zip(buttons, resList) :
      
      spinsys = res.solutions[self.selectedSolution-1]
      
      text = spinsys.getDescription(noSerialWhenSeqCodeIsPresent=False) #str(res.getSeqCode()) + ' ' + res.getCcpCode() + ': ' + spinsys.getDescription() #self.getStringDescriptionOfSpinSystem(spinsys)
      
      button.config(text=text)
      
  def updateResultsBottomRowButtons(self) :
    
    resList = self.getCurrentlyDisplayedResidues()
      
    buttons = self.sequenceButtonsB.buttons[::2]
    
    for button,  res in zip(buttons, resList) :
      
      if res.userDefinedSolution :
        
        selectedSpinSystem = res.userDefinedSolution
        text = selectedSpinSystem.getDescription(noSerialWhenSeqCodeIsPresent=False) #self.getStringDescriptionOfSpinSystem(selectedSpinSystem)
        
        if len(selectedSpinSystem.userDefinedSolutions) > 1 :
      
          button.config(text=text,  bg='red')                                                           # The red color signals that the spinssystem is used in more than 1 place in the sequence
          
        else :
          
          button.config(text=text,  bg='grey83')  
        
      else :
        
        text = '-' #str(res.getSeqCode()) + ' ' + res.getCcpCode() + ': -'
      
        button.config(text=text)
        
  def updateButtonHighLights(self):
    
    self.setAllButtonsToGrey()
    
    if self.selectedResidueA :
      
      buttons = [self.sequenceButtons.buttons[0], self.sequenceButtons.buttons[2], self.sequenceButtons.buttons[4], self.sequenceButtons.buttons[6],self.sequenceButtons.buttons[8]]
      buttons[self.selectedResidueA - 1].config(bg='yellow')
      
    
    elif self.selectedResidueB :
      
      buttons = [self.sequenceButtonsB.buttons[0], self.sequenceButtonsB.buttons[2], self.sequenceButtonsB.buttons[4], self.sequenceButtonsB.buttons[6],self.sequenceButtonsB.buttons[8]]
      buttons[self.selectedResidueB - 1].config(bg='yellow')  
    
    elif self.selectedLinkA :
      
      buttons = [self.sequenceButtons.buttons[1], self.sequenceButtons.buttons[3], self.sequenceButtons.buttons[5], self.sequenceButtons.buttons[7]]
      buttons[self.selectedLinkA - 1].config(bg='yellow')

    elif self.selectedLinkB :
      
      buttons = [self.sequenceButtonsB.buttons[1], self.sequenceButtonsB.buttons[3], self.sequenceButtonsB.buttons[5], self.sequenceButtonsB.buttons[7]]
      buttons[self.selectedLinkB - 1].config(bg='yellow')
      
  def setAllButtonsToGrey(self) :
    
    for button in self.sequenceButtons.buttons + self.sequenceButtonsB.buttons :
      
      button.config(bg='grey83')
    
  def sortDisplayResultsTable(self) :                                                       # Remove: not used any longer.
    
    for table in self.displayResultsTables :
      
      #table.sortDown = False
      table.sortLine(1, noUpdate=True)
      
  def getCurrentlyDisplayedResidues(self) :
    
    resNumber = self.resultsResidueNumber    

    residues = self.connector.results.myChain.residues[resNumber - 3 : resNumber + 2]
  
    return residues 

  def pickColorByPercentage(self, percentage):

    percentage = float(percentage)
    
    if percentage < 1 :
      
      return 'grey83'
    
    if percentage > 80 :

      red = 0
      green = int(percentage/100.0 * 255.0)
      blue = int(255.0 - (percentage/100.0 * 255.0))
      

    elif percentage > 50 :

      green = 0 
      blue = int(percentage/80.0 * 255.0)
      red = int(255.0 - (percentage/80.0 * 255.0))

    else :

      red = 255
      green = 0
      blue = 0

    red = self.rgbToHex(red)
    green = self.rgbToHex(green)
    blue = self.rgbToHex(blue)

    color = '#' + red + green + blue
    return color
    
  def rgbToHex(self,rgb):

    div = str(rgb/16)
    rem = str(rgb%16)

    for i in [['10','A'],['11','B'],['12','C'],['13','D'],['14','E'],['15','F']] :

      if div == i[0] :
        
        div = i[1]

      if rem == i[0] :
        
        rem = i[1]

    h = div + rem

    return h
   
  def toggleTab(self, index):
  
    pass
     
  def administerNotifiers(self, notifyFunc):

    notifyFunc(self.spectrumSelectionTab.update, 'ccp.nmr.Nmr.Experiment', 'setName')
    notifyFunc(self.spectrumSelectionTab.update, 'ccp.nmr.Nmr.DataSource', 'setName')
    notifyFunc(self.spectrumSelectionTab.update, 'ccp.nmr.Nmr.DataSource', 'delete')
    notifyFunc(self.spectrumSelectionTab.update, 'ccp.nmr.Nmr.DataSource', '__init__')
    
    notifyFunc(self.annealingSettingsTab.updateChains, 'ccp.molecule.MolSystem.Chain', 'delete')
    notifyFunc(self.annealingSettingsTab.updateChains, 'ccp.molecule.MolSystem.Chain', '__init__')
    
    
     
  def updateAfter(self, *opt):

    if self.waiting:
      return
    else:
      self.waiting = True
      self.after_idle(self.update)
    
  def update(self):
    
    self.toggleTab(self.tabbedFrame.selected) 
    self.waiting = False
      
  def destroy(self):                                                      
  
    self.administerNotifiers(self.unregisterNotify)
    
    BasePopup.destroy(self)
  
  def noAction(self):
      
    pass
    
  def updateInfoText(self,string):
    
    self.infoLabel.set(string)
    
    self.infoLabel.update()
    