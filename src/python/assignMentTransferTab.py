
from memops.gui.CheckButton import CheckButton
from memops.gui.Label import Label
from memops.gui.RadioButtons import RadioButtons
from memops.gui.Frame import Frame
from memops.gui.LabelFrame import LabelFrame
from memops.gui.FloatEntry import FloatEntry
from memops.gui.IntEntry import IntEntry
from memops.gui.PulldownList import PulldownList
from memops.gui.ButtonList import ButtonList

from ccpnmr.analysis.core.AssignmentBasic import assignResToDim, assignSpinSystemResidue, mergeSpinSystems

from Tkinter import VERTICAL

from assignmentFunctions import assignSpinSystemstoResidues

class AssignMentTransferTab(object) :

  def __init__(self, parent, frame) :
    
    self.guiParent = parent
  
    self.frame = frame
    
    self.minScore = 80.0
    
    self.dataModel = None
    
    self.spectrum = None
    
    self.selectedSolution = 1
    
    self.body()
    
    self.resonanceToDimension = True
    self.spinSystemToResidue = True
    self.intra = True
    self.sequential = True
    self.noDiagonal = True
    self.allSpectra = True
    
    self.spinSystemType = 0
    
  def body(self) :
    
    #self.frame.expandColumn(0)
    self.frame.expandGrid(8,0)
    self.frame.expandGrid(8,1)
    
    typeOfAssignmentFrame = LabelFrame(self.frame, text='type of assignment')
    typeOfAssignmentFrame.grid(row=0, column=0, sticky='nesw')
    #typeOfAssignmentFrame.expandGrid(0,5)
    
    peakSelectionFrame = LabelFrame(self.frame, text='which peaks to assign')
    peakSelectionFrame.grid(row=0, column=1, sticky='nesw',rowspan=2)
    
    spinSystemSelectionFrame = LabelFrame(self.frame, text='Which spin-systems to use')
    spinSystemSelectionFrame.grid(row=2, column=0, sticky='nesw')
    
    tipText='What to do when a residue has already a spin system assigned to it.'
    assignedResidueFrame = LabelFrame(self.frame, text='if residue already has spin-system', tipText=tipText)
    assignedResidueFrame.grid(row=2, column=1, sticky='nesw')
    
    spectrumSelectionFrame = LabelFrame(self.frame, text='spectra')
    spectrumSelectionFrame.grid(row=1, column=0, sticky='nesw')
    
    row = 0
    
    label = Label(typeOfAssignmentFrame, text='Resonances to Peak Dimensions', grid=(row,0)) 
    self.peaksCheckButton = CheckButton(typeOfAssignmentFrame, selected=True, grid=(row,1))
    
    row += 1
    
    label = Label(typeOfAssignmentFrame, text='SpinSystems to Residues', grid=(row,0)) 
    self.residuesCheckButton = CheckButton(typeOfAssignmentFrame, selected=True, grid=(row,1))
    
    row = 0
    
    label = Label(peakSelectionFrame, text='Intra-Residual', grid=(row,0)) 
    self.intraCheckButton = CheckButton(peakSelectionFrame, selected=True, grid=(row,1))
    
    row += 1
    
    label = Label(peakSelectionFrame, text='Sequential', grid=(row,0)) 
    self.sequentialCheckButton = CheckButton(peakSelectionFrame, selected=True, grid=(row,1))
    
    row += 1
    
    label = Label(peakSelectionFrame, text='Do not assign diagonal peaks', grid=(row,0)) 
    self.noDiagonalCheckButton = CheckButton(peakSelectionFrame, selected=True, grid=(row,1))
    
    entries = ['Only assigned spin systems','All that have a score of at least: ','User Defined', 'Solution number:']
    tipTexts = ['Only assign resonances of spin systems that already have a sequential assignment for the assignment of peak dimensions. Spin system to residue assignment is not relevant in this case.',
                'Assign all spin systems that have a score of at least a given percentage. 50% or lower is not possible, because than spin systems might have to be assigned to more than 1 residue, which is impossible.',
                "As defined in the lower row of buttons in the 'results' tab.",
                'One of the single solutions of the annealing.']
    self.spinSystemTypeSelect = RadioButtons(spinSystemSelectionFrame, entries=entries, grid=(0,0),
                                     select_callback=None, direction = VERTICAL, 
                                     gridSpan=(4,1), tipTexts=tipTexts)
    

    
    tipText = 'The minimal amount of colabelling the different nuclei should have in order to still give rise to a peak.'
    self.minScoreEntry = FloatEntry(spinSystemSelectionFrame, grid=(1,1), width=7, text=str(self.minScore),
                                      returnCallback=self.changeMinScore,
                                      tipText=tipText)
    self.minScoreEntry.bind('<Leave>', self.changeMinScore, '+')
    
    self.solutionNumberEntry = IntEntry(spinSystemSelectionFrame, grid=(3,1), width=7, text=1,
                                      returnCallback=self.solutionUpdate,
                                      tipText=tipText)
    self.solutionNumberEntry.bind('<Leave>', self.solutionUpdate, '+')
    
    #self.solutionPullDown = PulldownList(spinSystemSelectionFrame, None, grid=(3,1), sticky='w')
    
    entries = ['all spectra','only:']
    tipTexts = ['Assign peaks in all the spectra that where selected before the annealing ran.',
                'Only assign peaks in one particular spectrum. You can of course repeat this multiple times for different spectra.']
    self.spectrumSelect = RadioButtons(spectrumSelectionFrame, entries=entries, grid=(0,0),
                                     select_callback=None, direction = VERTICAL, 
                                     gridSpan=(2,1), tipTexts=tipTexts)
    
    self.spectraPullDown = PulldownList(spectrumSelectionFrame, self.changeSpectrum, grid=(1,1), sticky='w')
    
    
    
    
    entries = ['skip this residue','de-assign old spin system from residue','assign, but never merge','warn to merge']
    tipTexts = ["Don't assign the new spin system to the residue. The residue is not skipped when the old spin system does not contain any resonances",
                "De-assign old spin system from residue, unless the old spin system is a spin system without any resonances.",
                "Don't merge any spin systems, merging can be performed later if nescesary in the Resonance --> SpinSystems window.",
                "Ask to merge individually for each spin system, this might result in clicking on a lot of popups."]
    self.assignedResidueStrategySelect = RadioButtons(assignedResidueFrame, entries=entries, grid=(0,0),
                                     select_callback=None, direction = VERTICAL, 
                                     gridSpan=(2,1), tipTexts=tipTexts)
    
    
    
    
    
    texts    = ['Transfer Assignments']
    commands = [self.transferAssignments]
    self.transferButton = ButtonList(self.frame,commands=commands, texts=texts)
    self.transferButton.grid(row=5, column=0, sticky='nsew', columnspan=2) 
    
  def update(self) :
    
    self.dataModel = self.guiParent.connector.results
    self.updateSpectra()
  
  def setDataModel(self, dataModel) :
    
    self.dataModel = dataModel
    
    self.update()
    
  def updateSpectra(self, *opt):
    
    if not self.dataModel :
      
      return
    
    spectrum = self.spectrum 
    
    spectra = self.dataModel.getSpectra()
    
    if spectra :
      
      names = [spectrum.name for spectrum in spectra]
      index = 0
      
      if self.spectrum not in spectra :
        
        self.spectrum = spectra[0]
        
      else :
        
        index = spectra.index(self.spectrum)
        
    self.spectraPullDown.setup(names,spectra, index)
      
  def changeSpectrum(self, spectrum):
    
    self.spectrum = spectrum
    
  def solutionUpdate(self, event=None, value=None):

    if not self.dataModel :
      
      return
  
    Nsolutions = len(self.dataModel.chain.residues[0].solutions)

    if value is None :
      
      value = self.solutionNumberEntry.get()
      
    if value == self.selectedSolution:
      return 
    else :
      self.selectedSolution = value
    if value < 1:
      self.solutionNumberEntry.set(1)
      self.selectedSolution = 1
    elif value > Nsolutions:
      self.selectedSolution = Nsolutions
      self.solutionNumberEntry.set(self.selectedSolution)
    else :
      self.solutionNumberEntry.set(self.selectedSolution)
  
  def fetchOptions(self) :
    
    self.resonanceToDimension = self.peaksCheckButton.get()
    self.spinSystemToResidue = self.residuesCheckButton.get()
    self.intra = self.intraCheckButton.get()
    self.sequential = self.sequentialCheckButton.get()
    self.noDiagonal = self.noDiagonalCheckButton.get()
    self.spinSystemType = self.spinSystemTypeSelect.getIndex()
    self.strategy = ['skip','remove','noMerge',None][self.assignedResidueStrategySelect.getIndex()]
    self.allSpectra = [True,False][self.spectrumSelect.getIndex()]
      
  def changeMinScore(self, event=None) :
    
    newMinScore = self.minScoreEntry.get()
    
    if self.minScore != newMinScore :
      
      if newMinScore <= 50.0 :
        
        self.minScore = 51.0
        self.minScoreEntry.set(51.0)
        
      elif newMinScore > 100.0 :
        
        self.minScore = 100.0
        self.minScoreEntry.set(100.0)
      
      else :
        
        self.minScore = newMinScore
        
  def transferAssignments(self) :
    
    self.fetchOptions()
    
    if not self.dataModel or (not self.resonanceToDimension and not self.spinSystemToResidue) :
      
      return
    
    strategy = self.strategy
    
    lookupSpinSystem = [self.getAssignedSpinSystem, self.getBestScoringSpinSystem, self.getUserDefinedSpinSystem, self.getSelectedSolutionSpinSystem][self.spinSystemType]
    
    residues = self.dataModel.chain.residues
    
    spinSystemSequence = [lookupSpinSystem(res) for res in residues]
    
    #print spinSystemSequence
    
    ccpnSpinSystems = []
    ccpnResidues = []
    
    if self.spinSystemToResidue and not self.spinSystemType == 0 :          # if self.spinSystemType == 0 it means that it for sure already assigned like this
      
      for spinSys,res in zip(spinSystemSequence, residues) :
        
        if spinSys and res:
          
          ccpnSpinSystems.append(spinSys.getCcpnResonanceGroup())
          ccpnResidues.append(res.getCcpnResidue())
      
      assignSpinSystemstoResidues(ccpnSpinSystems, ccpnResidues, strategy=strategy, guiParent=self.guiParent)
      
      
          
    if self.resonanceToDimension :
      
      allSpectra = self.allSpectra
      
      if self.intra :
        
        for residue, spinSystem in zip(residues, spinSystemSequence) :
          
          if not spinSystem :
            
            continue
            
          intraLink = residue.getIntraLink(spinSystem)
          
          for pl in intraLink.getPeakLinks() :
            
            peak = pl.getPeak()
            
            if not allSpectra and peak.getSpectrum() is not self.spectrum :
              
              continue
            
            if not peak :
              
              continue
              
            resonances = pl.getResonances()
            
            if self.noDiagonal and len(set(resonances)) < len(resonances) :
              
              continue
            
            for resonance, dimension in zip(resonances, peak.getDimensions()) :
              
              ccpnResonance = resonance.getCcpnResonance()
              ccpnDimension = dimension.getCcpnDimension()
              
              assignResToDim(ccpnDimension, ccpnResonance)
        
      if self.sequential :
        
        for residue, spinSystemA, spinSystemB in zip(residues, spinSystemSequence, spinSystemSequence[1:]) :
          
          if not spinSystemA or not  spinSystemB:
            
            continue
            
          link = residue.getLink(spinSystemA, spinSystemB)
          
          for pl in link.getPeakLinks() :
            
            peak = pl.getPeak()
            
            if not allSpectra and peak.getSpectrum() is not self.spectrum :
              
              continue
            
            if not peak :
              
              continue
              
            resonances = pl.getResonances()
            
            if self.noDiagonal and len(set(resonances)) < len(resonances) :
              
              continue
            
            for resonance, dimension in zip(resonances, peak.getDimensions()) :
              
              ccpnResonance = resonance.getCcpnResonance()
              ccpnDimension = dimension.getCcpnDimension()
              
              assignResToDim(ccpnDimension, ccpnResonance)
    
    self.guiParent.resultsTab.update()
              
  def getAssignedSpinSystem(self, res) :
    
    ccpCode = res.ccpCode
    seqCode = res.getSeqCode()
    spinSystems= self.dataModel.getSpinSystems()[ccpCode]
    
    ccpnResidue = res.getCcpnResidue()
    if ccpnResidue:
      assignedResonanceGroups = ccpnResidue.getResonanceGroups()
      if len(assignedResonanceGroups) > 1:
        print 'There is more than one spin system assigned to residue %s, did not know which one to use to assign peaks. Therefor this residue is skipped.' %(seqCode)
        return
      
      assignedResonanceGroup = ccpnResidue.findFirstResonanceGroup()
      
      if assignedResonanceGroup:
        
        for spinSystem in spinSystems :
          
          if spinSystem.getCcpnResonanceGroup is assignedResonanceGroup :    
            if not self.skipResidue(res,spinSystem) :         # Just checking to make sure, analysis project could have changed
              
              return spinSystem
                   
   
  def getBestScoringSpinSystem(self, res) :
    
    ccpCode = res.ccpCode
    seqCode = res.getSeqCode()
    spinSystems= self.dataModel.getSpinSystems()[ccpCode]

    solutions = res.solutions
    
    weigth = 1.0/len(solutions)
    
    score, bestSpinSystem = max([(solutions.count(solution)*weigth*100.0, solution) for solution in solutions])
    
    if score >= self.minScore and not bestSpinSystem.getIsJoker() and not self.skipResidue(res,bestSpinSystem) :
      
      return bestSpinSystem
    
    return None
    
  def getUserDefinedSpinSystem(self, res) :
    
    userDefinedSpinSystem = res.userDefinedSolution
    
    if userDefinedSpinSystem and not userDefinedSpinSystem.getIsJoker() and not self.skipResidue(res,userDefinedSpinSystem) :
      
      return userDefinedSpinSystem
    
    return None
  
  def getSelectedSolutionSpinSystem(self, res) :
    
    solutions = res.solutions
    
    spinSystem = solutions[self.selectedSolution-1]
    
    if not spinSystem.getIsJoker() and not self.skipResidue(res,spinSystem) :
      
      return spinSystem
    
    return None
  
  def skipResidue(self,residue,spinSystem) :
    
    if self.strategy == 0 :
      
      assignedSpinSystems = residue.getCcpnResidue().getResonanceGroups()
      
      if assignedSpinSystems and spinSystem not in assignedSpinSystems :
        
        return True
      
    return False
  