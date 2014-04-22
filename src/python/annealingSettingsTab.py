
from memops.gui.Button import Button
from memops.gui.CheckButton import CheckButton
from memops.gui.FloatEntry import FloatEntry
from memops.gui.IntEntry import IntEntry
from memops.gui.Label import Label
from memops.gui.PulldownList import PulldownList
from memops.gui.ScrolledGraph import ScrolledGraph
from memops.gui.Entry import Entry
from ccpnmr.analysis.core.AssignmentBasic import getShiftLists

class AnnealingSettingsTab(object) :

  def __init__(self, parent, frame) :
    
    self.guiParent = parent
  
    self.frame = frame
    
    self.project = parent.project
    
    self.nmrProject = parent.nmrProject
    
    self.minIsoFrac = 0.1
    self.leavePeaksOutFraction = 0.0
    self.minTypeScore = 1.0
    self.chain         = None
    self.amountOfRepeats = 10
    self.amountOfSteps = 10000
    #self.acceptanceConstantList = [0, 0.1, 0.2,0.4, 0.8,1.0, 1.1, 1.2, 1.4, 1.6,2.0, 2.4,2.8, 3.2,4.0,4.5,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0]
    self.acceptanceConstantList = [0.0, 0.01, 0.015, 0.022, 0.033, 0.050, 0.075, 0.113, 0.170, 0.256, 0.384, 0.576, 0.864, 1.297, 1.946, 2.919, 4.378, 6.568, 9.852, 14.77, 22.16, 33.25]
    self.energyDataSets = [[]]
    
    self.body()

  def body(self) :
  
    frame = self.frame

    #frame.expandGrid(13,0)
    frame.expandGrid(13,1)
    row = 0

    text = 'Calculate Assignment Suggestions'
    command = self.runCalculations
    self.startButton = Button(frame,command=command, text=text)
    self.startButton.grid(row=row, column=0, sticky='nsew', columnspan=2)
    
    row += 1

    label = Label(frame, text='Amount of runs: ', grid=(row,0))      
    tipText = 'The amount of times the whole optimization procedure is performed, each result is safed'
    self.repeatEntry = IntEntry(frame, grid=(row,1), width=7, text=10,
                                      returnCallback=self.updateRepeatEntry,
                                      tipText=tipText, sticky= 'nsew')
    self.repeatEntry.bind('<Leave>', self.updateRepeatEntry, '+')
    
    row += 1
    
    label = Label(frame, text='Temperature regime: ', grid=(row,0))
    tipText = 'This list of numbers govern the temperature steps during the annealing, every number represents 1/(kb*t), where kb is the Boltzmann constant and t the temperature of one step.'
    self.tempEntry = Entry(frame, text=map(str, self.acceptanceConstantList), width=64,
                           grid=(row,1), isArray=True, returnCallback=self.updateAcceptanceConstantList,
                           tipText=tipText, sticky= 'nsew')
    
    row += 1    

    label = Label(frame, text='Amount of attempts per temperature:', grid=(row,0))      
    tipText = 'The amount of attempts to switch the position of two spinsystems in the sequence are performed for each temperature point'
    self.NAStepEntry = IntEntry(frame, grid=(row,1), width=7, text=10000,
                                      returnCallback=self.updateStepEntry,
                                      tipText=tipText, sticky= 'nsew')
    self.NAStepEntry.bind('<Leave>', self.updateStepEntry, '+')
    
    row += 1
        
    label = Label(frame, text='Fraction of peaks to leave out:', grid=(row,0))      
    tipText = 'In each run a fraction of the peaks can be left out of the optimization, thereby increasing the variability in the outcome and reducing false negatives. In each run this will be different randomly chosen sub-set of all peaks. 0.1 (10%) can be a good value.'
    self.leaveOutPeaksEntry = FloatEntry(frame, grid=(row,1), width=7, text=0.0,
                                      returnCallback=self.updateLeavePeaksOutEntry,
                                      tipText=tipText, sticky= 'nsew')
    self.leaveOutPeaksEntry.bind('<Leave>', self.updateLeavePeaksOutEntry , '+')
    
    row += 1
    
    label = Label(frame, text='Minmal amino acid typing score:', grid=(row,0))      
    tipText = 'If automatic amino acid typing is selected, a cut-off value has to set. Every amino acid type that scores higher than the cut-off is taken as a possible type. This is the same score as can be found under resonance --> spin systems --> predict type. Value should be between 0 and 100'
    self.minTypeScoreEntry = FloatEntry(frame, grid=(row,1), width=7, text=1.0,
                                      returnCallback=self.updateMinTypeScoreEntry,
                                      tipText=tipText, sticky= 'nsew')
    self.minTypeScoreEntry.bind('<Leave>', self.updateMinTypeScoreEntry, '+')
    

    row += 1

    label = Label(frame, text='Minimal colabelling fraction:', grid=(row,0))      
    tipText = 'The minimal amount of colabelling the different nuclei should have in order to still give rise to a peak.'
    self.minLabelEntry = FloatEntry(frame, grid=(row,1), width=7, text=0.1,
                                      returnCallback=self.updateMinLabelEntry,
                                      tipText=tipText, sticky= 'nsew')
    self.minLabelEntry.bind('<Leave>', self.updateMinLabelEntry, '+')
    
    row += 1
    
    label= Label(frame, text = 'Use sequential assignments:', grid=(row,0))
    tipText = 'When this option is select the present sequential assignments will be kept in place'
    self.useAssignmentsCheck = CheckButton(frame, selected=True, tipText=tipText, grid=(row,1))
    
    row += 1
    
    label= Label(frame, text = 'Use tentative assignments:', grid=(row,0))
    tipText = 'If a spin system has tentative assignments this can be used to narrow down the amount of possible sequential assignments.'
    self.useTentativeCheck = CheckButton(frame, selected=True, tipText=tipText, grid=(row,1))

    row += 1
    
    label = Label(frame, text='Use amino acid types:', grid=(row,0))      
    tipText = 'Use amino acid types of the spin systems. If this option is not checked the spin systems are re-typed, only resonance names and frequencies are used'
    self.useTypeCheck = CheckButton(frame, selected=True, tipText=tipText, grid=(row,1))
     
    row += 1
    
    label = Label(frame, text='Include untyped spin systems:', grid=(row,0))      
    tipText = 'Also include spin system that have no type information. Amino acid typing will be done on the fly.'
    self.useAlsoUntypedSpinSystemsCheck = CheckButton(frame, selected=True, tipText=tipText, grid=(row,1))

    row += 1

    label = Label(frame, text='Use dimensional assignments:', grid=(row,0))      
    tipText = 'If one or more dimensions of a peak is already assigned, assume that this assignment is the only option. If not the check the program will consider all possibilities for the assignment of the dimension.'
    self.useDimenionalAssignmentsCheck = CheckButton(frame, selected=True, tipText=tipText, grid=(row,1))
    
    row += 1

    self.shiftListLabel    = Label(frame, text ='Shift List:', grid=(row,0), sticky='w')
    self.shiftListPulldown = PulldownList(frame, self.changeShiftList, grid=(row,1), sticky='w')
    
    row += 1
    
    label = Label(frame, text='Chain:', grid=(row,0))
    self.molPulldown = PulldownList(frame, callback=self.changeMolecule, grid=(row,1))
    self.updateChains()  
    
    self.updateShiftLists()
    
    row += 1
    
    self.energyPlot = ScrolledGraph(frame,symbolSize=2, width=600,
                                       height=200, title='Annealing',
                                       xLabel='temperature step', yLabel='energy')
    self.energyPlot.grid(row=row, column=0, columnspan=2, sticky='nsew')
    
  def runCalculations(self):
    
    self.startButton.disable()
    self.disableIllegalButtonsAfterPrecalculations()
    self.guiParent.connector.runAllCalculations()
    self.startButton.configure(text='More runs')
    self.startButton.enable()
    
  def disableIllegalButtonsAfterPrecalculations(self):
    
    illegalButtons = [self.minTypeScoreEntry, self.minLabelEntry,
                      self.useAlsoUntypedSpinSystemsCheck, self.useAssignmentsCheck,
                      self.useTypeCheck, self.useDimenionalAssignmentsCheck,
                      self.useTentativeCheck]
    
    for illegalButton in illegalButtons :
      illegalButton.configure(state='disabled')
      
    self.molPulldown.disable()
    self.shiftListPulldown.disable()
      
     
  def getChainName(self, chain):
  
    return '%s:%s (%s)' % (chain.molSystem.code,chain.code,chain.molecule.molType)
    
  def getChains(self):
  
    chains = []
    if self.project:
      for molSystem in self.project.sortedMolSystems():
        for chain in molSystem.sortedChains():
          if chain.residues:
            chains.append(chain)
   
    return chains
    
  def updateChains(self, *opt):
  
    index  = 0
    texts  = []
    chains = self.getChains()
    chain  = self.chain

    if chains:
      if chain not in chains:
        chain = chains[0]

      texts = [self.getChainName(c) for c in chains]
      index = chains.index(chain)
        
    else:
      chain = None  
      
    self.molPulldown.setup(texts,chains,index)
  
    if chain is not self.chain:
      self.chain = chain
       
  def changeMolecule(self, chain):
  
    if chain is not self.chain:
      self.chain = chain
      
  def changeShiftList(self, shiftList):
    
    if self.shiftList is not shiftList:
      self.shiftList = shiftList

  def updateShiftLists(self):
  
    index = 0
    shiftLists = getShiftLists(self.nmrProject) + [None,]
    shiftListNames = self.getShiftListNames(shiftLists[:-1]) + ['None',]
    
    if shiftListNames:
      self.shiftList = shiftLists[0]
      index = 0
      
    else:
      shiftList = None
    
    self.shiftListPulldown.setup(shiftListNames,shiftLists,index)

  def getShiftListNames(self, shiftLists):
    
    shiftListNames = []
    for shiftList in shiftLists:
      if not shiftList.name:
        shiftList.name = "ShiftList "+ str(shiftList.serial)
      shiftListNames.append(shiftList.name)

    return shiftListNames
  
  def updateStepEntry(self,  event =None):

    value = self.NAStepEntry.get()
    if value == self.amountOfSteps:
      return 
    if value < 1:
      self.NAStepEntry.set(1)
      self.amountOfSteps = 1
    else :
      self.amountOfSteps = value
      self.NAStepEntry.set(value)
    
  def updateRepeatEntry(self,  event =None):

    value = self.repeatEntry.get()

    if value == self.amountOfRepeats:
      return 
    if value < 1:
      self.repeatEntry.set(1)
      self.amountOfRepeats = 1
    else :
      self.amountOfRepeats = value
      self.repeatEntry.set(value)
      
  def updateMinTypeScoreEntry(self, event = None) :
    
    value = self.minTypeScoreEntry.get()
    
    if value == self.minTypeScore :
      return
    if value < 0:
      self.minTypeScoreEntry.set(0.0)
      self.minTypeScore = 0.0
    elif value > 100 :
      self.minTypeScoreEntry.set(100.0)
      self.minTypeScore = 100.0 
    else :
      self.minTypeScoreEntry.set(value)
      self.minTypeScore = value
      
    
      
  def updateMinLabelEntry(self, event=None) :
    
    value = self.minLabelEntry.get()
    
    if value == self.minIsoFrac :
      return
    if value < 0 :
      self.minIsoFrac = 0.0
      self.minLabelEntry.set(0.0)
    elif value > 1 :
      self.minIsoFrac = 1.0
      self.minLabelEntry.set(1.0)
    else :
      self.minIsoFrac = value
      self.minLabelEntry.set(value)
      
  def updateLeavePeaksOutEntry(self, event=None) :
    
    value = self.leaveOutPeaksEntry.get()
    
    if value == self.leavePeaksOutFraction :
      return
    if value < 0 :
      self.leavePeaksOutFraction = 0.0
      self.leaveOutPeaksEntry.set(0.0)
    elif value > 1 :
      self.leavePeaksOutFraction = 1.0
      self.leaveOutPeaksEntry.set(1.0)
    else :
      self.leavePeaksOutFraction = value
      self.leaveOutPeaksEntry.set(value)  
    
       
  def updateAcceptanceConstantList(self,event=None) :
    
    acList = self.tempEntry.get()
    
    newList = []
    
    for constant in acList :
      
      try :
        
        number = float(constant)
        
        newList.append(number)
        
      except ValueError :

        string = constant + ' in temperature constants is not a number.'
        
        showWarning('Not A Number', string, parent=self.guiParent)
        
        return False
        
  
    self.acceptanceConstantList = newList
    
    return True
  
  def addEnergyPoint(self,energy,time) :
    
    point = (time,energy*-1)
    
    if len(self.energyDataSets[-1]) / (len(self.acceptanceConstantList) + 1) :                # This means one run has finished
      
      self.energyDataSets.append([point])
      
    else :
      
      self.energyDataSets[-1].append(point)
      
    colors = ['#993366','#000000','#FF0099','#33FF00','#003300','#999999','#FF6633','#000099','#33CCFF','#FFCC00']
    
    Ncolors = len(colors)
    
    NdataSets = len(self.energyDataSets)
    
    colorList = (NdataSets/Ncolors)*colors + colors[:NdataSets%Ncolors]
    
    self.energyPlot.update(dataSets=self.energyDataSets, dataColors=colorList)
    
    # Forcing the graph to draw, eventhough calculations
    # are still running. Only do this with high numbers of
    # steps, otherwise drawing takes longer than annealling.
    if self.amountOfSteps >= 100000 :
      
      self.energyPlot.draw()
  