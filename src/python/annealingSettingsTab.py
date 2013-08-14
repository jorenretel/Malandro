
from memops.gui.ButtonList import ButtonList
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
    self.chain         = None
    self.amountOfRepeats = 10
    self.amountOfSteps = 10000
    self.acceptanceConstantList = [0, 0.1, 0.2,0.4, 0.8,1.0, 1.1, 1.2, 1.4, 1.6,2.0, 2.4,2.8, 3.2,4.0,4.5,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0]
    self.energyDataSets = [[]]
    
    self.body()

  def body(self) :
  
    frame = self.frame

    frame.expandGrid(11,0)
    frame.expandGrid(11,1)
    row = 0

    texts    = ['Calculate Assignment Suggestions']
    commands = [self.guiParent.connector.runAllCalculations]
    self.startNAButton = ButtonList(frame,commands=commands, texts=texts)
    self.startNAButton.grid(row=row, column=0, sticky='nsew', columnspan=2)
    
    row += 1
    

    label = Label(frame, text='Minimal colabelling fraction:', grid=(row,0))      
    tipText = 'The minimal amount of colabelling the different nuclei should have in order to still give rise to a peak.'
    self.minLabelEntry = FloatEntry(frame, grid=(row,1), width=7, text=0.1,
                                      returnCallback=self.updateMinLabelEntry,
                                      tipText=tipText)
    self.minLabelEntry.bind('<Leave>', self.updateMinLabelEntry, '+')
    
    row += 1

    label = Label(frame, text='Amount of attempts per temperature point in the annealling:', grid=(row,0))      
    tipText = 'The amount of attempts to switch the position of two spinsystems in the sequence are performed for each temperature point'
    self.NAStepEntry = IntEntry(frame, grid=(row,1), width=7, text=10000,
                                      returnCallback=self.updateStepEntry,
                                      tipText=tipText)
    self.NAStepEntry.bind('<Leave>', self.updateStepEntry, '+')
    
    row += 1
    
    label = Label(frame, text=' Perform N runs of the optimsation, N', grid=(row,0))      
    tipText = 'The amount of times the whole optimization procedure is performed, each result is safed'
    self.repeatEntry = IntEntry(frame, grid=(row,1), width=7, text=10,
                                      returnCallback=self.updateRepeatEntry,
                                      tipText=tipText)
    self.repeatEntry.bind('<Leave>', self.updateRepeatEntry, '+')
    
    row += 1
    
    label = Label(frame, text='Temperature constants: ', grid=(row,0))   
    self.tempEntry = Entry(frame, text=map(str, self.acceptanceConstantList), width=64, grid=(row,1), isArray=True, returnCallback=self.updateAcceptanceConstantList)
    
    row += 1
    
    label = Label(frame, text='Re-type all typed spin systems:', grid=(row,0))      
    tipText = 'The algorithm does not take into account any previously made assignment and re-types all spin systems, only resonance names and frequencies are used'

    self.reTypeSpinSystemsCheck = CheckButton(frame, selected=False, grid=(row,1))
     
    row += 1
    
    label = Label(frame, text='Type untyped spin systems on the fly:', grid=(row,0))      
    tipText = 'Spin system typing can be carried out so also untyped spin systems can be used. Only amino acid types that have a reasonable score (i.e. higher than the random chance) will be considered'

    self.typeSpinSystemsCheck = CheckButton(frame, selected=False, grid=(row,1))

    row += 1

    label = Label(frame, text='If a peak dimension has resonances assigned, only consider those:', grid=(row,0))      
    tipText = 'If one or more dimensions of a peak is already assigned, assume that this assignment is the only option. If not the program will assume the peak can have other assignment options as well.'
    self.useDimenionalAssignmentsCheck = CheckButton(frame, selected=True, grid=(row,1))
    
    row += 1

    label= Label(frame, text = 'Keep sequentially assigned spin systems in place:', grid=(row,0))
    self.useAssignmentsCheck = CheckButton(frame, selected=True, grid=(row,1))
    
    row += 1
    
    label= Label(frame, text = 'If a spin system has tentative residue assignments, only consider those:', grid=(row,0))
    self.useTentativeCheck = CheckButton(frame, selected=True, grid=(row,1))

    row += 1

    self.shiftListLabel    = Label(frame, text ='Shift List:', grid=(row,0), sticky='w')
    self.shiftListPulldown = PulldownList(frame, self.changeShiftList, grid=(row,1), sticky='w')
    
    row += 1
    
    label = Label(frame, text='Chain:', grid=(row,0))
    self.molPulldown = PulldownList(frame, callback=self.changeMolecule, grid=(row,1))
    self.updateChains()  
    
    self.updateShiftLists()
    
    row += 1
    
    self.energyPlot = ScrolledGraph(frame,symbolSize=2, width=500,
                                       height=300, title='Annealing',
                                       xLabel='time', yLabel='energy')
    self.energyPlot.grid(row=row, column=0, columnspan=2, sticky='nsew')
     
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
      
  def updateMinLabelEntry(self, event=None) :
    
    value = self.minLabelEntry.get()
    
    if value == self.minIsoFrac :
      return
    if value < 0 :
      self.minIsoFrac = 0.0
      self.minLabelEntry.set(0.0)
    else :
      self.minIsoFrac = value
      self.minLabelEntry.set(value)
       
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
    
    point = (time,energy)
    
    if len(self.energyDataSets[-1]) / len(self.acceptanceConstantList) :                # This means one run has finished
      
      self.energyDataSets.append([point])
      
    else :
      
      self.energyDataSets[-1].append(point)
      
    colors = ['#993366','#000000','#FF0099','#33FF00','#003300','#999999','#FF6633','#000099','#33CCFF','#FFCC00']
    
    Ncolors = len(colors)
    
    NdataSets = len(self.energyDataSets)
    
    colorList = (NdataSets/Ncolors)*colors + colors[:NdataSets%Ncolors]
    
    self.energyPlot.update(dataSets=self.energyDataSets, dataColors=colorList)  
  