
from memops.gui.PulldownList import PulldownList
from memops.gui.ScrolledMatrix  import ScrolledMatrix

class SpectrumSelectionTab(object) :

  def __init__(self, parent, frame) :
    
    self.guiParent = parent
  
    self.frame = frame
    
    self.project = parent.project
    
    self.nmrProject = parent.nmrProject
    
    self.specConfigList = []
    
    self.body()
    self.updateSpecSelection()
    #self.setupSpectrumSettings()
    self.updateAutoMatrix()


  def body(self) :
    
    frame = self.frame

    frame.grid_rowconfigure(2, weight=1)
    frame.grid_columnconfigure(0,  weight=1)
    

    headingList = ['#','Spectrum', 'Peak List','use?','labelling scheme']

    tipTexts = [ 'Row number','spectrum name', 'which peak list to use',  'use this spectrum?', 'Which labelling scheme belongs to this spectrum?']


    self.autoLabellingPulldown = PulldownList(self.guiParent, self.setAutoLabellingScheme)
    self.autoPeakListPulldown = PulldownList(self.guiParent, self.setAutoPeakList)


    editWidgets = [None, None,self.autoPeakListPulldown,  None, self.autoLabellingPulldown]

    editGetCallbacks = [None, None, self.getAutoPeakLists,  self.changeUse, self.getAutoLabellingSchemes]

    editSetCallbacks = [None, None, None, None, None]

    self.displayTable = ScrolledMatrix(frame,headingList=headingList,
                                       callback=self.selectAutoSpec,
                                       editWidgets=editWidgets, multiSelect=False,
                                       editGetCallbacks=editGetCallbacks,
                                       editSetCallbacks=editSetCallbacks,
                                       tipTexts=tipTexts)
    self.displayTable.grid(row=2, column=0, sticky='nsew')

  def updateSpecSelection(self):
  
    presentSpectrumSettings = set([specSetting.ccpnSpectrum  for specSetting in self.specConfigList])
  
    for expt in self.nmrProject.sortedExperiments():
      for spec in expt.sortedDataSources():
        if spec.dataType == 'processed' and spec not in presentSpectrumSettings :
          
          newSpectrumSetting = SpectrumSettings()
          
          newSpectrumSetting.ccpnSpectrum = spec
          
          newSpectrumSetting.peakList = spec.getActivePeakList()
    
          self.specConfigList.append(newSpectrumSetting)
      
  def getLabellingSchemes(self):
    
    return [True, None,] + self.project.sortedLabelingSchemes()
  
  def setAutoLabellingScheme(self, scheme):
    
    self.selectedAutoSpec.setupLabellingScheme(scheme)

    self.updateAutoMatrix()
    
  def getAutoLabellingSchemes(self, notifyScheme=None):
  
    #names = []
    #index = 0
    #
    #scheme = self.selectedAutoSpec.labellingScheme
    #schemes = self.getLabellingSchemes()
    #
    #
    #if schemes:
    #  names = ['Automatic from sample','<None>'] + [sc.name for sc in schemes[2:]]
    #
    #  index = schemes.index(scheme)
    
    names = ['Automatic from sample']
    schemes = [True]
    index = 0
    
    self.autoLabellingPulldown.setup(names, schemes,  index)
    
  def setAutoPeakList(self, peakList):
    
    self.selectedAutoSpec.setupPeakList(peakList)

    self.updateAutoMatrix()

  def selectAutoSpec(self, obj, row, col):
    
    self.selectedAutoSpec = obj

  def setAutoPeakList(self, peakList):
    
    self.selectedAutoSpec.setupPeakList(peakList)

    self.updateAutoMatrix()
       
  def getAutoPeakLists(self,  notifyScheme=None):
  
    names = []
    index = 0
    
    peakList = self.selectedAutoSpec.peakList
    
    peakLists = self.selectedAutoSpec.ccpnSpectrum.sortedPeakLists()


    if peakLists:
      names = [str(pl.serial) for pl in peakLists]

      index = peakLists.index(peakList)
    
    self.autoPeakListPulldown.setup(names, peakLists,  index)

  def changeUse(self, x):
    
    if self.selectedAutoSpec.used is False :
      self.selectedAutoSpec.changeSpectrumUse(True)
    elif self.selectedAutoSpec.used is True :
      self.selectedAutoSpec.changeSpectrumUse(False)
      
    self.updateAutoMatrix()
    
  def updateAutoMatrix(self):
    
    textMatrix = []
    objectMatrix = []
    colorMatrix = []
    
    spectra = self.specConfigList                                     
    
    
    for i, spectrum in enumerate(spectra) :
      
      dataSource = spectrum.ccpnSpectrum
      expt = dataSource.experiment
      
      name = '%s:%s' % (expt.name, dataSource.name)
      
      peakListName = str(spectrum.peakList.serial)
      

      if spectrum.labellingScheme is None :
        schemeName = 'None'
        
      elif spectrum.labellingScheme is True :
        
        schemeName = 'Automatic from sample'
      
      else :
        
        schemeName = spectrum.labellingScheme.name
    
      datum = [i+1, name, peakListName, spectrum.used and 'Yes' or 'No', schemeName]

      textMatrix.append(datum)
      
      
      hexColors = dataSource.analysisSpectrum.posColors
      hexColor = hexColors[int(0.7*len(hexColors))]
      
      if spectrum.used :
      
        colorMatrix.append([hexColor,  hexColor,  hexColor,  hexColor,  hexColor])
        
      else :
        
        colorMatrix.append([hexColor,  None,  None,  None,  None])
      
    self.displayTable.update(textMatrix=textMatrix,objectList=spectra,  colorMatrix=colorMatrix)
    
  def update(self, *x) :
    
    self.updateSpecSelection()
    self.updateAutoMatrix()
    
class SpectrumSettings(object):
  
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
 
