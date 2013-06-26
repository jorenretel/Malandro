
from ccp.api.nmr.Nmr import ResonanceGroup as CCPNResonanceGroup


class pyDataModel(object):

  def __init__(self) :
      
    self.spectra = []

    self.chain = None

    self.spinSystems = {}

    self.amountOfRepeats = None
    

    
  def __getstate__(self):
    
    state = dict(self.__dict__)
    
    if 'spectra' in state :
      print 'kiiii'
      del state['spectra']
      
#    if 'chain' in state :
#      print 'k'
#      del state['chain']
#      
#    if 'spinSystems' in state :
#      print 'k'
#      del state['spinSystems']
      
      
      
    
      
   # if 'spinSystems' in state :
   #   del state['spinSystems']
      
    return state
    
    
    
class pyChain(object):

  def __init__(self, ccpnChain) :
    
    self.ccpnChain = ccpnChain
    
    self.key = (ccpnChain.molSystem.code, ccpnChain.code)
  
    self.residues = []
    
    
  def __getstate__(self):
    
    print 'pychain'
    
    self.ccpnChain = None
    
    state = dict(self.__dict__)
      
    return state
      
  
  def connectToProject(self, project) :
    
    molSystemCode, chainCode = self.key
    
    molSystem = project.findFirstMolSystem(code=molSystemCode)
    
    if molSystem :
    
      self.ccpnChain = molSystem.findFirstChain(code=chainCode)
      
      if self.ccpnChain :
        
        for res in self.residues :
          
          pass
          #res.connectToProject()
        
        return
      
    print 'It seems like the ' + molSystemCode + ' ' + chainCode + ' chain was removed.'
    
    


class pyResidue(object):

  def __init__(self) :
  
    self.seqCode = None
    
    self.ccpCode = None
    
    self.solutions = []                                                     # !!!!
    
    self.linkDict = {}                                                      # !!!!
    
    self.intraDict = {}
    
    self.userDefinedSolution = None
    
  def __getstate__(self):
    
    print 'pyResidue'
    
    state = dict(self.__dict__)
      
    return state
    
class pySpectrum(object):
  
  def __init__(self) :
    
    self.name = None
    
    self.peaks = []
    
  def __getstate__(self):
    
    state = dict(self.__dict__)
    
    if 'peaks' in state :
      
      print 'j'
      
      del state['peaks']
      
    return state
    
    
class pyPeak(object):
  
  def __init__(self) :
    
    self.dimensions = []
    
    self.serial = None
    
    self.peakListSerial = None
    
    self.spectrum = None                                                        # This could be changed to just the name of the spectrum
    
  def __getstate__(self):
    
    print 'pyPeak'
    
    state = dict(self.__dict__)
      
    return state
    

class pyDimension(object):
  
  def __init__(self) :
    
    self.ppmValue = None
    
    self.dimNumber = None
    
  def __getstate__(self):
    
    print 'pyDimension'
    
    state = dict(self.__dict__)
      
    return state
    
class pySpinSystemLink(object):
  
  def __init__(self):
    
    self.spinSystem1 = None                                                   #!!!                           
    self.spinSystem2 = None                                                   #!!!    
    self.realPeaks = []                                                               #!!!    
    self.simulatedPeaks = []                                                      #!!!    
    
  def __getstate__(self):
    
    print 'link'
    
    state = dict(self.__dict__)
    
    if 'spinSystem1' in state :
      
      del state['spinSystem1']
      
    if 'spinSystem2' in state :
      
      del state['spinSystem2']
    

      
    return state
    
class pySimulatedPeak(object):

  def __init__(self):
    
    self.spectrum = None
    
    self.colabelling = None
    
    self.simulatedPeakContribs = []
    
  def __getstate__(self):
    
    print 'simulatedpeak'
    
    state = dict(self.__dict__)
      
    return state
    
class pySimulatedPeakContrib(object):
  
  def __init__(self):
      
    self.ccpCode = None
    
    self.atomName = None
    
    self.dimNumber = None
    
    self.firstOrSecondResidue = None
    
  def __getstate__(self):
    
    print 'simpeakcontrib'
    
    state = dict(self.__dict__)
      
    return state
    
class pySpinSystem(object):
  
  def __init__(self):
    
    self.spinSystemNumber = None
    
    self.ccpCode = None
    
    self.resonanceDict = {}
        
    self.ccpnSeqCode = None
    
    self.isJoker = None
    
    self.solutions = []                                                                 
    
    self.userDefinedSolutions = []             

    self.tentativeCcpCodes = []
  
    self.tentativeSeqCodes = []
    
    self.allowedResidues = set()
    
  def getCCPNObject(self, nmrProject) :
    
    return nmrProject.findFirstResonanceGroup(serial=self.spinSystemNumber)
  
  def __getstate__(self):
    
    state = dict(self.__dict__)
    
    print 'lllllll'
    
#    if 'solutions' in state :
#      
#      del state['solutions']
#      
#    if 'allowedResidues' in state :
#      
#      del state['allowedResidues']
      
    if 'resonanceDict' in state :
      
      del state['resonanceDict']
      
    return state
    
    

    
class pyResonance(object):
  
  def __init__(self):

    self.CS = None
    
    self.atomName = None
    
  def __getstate__(self):
    
    print 'resonance'
    
    state = dict(self.__dict__)
      
    return state

    
    
  
  
    
    

  
