
cdef class aPeak :

  cdef int degeneracy
  
  cdef int degeneracyTemp
  
  cdef object ccpnPeak
  
  cdef list dimensions
  
  cdef aSpectrum spectrum
  
  cdef int serial
  
  cdef int peakListSerial

  cdef object intraResidual
  
  cdef object pyPeak
  
  
  def __init__(self, spectrum, ccpnPeak):
    
    self.ccpnPeak = ccpnPeak
    self.spectrum = spectrum
    self.dimensions = []
    self.degeneracy = 0
    self.degeneracyTemp = 0
    self.intraResidual = False
    self.setupDimensions()
    self.checkForIntraResidualAssignment()
    self.serial = ccpnPeak.serial
    self.peakListSerial = ccpnPeak.peakList.serial
    
  def setupDimensions(self):
    
    ccpnDims = self.ccpnPeak.sortedPeakDims()
    
    self.dimensions = [None]*len(ccpnDims)            # We want the dimensions sorted by .dataDim.expDim.refExpDim.dim (this dimNumber is also used everywhere in the simulation, and it is handy during the matching procedure to have the contributions to the simulated peaks in the same order as dimensions of the real peaks in the spectra)
    
    for dim in ccpnDims :

      dimension = aDimension(self,dim)

      self.dimensions[dimension.dimNumber - 1] = dimension
      
      
  cdef void checkForIntraResidualAssignment(self):
    
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
      
  def __reduce__(self) :

    return (generalFactory, (type(self),), self.__getstate__())
    
  def __getstate__(self) :
    
    return (self.dimensions, self.spectrum, self.serial, self.peakListSerial)
  
  def __setstate__(self, state) :

    self.dimensions, self.spectrum, self.serial, self.peakListSerial = state
    
  def getSpectrum(self):
    
    return self.spectrum
  
  def getDimensions(self) :
    
    return self.dimensions
 