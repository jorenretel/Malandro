
cdef class myResonance :
  
  cdef mySpinSystem mySpinSystem
  
  cdef double CS
  
  cdef str isotopeCode
  
  cdef str atomName 
  
  cdef dict peakDimsLib
  
  cdef dict peakDimsLibIntra
  
  #cdef dict peakDimsLibUnlabelled
  
  cdef object ccpnResonance
  
  cdef int serial
  
  
  
  def __init__(self, mySpinSystem, ccpnResonance):

    self.mySpinSystem = mySpinSystem
    
    self.ccpnResonance = ccpnResonance
    
    self.serial = ccpnResonance.serial
    
    self.isotopeCode = ccpnResonance.isotopeCode
    
    self.atomName = ccpnResonance.assignNames[0]
    
    self.CS = ccpnResonance.findFirstShift().value
    
    self.peakDimsLib = {}
    
    self.peakDimsLibIntra = {}
    
  cdef void addPeakToPeakDimsLib(self,Peak peak, aDimension dim) :
    
    cdef Spectrum spectrum
    
    cdef dict dimsLib
    cdef dict entryForSpectrum
    
    spectrum = peak.spectrum
    
    if peak.intraResidual :
      
      dimsLib = self.peakDimsLibIntra
    
    else :
      
      dimsLib = self.peakDimsLib
    
    if spectrum.name in dimsLib :
      
      entryForSpectrum = dimsLib[spectrum.name]
      
      if dim.dimNumber in entryForSpectrum :
        
        listWithPeaks = entryForSpectrum[dim.dimNumber]
        listWithPeaks.append(peak)
      
      else :
        
        entryForSpectrum[dim.dimNumber] = [peak]
        
    else : 
      
      newlib = {}
      newlib[dim.dimNumber] = [peak]
      
      dimsLib[spectrum.name] = newlib
      
  cdef list getPeaksForSpectrumDim(self, Spectrum spectrum, int dimNumber, intraResidual=False) :
    
    cdef dict entryForSpectrum
    
    peaks = []
    
    if spectrum.name in self.peakDimsLib :
      
      entryForSpectrum = self.peakDimsLib[spectrum.name]
      
      if dimNumber in entryForSpectrum :
        
        peaks.extend( entryForSpectrum[dimNumber] )
      
    if intraResidual is True and spectrum.name in self.peakDimsLibIntra :
      
      entryForSpectrum = self.peakDimsLibIntra[spectrum.name]
      
      if dimNumber in entryForSpectrum :
        
        peaks.extend( entryForSpectrum[dimNumber] )
      
      
    return peaks
        
  def __reduce__(self) :

    return (generalFactory, (type(self),), self.__getstate__())
    
  def __getstate__(self) :

    return (self.mySpinSystem, self.CS, self.isotopeCode, self.atomName, self.serial)
  
  def __setstate__(self, state) :

    self.mySpinSystem, self.CS, self.isotopeCode, self.atomName, self.serial = state
    
  def connectToProject(self) :
    
    self.ccpnResonance = self.mySpinSystem.ccpnResonanceGroup.findFirstResonance(serial=self.serial)
    
    if not self.ccpnResonance :
      
      print 'Error: could not find resonance %s' %str(self.serial)
      return
  
  def getChemicalShift(self) :
    
    return self.CS
  
  def getSpinSystem(self) :
    
    return self.mySpinSystem
  
  def getAtomName(self) :
    
    return self.atomName
  
  def getCcpnResonance(self) :
    
    #return self.ccpnResonance
    
    ccpnResonanceGroup = self.mySpinSystem.getCcpnResonanceGroup()
    
    if ccpnResonanceGroup :
      
      ccpnResonance = ccpnResonanceGroup.findFirstResonance(serial=self.serial)
      
      if ccpnResonance :
        
        return ccpnResonance
      
      else :
      
        print 'Error: could not find resonance %s, might be deleted or merged with other resonance' %str(self.serial)
        
        return None

