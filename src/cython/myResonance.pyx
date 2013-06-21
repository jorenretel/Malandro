
cdef class myResonance :
  
  cdef mySpinSystem mySpinSystem
  
  cdef double CS
  
  cdef str isotopeCode
  
  cdef str atomName 
  
  cdef dict peakDimsLib
  
  cdef dict peakDimsLibIntra
  
  #cdef dict peakDimsLibUnlabelled
  
  cdef object ccpnResonance
  
  cdef object pyResonance
  
  
  
  def __init__(self, mySpinSystem, ccpnResonance):

    self.mySpinSystem = mySpinSystem
    
    self.ccpnResonance = ccpnResonance
    
    self.isotopeCode = ccpnResonance.isotopeCode
    
    self.atomName = ccpnResonance.assignNames[0]
    
    self.CS = ccpnResonance.findFirstShift().value
    
    self.peakDimsLib = {}
    
    self.peakDimsLibIntra = {}
    
  cdef void addPeakToPeakDimsLib(self, aPeak peak, aDimension dim) :
    
    cdef aSpectrum spectrum
    
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
      
  cdef list getPeaksForSpectrumDim(self, aSpectrum spectrum, int dimNumber, intraResidual=False) :
    
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
    
  cdef void createPythonStyleObject(self) :
    
    self.pyResonance = pyResonance()
    
    self.pyResonance.CS = self.CS
    
    self.pyResonance.atomName = self.atomName
    