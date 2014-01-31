
cdef class dictContainingLists(dict) :
  
  @classmethod  
  def fromDict(cls,dictionary) :
    
    newDictContainingLists = cls()
    
    for key, value in dictionary.items() :
      
      newDictContainingLists[key] = value
      
    return newDictContainingLists
    
  def __init__(self):
    
    pass

  cdef void addValue(self,key, value):  
    
    if key in self :
      self[key].append(value)      
    else :  
      self[key] = [value]
      
  cdef void addList(self, key, value) :
    
    if key in self :
      
      self[key].extend(value)
      
    else :
      
      self[key] = value
      
  def __add__(self,dictContainingLists otherDict):
    
    newDict = dictContainingLists()
    
    for key, value in self.items() :
      
      newDict.addList(key, value[:])

    for key, value in otherDict.items():
      
      newDict.addList(key, value[:])
          
    return newDict
  
  cdef dictContainingLists myCopy(self) :
    
    cdef dictContainingLists newDict
    
    newDict = dictContainingLists()
    
    for key, value in self.items() :
      
      newDict.addList(key, value[:])
      
    return newDict  
  
  cdef set getValueSet(self) :
    
    cdef list alist
    
    newList = []
    
    for alist in self.values() :
      
      newList.extend(alist)
    
    return set(newList)

cdef inline list commonElementInLists(list listOfLists):                                                         # Is only called in matchSpectrum, maybe should be function instead of method
     
  if listOfLists :
      
    return list(set(listOfLists[0]).intersection(*listOfLists))
            
  return []

cdef dict makePrivateCopyOfDictContainingLists(dict dictio):                                                   # Make a copy of the spinsystemDict that can be modified without changing the original dictionary
  
  cdef dict copiedDict
  cdef list copyOfList
  
  copiedDict = {}
  
  for key in dictio :
    
    copyOfList = dictio[key][:]                                 #a copy of each list is made by slicing the list, trick
    copiedDict[key] = copyOfList
    
  return copiedDict

cdef dict mergeDictionariesContainingLists(list dictionaries):
  
  
  cdef list dicts
  cdef dict dictio
  cdef dict newDict
  cdef list newList
  
  
  keys = []
  
  dicts = []
    
  for dictio in dictionaries :
    
    dicts.append(makePrivateCopyOfDictContainingLists(dictio))
    
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

cdef SpinSystemLink emptyLink = SpinSystemLink()

@cython.nonecheck(False)
@cython.wraparound(False)
@cython.boundscheck(False)
#@cython.cdivision(True)  
cdef inline double CcalcDeltaPeakScore(list peakSet,list oldPeaks,list newPeaks):
    
  cdef double peakScoreNew = 0
  cdef double peakScoreOld = 0
  
  #cdef Peak peak
  
  cdef PeakLink pl

  for pl in peakSet :
    
    pl.peak.degeneracyTemp = pl.peak.degeneracy
  
  for pl in oldPeaks :
    
    deg = pl.peak.degeneracy
    peakScoreOld = peakScoreOld + 1.0/pl.peak.degeneracy * pl.preMultipliedScore
    pl.peak.degeneracyTemp -= 1

  for pl in newPeaks :
    
    pl.peak.degeneracyTemp += 1
    
  for pl in newPeaks : 
    
    peakScoreNew += 1.0/pl.peak.degeneracyTemp * pl.preMultipliedScore
    
  return peakScoreNew - peakScoreOld

cdef inline double scorePeak(list peakDimensions,list resonances) :
  
  cdef PeakDimension dim
  cdef Resonance resonance
  cdef double k, z, top, summation, delta, tolerance
  
  k = 0.4
  top = 1.0/(k**2-1)
  #top = -1.0/(k**2-1)

  summation = 0.0
  
  for dim, resonance in zip(peakDimensions, resonances) :

    delta = dim.ppmValue - resonance.CS
    tolerance = dim.tolerance
    summation += ((delta/tolerance)**2-1)*top
    
  z = summation/len(peakDimensions)
  
  return z

cdef dict runAminoAcidTyping(resonanceGroup, shiftList, ccpnChain, double minTypeScore) :
  '''Find all residue types the spin system could be assigned to.
     The algorithm in CCPN analysis is used for this task.
     kwargs:  minTypeScore: cut-off value (percentage). Everything
                           scoring higher is consider a possible
                           residue type assignment.
     
  '''
  
  shifts = []
  for resonance in resonanceGroup.resonances:
    if resonance.isotopeCode in ('1H',  '13C',  '15N') :
      shift = resonance.findFirstShift(parentList=shiftList)
      if shift:
        shifts.append(shift)

  scores = getShiftsChainProbabilities(shifts, ccpnChain)
  total = sum(scores.values())
  
  scoreDict = {}
  
  if total:
      
    for ccpCode, score in scores.items() :
      
      relativeScore = score/total * 100.0         # minTypeScore is a percentage, therefor so should relativeScore 
      
      if relativeScore >= minTypeScore :
        
        scoreDict[ccpCode] = relativeScore
        
  return scoreDict