
from memops.gui.MessageReporter import showMulti
from ccpnmr.analysis.core.AssignmentBasic import mergeSpinSystems, assignResonanceResidue, findConnectedSpinSystem, clearSeqSpinSystemLinks


def assignSpinSystemResidueMinimal(spinSystem,residue=None):
  
  for resonanceProb in spinSystem.resonanceProbs:
    # NBNB TBD should this happen??? Rasmus Jan 2010
    resonanceProb.delete()
    
  if residue :
    
    molResidue  = residue.molResidue
    #ccpCode     = molResidue.ccpCode
    molType     = molResidue.molType
    
    spinSystem.setResidue(residue)
    spinSystem.setMolType(molType)
    
    for resonance in spinSystem.resonances:
      assignResonanceResidue(resonance,residue)
      
  else:
    spinSystem.setResidue(None)

    for resonance in spinSystem.resonances:
      assignResonanceResidue(resonance,None)
      
def assignSpinSystemstoResidues(spinSystems, residues, strategy=None, guiParent=None) :
  
  ''' input: ccpn objects
     
      output: None
      
      This method is somewhat different than multiple calls to assignmentBasic.assignSpinSystemResidue does.
      If you want to assign more spin systems in one go, you should prevent merging 2 spin systems
      A and B if B possibly needs to be assigned in a next iteration to another residue.
      Also multiple calls to assignmentBasic.assignSpinSystem would cause unnescesarry breaking of sequential links,
      because they don't fit the assignment temporarely. In this funtion, non-valid sequential links are only corrected after
      all assignments have been made. Also, present sequential links are ignored in this function, in the sense that linked
      spin systems are not automatically assigned to neighboring residues because they might conflict with the answer the
      annealing presented. Also no new sequential links are 
  '''
  
  #proposedSRdict = {}
  proposedRSdict = {}
  
  visitedSpinSystems = set()
  visitedResidues = set()
  
  deAssignedSpinSystems = set()
  
  
  # Determine the proposed mapping between spin systems and residues and the other way around. Only a unique mapping between residues and spin systems
  # is used, because it is impossible to assign one spin system to different residues. Although it is technically possible to assign two spin systems
  # to one residue that is not what I want to do (unless there is already a spin system assigned to the residue, in which case I'll ask what to do with it.).
  
  for spinSystem, residue in set(zip(spinSystems,residues)) :
    
    if not spinSystem or not residue :
      continue
    
    if spinSystem in visitedSpinSystems or residue in visitedResidues:
      
      proposedRSdict.pop(residue, None)
      continue
    
    proposedRSdict[residue] = spinSystem
    visitedSpinSystems.add(spinSystem)
    visitedResidues.add(residue)
 
  
  # Figuring out which of the spin systems assigned to the targetted residues are not going to be shuffled around, they potentially have to get
  # merged or correct assignment might already be present
  
  assignSpinSystems = set(proposedRSdict.values())
  
  notMovingRSdict = {}
  
  for residue, spinSystem in proposedRSdict.items() :
    
    nonMovingSpinSystems = []
    placeHolderSpinSystems = []
    
    oldSpinSystems = residue.getResonanceGroups()
    
    for oldSpinSystem in oldSpinSystems :
      
      if oldSpinSystem is spinSystem or oldSpinSystem not in assignSpinSystems:
        
        if oldSpinSystem.resonances :
          
          nonMovingSpinSystems.append(oldSpinSystem)
        
        else :
          
          placeHolderSpinSystems.append(oldSpinSystem)
        
      
    if nonMovingSpinSystems :
      
      if spinSystem in nonMovingSpinSystems :               # The spin system is already assigned to the residue
        
        continue
      
      else :
        
        tempStrategy = strategy or askToMergeSpinSystems(residue,spinSystem,nonMovingSpinSystems+placeHolderSpinSystems,guiParent=guiParent)
        # options are: 'merge','remove','noMerge','skip'
        
        if tempStrategy == 'skip' :
          
          continue
        
        elif tempStrategy == 'noMerge' :
          
          assignSpinSystemResidueMinimal(spinSystem,residue)
          
        elif tempStrategy == 'remove' :
          
          for oldSpinSystem in nonMovingSpinSystems :
            
            assignSpinSystemResidueMinimal(oldSpinSystem, None)
            deAssignedSpinSystems.add(oldSpinSystem)
            
          assignSpinSystemResidueMinimal(spinSystem,residue)
          
          for placeHolderSpinSystem in placeHolderSpinSystems :
            
            mergeSpinSystems(placeHolderSpinSystem,spinSystem)
          
        elif tempStrategy == 'merge' :
          
          assignSpinSystemResidueMinimal(spinSystem,residue)
          
          for oldSpinSystem in placeHolderSpinSystems + nonMovingSpinSystems:
            
            mergeSpinSystems(oldSpinSystem,spinSystem)
             
    else :        
            
      assignSpinSystemResidueMinimal(spinSystem,residue)
      
      for placeHolderSpinSystem in placeHolderSpinSystems :
        
        mergeSpinSystems(placeHolderSpinSystem,spinSystem)
          
  #Removing sequentials links that do not make sense any longer        
  for spinSystem in assignSpinSystems.union(deAssignedSpinSystems) :
    
    connectedSpinSystems = [findConnectedSpinSystem(spinSystem, delta=i) for i in (-1,1)]
    
    residue = spinSystem.getResidue()
    
    if residue :
      
      connectedResidues = [residue.chain.findFirstResidue(seqId = residue.seqId + i) for i in (-1,1)]
      
      for connectedResidue, connectedSpinSystem, i in zip(connectedResidues, connectedSpinSystems, [-1,1]) :
      
        if connectedResidue :
          
          assignedSpinSystems = connectedResidue.getResonanceGroups()
          
          if connectedSpinSystem not in assignedSpinSystems :
            
            clearSeqSpinSystemLinks(spinSystem, delta=i)
            
    else :
      
      for connectedSpinSystem,i in zip(connectedSpinSystems, [-1,1]) :
        
        if connectedSpinSystem :
          
          connectedResidue = connectedSpinSystem.getResidue()
          
          if connectedResidue :
            
            residue = connectedResidue.chain.findFirstResidue(seqId = residue.seqId + (i*-1))
            
            if residue and residue.getResonanceGroups() :
              
              clearSeqSpinSystemLinks(spinSystem, delta=i)
            

    
def askToMergeSpinSystems(residue,newSpinSystem,oldSpinSystems,guiParent=None) :
  
  title = 'Merge'
  amountOfOld = len(oldSpinSystems)
  
  if amountOfOld == 1 :
    
    oldSpinSystem = oldSpinSystems[0]
    
    message = 'You want to assign spin system %s to residue %s , \nbut spin system %s is already assigned to this residue. \nWhat to do?' %(str(newSpinSystem.serial),str(residue.seqCode),str(oldSpinSystem.serial)) 
    texts = ['Merge them','De-assign old spin system', 'Assign, no merge', 'Skip']
  
  else :
    
    message = 'You want to assign spin system %s to residue %s , \nbut there are already %s spin systems assigned to this residue. \nWhat to do?' %(str(newSpinSystem.serial),str(residue.seqCode),str(amountOfOld)) 
    texts = ['Merge them','De-assign old spin system', 'Assign, no merge', 'Skip']  

  objects = ['merge','remove','noMerge','skip']
  
  strategy = showMulti(title, message, texts, objects=objects, parent=guiParent)
          
  return strategy
