
cdef class Chain:

    cdef public list residues
    cdef dict residuesByCcpCode
    cdef object ccpnChain
    cdef str ccpnChainCode, molSystemCode

    def __init__(self, ccpnChain):
        '''Init chain. Creates a wrapper around the ccpn chain.
           kwargs: ccpn chain instance

        '''
        self.ccpnChain = ccpnChain
        self.ccpnChainCode = ccpnChain.code
        self.molSystemCode = ccpnChain.molSystem.code
        self.residues = []
        self.residuesByCcpCode = {}
        self.setupResidues()
        self.linkResiduesTogether()
        self.addDummyResiduesAtEnds()

    cdef void setupResidues(self):
        '''Sets up all residues in the chain and stores them in
           self.residues and in a dictionary self.residuesByCcpCode,
           where the keys are ccpCodes and the values lists with
           residues.

        '''

        for res in self.ccpnChain.sortedResidues():

            newresidue = Residue(self, res)
            self.residues.append(newresidue)

            if res.ccpCode in self.residuesByCcpCode:

                self.residuesByCcpCode[res.ccpCode].append(newresidue)

            else:

                self.residuesByCcpCode[res.ccpCode] = [newresidue]

    cdef void addDummyResiduesAtEnds(self):
        '''Put a residue before the beginning and after the end. Why
           do this: this saves a lot of checks during the annealing
           procedure for not getting a None-type object for
           .previousResidue or .nextResidue. It is not added to
           self.residues, so normally you don't even notice those
           dummy residues are there at all.

        '''
        cdef Residue res
        cdef SpinSystem spinSystem
        cdef Residue firstResidue
        cdef Residue lastResidue

        firstResidue = self.residues[0]
        lastResidue = self.residues[-1]

        res = Residue(self, None)
        spinSystem = SpinSystem()
        #spinSystem.isJoker = True

        res.currentSpinSystemAssigned = spinSystem

        firstResidue.previousResidue = res
        lastResidue.nextResidue = res

    cdef void linkResiduesTogether(self):
        '''This sets up the links between a residues and its
           neighbors in the chain in the instance variables
           'nextResidue' and 'previousResidue'.

        '''
        cdef list residues
        cdef int i
        cdef Residue res, nextResidue

        residues = self.residues

        for i,  res in enumerate(residues[:-1]):

            nextResidue = residues[i + 1]

            res.nextResidue = nextResidue

            nextResidue.previousResidue = res

    def __reduce__(self):

        return (generalFactory, (type(self),), self.__getstate__())

    def __getstate__(self):

        return (self.residues, self.molSystemCode, self.ccpnChainCode)

    def __setstate__(self, state):

        self.residues, self.molSystemCode, self.ccpnChainCode = state

    def connectToProject(self, project):
        '''(Re)set the connection to the chain in the analysis project,
            i.e. set self.ccpnChain atribute.
        '''
        cdef Residue residue

        molSystem = project.findFirstMolSystem(code=self.molSystemCode)

        if molSystem:

            ccpnChain = molSystem.findFirstChain(code=self.ccpnChainCode)

            if ccpnChain:

                self.ccpnChain = ccpnChain

        if not self.ccpnChain:

            print ('It seems like the %s %s chain was removed.'
                   % (self.molSystemCode, self.ccpnChainCode))
            return

        for residue in self.residues:

            residue.connectToProject()

    def getResidues(self):
        '''Returns the list of residues.
           Necessary for access from python.
        '''

        return self.residues
