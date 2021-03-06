
cdef class SimulatedPeakContrib:

    cdef str ccpCode, atomName, isotopeCode
    cdef int dimNumber, firstOrSecondResidue
    cdef Residue residue

    def __init__(self):

        self.ccpCode = None

        self.atomName = None

        self.isotopeCode = None

    def __reduce__(self):

        return (generalFactory, (type(self),), self.__getstate__())

    def __getstate__(self):

        return (self.ccpCode, self.atomName,
                self.dimNumber, self.firstOrSecondResidue)

    def __setstate__(self, state):

        (self.ccpCode,
         self.atomName,
         self.dimNumber,
         self.firstOrSecondResidue) = state

    def getCcpCode(self):

        return self.ccpCode

    def getAtomName(self):

        return self.atomName

    def getDimNumber(self):

        return self.dimNumber

    def getResidue(self):

        return self.residue
