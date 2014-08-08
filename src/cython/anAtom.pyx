

cdef class Atom:

    cdef object ccpnAtom
    cdef str atomName
    cdef Residue residue
    cdef list assignmentPossibilityDimensions
    cdef dict labelInfoTemp

    def __init__(self, residue, ccpnAtom):

        self.residue = residue
        self.ccpnAtom = ccpnAtom
        self.atomName = ccpnAtom.chemAtom.name
        self.assignmentPossibilityDimensions = []
        self.labelInfoTemp = {}

    def __reduce__(self):

        return (generalFactory, (type(self),), self.__getstate__())

    def __getstate__(self):

        return (self.atomName, self.residue)

    def __setstate__(self, state):

        self.atomName, self.residue = state

    def connectToProject(self):

        self.ccpnAtom = self.residue.ccpnResidue.findFirstAtom(
            name=self.atomName)
