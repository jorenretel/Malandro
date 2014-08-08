
cdef class SpinSystemLink:

    cdef int score
    cdef list realPeaks, simulatedPeaks, notFoundSimulatedPeaks
    cdef list peakLinks, activePeakLinks, notFoundPeakLinks
    cdef Residue residue1, residue2
    cdef SpinSystem spinSystem1, spinSystem2

    def __init__(self, residue1=None, residue2=None,
                 spinSystem1=None, spinSystem2=None):

        self.score = 0
        self.residue1 = residue1
        self.residue2 = residue2
        self.spinSystem1 = spinSystem1
        self.spinSystem2 = spinSystem2
        self.simulatedPeaks = []
        self.realPeaks = []
        self.notFoundSimulatedPeaks = []
        self.peakLinks = []
        self.activePeakLinks = []
        self.notFoundPeakLinks = []

    cdef void addPeak(self, Peak realPeak, SimulatedPeak simPeak,
                      list resonances, double score):

        if realPeak is not None:

            newPeakLink = PeakLink(realPeak, simPeak, resonances, score)
            self.peakLinks.append(newPeakLink)

        else:

            self.notFoundSimulatedPeaks.append(simPeak)

    def __reduce__(self):

        return (generalFactory, (type(self),), self.__getstate__())

    def __getstate__(self):

        # , self.notFoundPeakLinks)#, self.notFoundPeakLinks)
        return (self.peakLinks, self.notFoundSimulatedPeaks,
                self.spinSystem1, self.spinSystem2)

    def __setstate__(self, state):

        # , self.notFoundPeakLinks = state #self.notFoundPeakLinks
        (self.peakLinks,
         self.notFoundSimulatedPeaks,
         self.spinSystem1,
         self.spinSystem2) = state

    def getContributingResonances(self):
        '''Return the resonances that contribute to dimensions
           of peaks of all peakLinks.

        '''

        resonances = []

        for pl in self.peakLinks:

            resonances.extend(pl.getResonances())

        return set(resonances)

    def getResonancesContributingToActivePeakLinks(self):
        '''Return the resonances that contribute to the
           activePeakLinks. Because a percentage of the peaks, and
           thereby a percentage of peakLinks might be left out in a run
           of the annealing, this function returns only those resonance
           contributing to dimensions of the used (active) peaks.

        '''

        resonances = []

        for pl in self.activePeakLinks:

            resonances.extend(pl.getResonances())

        return set(resonances)

    def determineScore(self):
        '''Sets self.score to the amount of resonances contributing
           to the link. Also sets peakLink.preMultipliedScore so
           that does not have to be done during the annealing any longer.

        '''

        cdef PeakLink peakLink

        self.score = len(self.getContributingResonances())
        activeScore = len(self.getResonancesContributingToActivePeakLinks())

        for peakLink in self.activePeakLinks:

            peakLink.preMultipliedScore = peakLink.score * activeScore

    # TODO: change this, first of all it is not used and second of all
    # self.notFoundPeakLinks is empty now.
    def getAllResonances(self):

        resonances = []

        for pl in self.peakLinks:

            resonances.extend(pl.getResonances())

        for pl in self.notFoundPeakLinks:

            resonances.extend(pl.getResonances())

        return set(resonances)

    # TODO: update, simulatedPeaks is not used anymore
    def getSimulatedPeaks(self):

        return self.simulatedPeaks

    def getRealPeaks(self):

        return self.realPeaks

    def getPeakLinks(self):

        return self.peakLinks

    def getNotFoundPeakLinks(self):

        cdef SimulatedPeak simPeak
        cdef SimulatedPeakContrib contrib
        cdef SpinSystem spinSystem

        notFoundPeakLinks = []

        spinSystems = {1: self.spinSystem1, 2: self.spinSystem2}

        for simPeak in self.notFoundSimulatedPeaks:

            resonances = []

            for contrib in simPeak.simulatedPeakContribs:

                spinSystem = spinSystems[contrib.firstOrSecondResidue]

                resonances.append(
                    spinSystem.getResonanceForAtomName(contrib.atomName))

            notFoundPeakLinks.append(PeakLink(None, simPeak, resonances, 0.0))

        return notFoundPeakLinks

    def getAllPeakLinks(self):

        # return self.peakLinks + self.notFoundPeakLinks

        return self.getPeakLinks() + self.getNotFoundPeakLinks()

    def getResidues(self):

        return (self.residue1, self.residue2)

    def getSpinSystems(self):

        return (self.spinSystem1, self.spinSystem2)
