
cdef class Peak:

    cdef int degeneracy, degeneracyTemp
    cdef object ccpnPeak, pyPeak
    cdef list dimensions
    cdef Spectrum spectrum
    cdef int serial, peakListSerial
    cdef bint intraResidual

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
        # We want the dimensions sorted by .dataDim.expDim.refExpDim.dim (this
        # dimNumber is also used everywhere in the simulation, and it is handy
        # during the matching procedure to have the contributions to the
        # simulated peaks in the same order as dimensions of the real peaks in
        # the spectra)
        self.dimensions = [None] * len(ccpnDims)

        for dim in ccpnDims:

            dimension = PeakDimension(self, dim)
            self.dimensions[dimension.dimNumber - 1] = dimension

    cdef void checkForIntraResidualAssignment(self):

        lists = []

        for dim in self.ccpnPeak.peakDims:

            spinSystems = []

            for contrib in dim.peakDimContribs:

                if contrib.resonance.resonanceGroup:

                    spinSystems.append(contrib.resonance.resonanceGroup)

            lists.append(spinSystems)

        intra = False
        for element in lists[0]:

            inAllLists = True

            for list in lists[1:]:

                if element not in list:

                    inAllLists = False

            if inAllLists:

                intra = True

        self.intraResidual = intra

    def __reduce__(self):

        return (generalFactory, (type(self),), self.__getstate__())

    def __getstate__(self):

        return (self.dimensions, self.spectrum, self.serial, self.peakListSerial)

    def __setstate__(self, state):

        self.dimensions, self.spectrum, self.serial, self.peakListSerial = state

    def getSpectrum(self):

        return self.spectrum

    def getDimensions(self):

        return self.dimensions

    def getCcpnPeak(self):

        return self.ccpnPeak

    def getSerial(self):

        return self.serial

    def connectToProject(self):

        try:

            self.ccpnPeak = self.spectrum.ccpnSpectrum.findFirstPeakList(
                serial=self.peakListSerial).findFirstPeak(serial=self.serial)

        except AttributeError:

            print 'Error: Cannot find peak list %s in spectrum %s' % (str(self.peakListSerial), self.spectrum.name)

            return

        if not self.ccpnPeak:

            print 'Error: Cannot find peak %s from peak list %s in spectrum %s' % (str(self.serial), str(self.peakListSerial), self.spectrum.name)

            return

        for dimension in self.dimensions:

            dimension.connectToProject()
