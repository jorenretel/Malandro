
cdef class Resonance:

    cdef SpinSystem spinSystem
    cdef str isotopeCode, atomName
    cdef dict peakDimsLib, peakDimsLibIntra, shifts
    cdef object ccpnResonance
    cdef int serial

    def __init__(self, mySpinSystem, ccpnResonance):

        self.spinSystem = mySpinSystem
        self.ccpnResonance = ccpnResonance
        self.serial = ccpnResonance.serial
        self.isotopeCode = ccpnResonance.isotopeCode
        self.atomName = ccpnResonance.assignNames[0]
        #self.CS = ccpnResonance.findFirstShift().value
        self.peakDimsLib = {}
        self.peakDimsLibIntra = {}
        self.shifts = self.getShifts()

    cdef dict getShifts(self):

        cdef dict shiftDict

        shiftDict = {}
        shifts = self.ccpnResonance.sortedShifts()
        for shift in shifts:
            shiftListSerial = shift.parentList.serial
            shiftDict[shiftListSerial] = shift.value
        return shiftDict

    cdef double getShift(self, int serial):

        cdef list shifts

        if not serial:
            shifts = self.shifts.values()
            if shifts:
                return shifts[0]
            #else:
            #    return None

        elif serial in self.shifts:
            return self.shifts[serial]

        #else:
        #    return None

    cdef void addPeakToPeakDimsLib(self, Peak peak, PeakDimension dim):

        cdef Spectrum spectrum
        cdef dict dimsLib, entryForSpectrum

        spectrum = peak.spectrum

        if peak.intraResidual:

            dimsLib = self.peakDimsLibIntra

        else:

            dimsLib = self.peakDimsLib

        if spectrum.name in dimsLib:

            entryForSpectrum = dimsLib[spectrum.name]

            if dim.dimNumber in entryForSpectrum:

                listWithPeaks = entryForSpectrum[dim.dimNumber]
                listWithPeaks.append(peak)

            else:

                entryForSpectrum[dim.dimNumber] = [peak]

        else:

            newlib = {}
            newlib[dim.dimNumber] = [peak]

            dimsLib[spectrum.name] = newlib

    cdef list getPeaksForSpectrumDim(self, Spectrum spectrum,
                                     int dimNumber,
                                     intraResidual=False):

        cdef dict entryForSpectrum

        peaks = []

        if spectrum.name in self.peakDimsLib:

            entryForSpectrum = self.peakDimsLib[spectrum.name]

            if dimNumber in entryForSpectrum:

                peaks.extend(entryForSpectrum[dimNumber])

        if intraResidual is True and spectrum.name in self.peakDimsLibIntra:

            entryForSpectrum = self.peakDimsLibIntra[spectrum.name]

            if dimNumber in entryForSpectrum:

                peaks.extend(entryForSpectrum[dimNumber])

        return peaks

    def __reduce__(self):

        return (generalFactory, (type(self),), self.__getstate__())

    def __getstate__(self):

        return (self.spinSystem, self.shifts, self.isotopeCode,
                self.atomName, self.serial)

    def __setstate__(self, state):

        self.spinSystem, self.shifts, self.isotopeCode, self.atomName, self.serial = state

    def connectToProject(self):

        self.ccpnResonance = self.getCcpnResonance()

    def getChemicalShift(self, shiftListSerial=0):

        return self.getShift(shiftListSerial)

    def getSpinSystem(self):

        return self.spinSystem

    def getAtomName(self):

        return self.atomName

    def getCcpnResonance(self):

        # return self.ccpnResonance

        ccpnResonanceGroup = self.spinSystem.getCcpnResonanceGroup()

        if ccpnResonanceGroup:

            ccpnResonance = ccpnResonanceGroup.findFirstResonance(
                serial=self.serial)

            if ccpnResonance:

                return ccpnResonance

            else:

                print 'Error: could not find resonance %s, might be deleted or merged with other resonance. It is advisable to re-run the algorithm since the state of the project has changed.' % str(self.serial)
