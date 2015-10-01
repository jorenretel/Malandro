
# Is only called in matchSpectrum, maybe should be function instead of method
cdef inline list commonElementInLists(list listOfLists):

    if listOfLists:

        return list(set(listOfLists[0]).intersection(*listOfLists))

    return []

cdef SpinSystemLink emptyLink = SpinSystemLink()


@cython.nonecheck(False)
@cython.wraparound(False)
@cython.boundscheck(False)
#@cython.cdivision(True)
cdef inline double CcalcDeltaPeakScore(list peakSet,
                                       list oldPeaks,
                                       list newPeaks):

    cdef double peakScoreNew = 0
    cdef double peakScoreOld = 0

    # cdef Peak peak

    cdef PeakLink pl

    for pl in peakSet:

        pl.peak.degeneracyTemp = pl.peak.degeneracy

    for pl in oldPeaks:

        deg = pl.peak.degeneracy
        peakScoreOld = peakScoreOld + 1.0 / \
            pl.peak.degeneracy * pl.preMultipliedScore
        pl.peak.degeneracyTemp -= 1

    for pl in newPeaks:

        pl.peak.degeneracyTemp += 1

    for pl in newPeaks:

        peakScoreNew += 1.0 / pl.peak.degeneracyTemp * pl.preMultipliedScore

    return peakScoreNew - peakScoreOld

cdef inline double scorePeak(list peakDimensions, list resonances):

    cdef PeakDimension dim
    cdef Resonance resonance
    cdef double k, z, top, summation, delta, tolerance
    cdef int shiftListSerial

    dim = peakDimensions[0]
    shiftListSerial = dim.peak.spectrum.shiftListSerial

    k = 0.4
    top = 1.0 / (k ** 2 - 1)
    #top = -1.0/(k**2-1)

    summation = 0.0

    for dim, resonance in zip(peakDimensions, resonances):

        delta = dim.ppmValue - resonance.getShift(shiftListSerial)
        tolerance = dim.tolerance
        summation += ((delta / tolerance) ** 2 - 1) * top

    z = summation / len(peakDimensions)

    return z

cdef dict runAminoAcidTyping(resonanceGroup, ccpnChain, double minTypeScore):
    '''Find all residue types the spin system could be assigned to.
       The algorithm in CCPN analysis is used for this task.
       kwargs:  minTypeScore: cut-off value (percentage). Everything
                             scoring higher is consider a possible
                             residue type assignment.

    '''

    shifts = []
    for resonance in resonanceGroup.resonances:
        if resonance.isotopeCode in ('1H',  '13C',  '15N'):
            shift = resonance.findFirstShift()
            if shift:
                shifts.append(shift)

    scores = getShiftsChainProbabilities(shifts, ccpnChain)
    total = sum(scores.values())

    scoreDict = {}

    if total:

        for ccpCode, score in scores.items():

            # minTypeScore is a percentage, therefor so should relativeScore
            relativeScore = score / total * 100.0

            if relativeScore >= minTypeScore:

                scoreDict[ccpCode] = relativeScore

    return scoreDict
