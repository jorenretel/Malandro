'''Contains the PeakLink class.
'''

cdef class PeakLink:
    '''peakLinks record the match between expected/simualted peaks
       and real peaks in the spectrum.
       attributes:
         peak: (Peak) the real peak from the spectrum.
         simulatedPeak: (SimulatedPeak) the expected peak.
         score: the score of how well the peak matches to the chemical
             shifts of the resonances making up the dimensional
             contributions. This is scored (with a flat-bottom
             potential, in the case of this algorithm).
         preMultipliedScore: exactly the same as peak score * link score.
             this attribute is only there to make the monte carlo
             procedure faster and easier. Normally we would do during
             Monte Carlo:
             (peakLink1.score + peakLink2.score + ....)*spinSystemLink.score
             but we can already do this multiplication before because:
             (a+b)*c == a*c + b*c. This does not only save a
             multiplication during the run but also makes the loop
             simpler since the lookup of spinSystemLink.score can be
             ommited, and therefor it is also not nesesary to keep.
         resonances: The resonances corresponding to the dimensions
             of the (simulated) peak. This should not be confused
             with the assignment of the peak.
             Only if the spinSystemLink this peakLink belongs to is
             correct (i.e. the two spin systems connected by the
             spinSystemLink are really sequential spin systems
             representing a specific pair of sequential residues in the
             sequence) these resonances can be seen as the correct
             assignments for the peak dimensions.

    '''

    cdef Peak peak
    cdef SimulatedPeak simulatedPeak
    cdef double score, preMultipliedScore
    cdef list resonances

    def __init__(self, Peak peak, SimulatedPeak simulatedPeak,
                 list resonances, double score):
        '''Init
             args:
               peak: (Peak) the real peak
               simulatedPeak: (SimulatedPeak) the simulatedPeak
               score: how well the peak fits the simualted peak
               resonances: the resonances corresponding to the
                           peak dimensions of the peak.
        '''

        self.peak = peak
        self.simulatedPeak = simulatedPeak
        self.score = score
        self.preMultipliedScore = score
        self.resonances = resonances

    def __reduce__(self):

        return (generalFactory, (type(self),), self.__getstate__())

    def __getstate__(self):

        return (self.peak, self.simulatedPeak, self.resonances, self.score)

    def __setstate__(self, state):

        self.peak, self.simulatedPeak, self.resonances, self.score = state

    def getPeak(self):
        '''Returns the real peak'''

        return self.peak

    def getSimulatedPeak(self):
        '''Returns the expected/simualted peak.'''
        return self.simulatedPeak

    def getResonances(self):
        '''Returns the resonances corresponding to the
           (simulated) peak dimensions.
        '''
        return self.resonances
