
cdef class PeakDimension:
    '''Wraps peak dimension
       attributes:
           peak: parent, the peak this is a dimension of.
           ccpnDim: (ccp.api.nmr.Nmr.PeakDim) peakDim that is wrapped.
           ppmValue: (double) chemical shift in ppm
           tolerance: (double) tolerance of dimension in ppm
           dimNumber: (int) number corresponding to the reference
                      experiment dimension. This is used to connect
                      expected peaks (which get calculated by
                      recursively walking the reference experiment).
                      This attribute is nescesarry to connect the
                      dimensions of expected and real peaks together.
           dim: (int) number connected to this peak dimension. Is only
                used to retrieve the peakDim object (self.ccpnDim)from
                the ccpn analysis project if that connection is broken.
           possibleContributions: (list) Is used to list all resonances
                                  that have a chemical shift within the
                                  tolerance of the dimension.



    '''
    cdef object ccpnDim
    cdef double ppmValue, tolerance

    # I always use dimNumber, which comes from the reference Experiment
    # because this is also the basis of the simulation of the spectra. 'dim'
    # just serves as a key.
    cdef int dimNumber, dim
    cdef list possibleContributions
    cdef Peak peak

    def __init__(self, peak, ccpnDim):
        '''Init.
           args:
               peak: parent, the peak this is
                     a dimension of.
               ccpnDim: (ccp.api.nmr.Nmr.PeakDim) the
                        ccpn analysis peak dimension
        '''

        self.peak = peak
        self.ccpnDim = ccpnDim
        self.ppmValue = ccpnDim.value
        self.dimNumber = ccpnDim.dataDim.expDim.refExpDim.dim
        #self.dim = ccpnDim.dim
        self.tolerance = getAnalysisDataDim(ccpnDim.dataDim).assignTolerance

        # All resonances in the resonanceList that could potentially contribute
        # to this dimension of the peak
        self.possibleContributions = []

    def __reduce__(self):

        return (generalFactory, (type(self),), self.__getstate__())

    def __getstate__(self):

        return (self.dimNumber, self.dim, self.ppmValue, self.peak)

    def __setstate__(self, state):

        self.dimNumber, self.dim, self.ppmValue, self.peak = state

    def connectToProject(self):
        '''Connect this instance again to its counterpart
           in ccpn analysis by setting the ccpnDim attribute.
        '''
        self.ccpnDim = self.peak.ccpnPeak.findFirstPeakDim(dim=self.dim)

    def getChemicalShift(self):
        '''Returns the chemical shift in
           ppm of the peak dimension
        '''

        return self.ppmValue

    def getDimNumber(self):
        '''Get the dimension number of the
           reference experiment dimension.
        '''

        return self.dimNumber

    def getCcpnDimension(self):
        '''Get the ccpn dimension this instance wraps.
        '''
        return self.ccpnDim
