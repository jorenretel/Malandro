
cdef class SpinSystem:
    '''Wrapper around resonceGroup object in CCPN.
       attributes:
         DataModel: parent, root of the model
         ccpnResonanceGroup: (ccp.api.nmr.Nmr.ResonanceGroup) ccpn
             object that is wrapped.
         allowedResidues: (set) all residues this spin system can be
              assigned to according to the multiple levels of assignment
              in the ccpn analysis project and the settings configured
              by the user.
         ccpCodes: (set) The set of ccpCodes fitting to the
              allowedResidues
         resonanceDict: (dict) All the resonances that make up the
              spin system as values and their assignName as the keys.
         userDefinedSolutions: (list) The user can shuffle around
              assignments in the GUI. This is the attribute that holds
              this user generated assignment.
         solutions: (list) List of the residue assignment in each one of
              the runs.
         currentResidueAssignment: (Residue) This is the attribute that
              gets varied during the monte Carlo procedure.
         isJoker: (boolean) True if the spin system is a Joker. Jokers
              are necessary to guarantee that in each step of the
              optimization every residue can be assigned to a
              spin system. This is important. Because two spin
              systems rather than two residues are picked randomly
              to exchange assignment. If a residue would be left
              without a spin system assignment, it would be
              excluded during the rest of the run.
          exchangeSpinSystems: (list) All other spin spin systems this
              spin system can ever be exchanged with during the monte
              carlo procedure.
          resonancesByAtomSiteName: (dict) AtomSites are groups of types
              of nuclei used to described a part of magnetization
              transfer pathway. For instance 'Caro' for aromatic
              carbons. This is the dict
              {'Caro':[resonance1, resonance2....]}
          spinSystemNumber: (int) is the same as the serial of the
              ccpn resonanceGroup. Is used to connect back to this
              resonanceGroup.



    '''

    cdef myDataModel DataModel
    cdef Residue currentResidueAssignment
    cdef int spinSystemNumber
    cdef list exchangeSpinSystems
    cdef object ccpnResonanceGroup
    cdef bint isJoker
    cdef public list solutions, userDefinedSolutions
    cdef public set allowedResidues, ccpCodes
    cdef dict resonanceDict, resonancesByAtomSiteName
    cdef int[:] allowedResidueView

    def __init__(self, DataModel=None, ccpnResonanceGroup=None,
                 allowedResidues=None, ccpCodes=None):
        '''kwargs: DataModel, parent, root of the model
                   ccpnResonanceGroup:  resonanceGroup this spin system
                                        corresponds to in CCPN.


        '''
        self.DataModel = DataModel
        self.ccpnResonanceGroup = ccpnResonanceGroup
        self.allowedResidues = allowedResidues or set()
        self.ccpCodes = ccpCodes or set()
        self.resonanceDict = {}
        self.userDefinedSolutions = []
        self.currentResidueAssignment = None
        self.isJoker = False
        self.solutions = []
        self.exchangeSpinSystems = []
        self.resonancesByAtomSiteName = {}
        self.setupAllowedResidueView()
        if not ccpnResonanceGroup:
            self.isJoker = True
            return
        self.spinSystemNumber = ccpnResonanceGroup.serial
        self.setupResonances()
        self.groupResonancesByAtomSite()

    cdef set getCcpCodes(self):
        '''Return three-letter amino acid code this spin system could
          be assigned to.

        '''
        return self.ccpCodes

    cdef void setupResonances(self):
        '''Sets up all resonances belonging in the spin system based on
           the resonances in the resonanceGroup in CCPN. Resonances are
           stored by their atomName in a dict.

        '''
        for resonance in self.ccpnResonanceGroup.resonances:

            if resonance.getAssignNames() and resonance.findFirstShift():

                newResonance = Resonance(self, resonance)
                self.resonanceDict[newResonance.atomName] = newResonance

    cdef void setupAllowedResidueView(self):
        ''' Sets up a memory view containing bints. Every bint
            corresponts to a residue in the sequence. If the spin
            system could, based on its residue type, be assigned
            to this residue the bint is set to True.

        '''

        cdef Chain chain
        cdef int resNumber
        cdef Residue residue

        if not self.DataModel or not self.DataModel.chain:

            return

        chain = self.DataModel.chain
        cythonArray = cvarray(shape=(len(chain.residues) + 1,),
                              itemsize=sizeof(bint), format="i")
        self.allowedResidueView = cythonArray
        self.allowedResidueView[:] = False

        for residue in self.allowedResidues:

            resNumber = residue.seqCode
            self.allowedResidueView[resNumber] = True

    # Used in matchSpectrum()
    cdef Resonance getResonanceForAtomName(self, str atomName):
        '''Returns a resonace for a atomName.
           kwargs: atomName: name of the atom, for instance 'CA'.

        '''

        if atomName in self.resonanceDict:

            return self.resonanceDict[atomName]

        else:

            return None

    cdef list getResonancesForAtomSiteName(self, str atomSiteName):
        '''Returns resonances for AtomSiteName.
           kwargs: atomSiteName: string describing atomSite,
                                 for instance Cali.

        '''

        return self.resonancesByAtomSiteName.get(atomSiteName, [])

    cdef void setupExchangeSpinSystems(self):
        '''Sets up a list of all spin system this spin system could
           possibly swap out residue assignment with during the Monte
           Carlo procedure. This list is later used in the Monte Carlo
           procedure to pull a random member out off, to attempt a
           switch of residue assignment.

        '''

        cdef SpinSystem spinSys
        cdef dict spinSystemDict

        spinSystems = self.DataModel.getSpinSystemSet()

        for spinSys in spinSystems:
            if (not spinSys is self and
                not (self.isJoker and spinSys.isJoker) and
                self.allowedResidues & spinSys.allowedResidues):

                self.exchangeSpinSystems.append(spinSys)

    cdef void groupResonancesByAtomSite(self):  # TODO: finish
        ''' Sets up resonancesByAtomSiteName dict. AtomSites are used in
            CCPN in de description of magnetization transfer pathways.
            Used later to find possibly correlated resonances in an
            experiment,so the corresponding spectrum can be searched
            for peaks.

        '''
        cdef Resonance resonance

        resonancesByAtomSiteName = self.resonancesByAtomSiteName
        resonancesByName = self.resonanceDict

        ccpCodes = self.getCcpCodes()

        HAs = []
        HBs = []
        CXs = []
        Calis = []
        Caros = []
        COs = []

        resonancesByAtomSiteName['CA'] = []
        resonancesByAtomSiteName['CO'] = []

        # C,F,Br,H,P,N, CX
        for name, resonance in resonancesByName.items():

            elementSymbol = resonance.ccpnResonance.isotope.chemElement.symbol

            resonancesByAtomSiteName[elementSymbol] = resonancesByAtomSiteName.get(elementSymbol, []) + [resonance]

            if elementSymbol == 'C':

                CXs.append(resonance)

        # CA
        if 'CA' in resonancesByName:
            resonancesByAtomSiteName['CA'] = [resonancesByName['CA']]

        # CB
        if 'CB' in resonancesByName:
            resonancesByAtomSiteName['CB'] = [resonancesByName['CB']]

        # CO
        if 'C' in resonancesByName:
            COs.append(resonancesByName['C'])

        if ('Gln' in ccpCodes or 'Glu' in ccpCodes) and 'CD' in resonancesByName:
            COs.append(resonancesByName['CD'])

        if ('Asp' in ccpCodes or 'Asn' in ccpCodes) and 'CG' in resonancesByName:
            COs.append(resonancesByName['CG'])

        # HA
        for name in ['HA', 'HA1', 'HA2', 'HA3']:

            if name in resonancesByName:
                HAs.append(resonancesByName[name])

        # HB
        for name in ['HB', 'HB1', 'HB2', 'HB3', 'HB*']:

            if name in resonancesByName:
                HBs.append(resonancesByName[name])

        # Caros
        if ccpCodes.intersection(set(['Phe', 'Tyr', 'Trp', 'His'])):
        #if ('Phe' in ccpCodes or
        #    'Tyr' in ccpCodes or
        #    'Trp' in ccpCodes or
        #    'His' in ccpCodes):

            resonancesByAtomSiteName['Caro'] = []

            for name in ['CG', 'CD1', 'CD2', 'CE1',
                         'CE2', 'CE3', 'CZ', 'CZ2',
                         'CZ3', 'CH2', 'CE*', 'CD*']:

                if name in resonancesByName:

                    resonance = resonancesByName[name]

                    if 100 < resonance.CS < 170:

                        Caros.append(resonancesByName[name])

        # Calis
        else:

            for resonance in resonancesByAtomSiteName.get('C', []):

                if 0 < resonance.CS < 100:
                    # if resonance.atomName != 'C' :
                    Calis.append(resonance)

        resonancesByAtomSiteName['HA'] = HAs
        resonancesByAtomSiteName['HB'] = HBs
        resonancesByAtomSiteName['CX'] = CXs
        resonancesByAtomSiteName['Cali'] = Calis
        resonancesByAtomSiteName['Caro'] = Caros
        resonancesByAtomSiteName['CO'] = COs

    def __reduce__(self):

        return (generalFactory, (type(self),), self.__getstate__())

    def __getstate__(self):

        return (self.DataModel, self.spinSystemNumber,
                self.isJoker, self.resonanceDict,
                self.allowedResidues, self.ccpCodes,
                self.solutions, self.userDefinedSolutions)

    def __setstate__(self, state):

        (self.DataModel, self.spinSystemNumber,
         self.isJoker, self.resonanceDict,
         self.allowedResidues, self.ccpCodes,
         self.solutions, self.userDefinedSolutions) = state

    def connectToProject(self, nmrProject):
        ''' Connect the spin system to the correct resonanceGroup
            in CCPN. Can be ran after unpickling.
        '''

        if self.isJoker:

            return

        self.ccpnResonanceGroup = nmrProject.findFirstResonanceGroup(
            serial=self.spinSystemNumber)

        if not self.ccpnResonanceGroup:

            print 'Error: could not find spin system %s' % str(self.spinSystemNumber)
            return

        for resonance in self.resonanceDict.values():

            resonance.connectToProject()

    def getDescription(self, noSerialWhenSeqCodeIsPresent=False):
        '''Return a string describing the spin system, for instance
          '{22} Asp 100'. Depends on assignment status.
          kwargs: noSerialWhenSeqCodeIsPresent: if set to True, the
                                                serial will not be
                                                part of the string.
        '''

        if self.isJoker:

            return 'Joker'

        info = []

        resonanceGroup = self.getCcpnResonanceGroup()

        if not resonanceGroup:

            return '{%s}' % self.spinSystemNumber

        if resonanceGroup.residue:

            serial = '{%s}-' % resonanceGroup.serial
            resCode = resonanceGroup.residue.seqCode
            ccpCode = resonanceGroup.residue.ccpCode
            if noSerialWhenSeqCodeIsPresent:

                serial = ''

            return '%s%s%s' % (serial, resCode, ccpCode)

        # SpinSystem has one or more tentative assignments. Got this piece out
        # of EditSpinSystem.py in popups.
        elif resonanceGroup.residueProbs:

            ccpCodes = []
            seqCodes = []

            spinSystemInfo = ''

            for residueProb in resonanceGroup.residueProbs:
                if not residueProb.weight:
                    continue

                residue = residueProb.possibility

                if residue:

                    spinSystemInfo += '%s%s?/' % (residue.seqCode,
                                                  residue.ccpCode)

            if spinSystemInfo:

                spinSystemInfo = spinSystemInfo[:-1]

                if noSerialWhenSeqCodeIsPresent:

                    return spinSystemInfo

            return '{%s}-%s' % (resonanceGroup.serial, spinSystemInfo)

        elif resonanceGroup.ccpCode:

            return '{%s}-%s' % (resonanceGroup.serial, resonanceGroup.ccpCode)

        elif resonanceGroup.residueTypeProbs:

            typeProbCcpCodes = [residueTypeProb.possibility.ccpCode for residueTypeProb in resonanceGroup.residueTypeProbs]
            ccpCodeString = (('%s/' * len(typeProbCcpCodes)) %tuple(typeProbCcpCodes))[:-1]

            return '{%s}-%s' % (resonanceGroup.serial, ccpCodeString)

        else:

            return '{%s}' % resonanceGroup.serial

    def getIsJoker(self):
        '''Returns Bool indicating whether spin system is a joker.'''

        return self.isJoker

    def getSerial(self):
        '''Returns serial. Serial of spin system is equal
           to serial of resonanceGroup in CCPN.
        '''

        return self.spinSystemNumber

    def getCcpnResonanceGroup(self):
        '''Returns corresponding resonanceGroup. This is determined
           every time from scratch because during the assignment process
           the resonanceGroup might have been removed.

        '''

        if self.isJoker:

            return None

        ccpnResonanceGroup = self.DataModel.nmrProject.findFirstResonanceGroup(serial=self.spinSystemNumber)

        if ccpnResonanceGroup:

            return ccpnResonanceGroup

        else:

            print 'Error: could not find spin system %s , it possibly does not exist anylonger because it was deleted or merged with another spin-system. I advise to re-run the algorithm because the state of the project has changed.' % str(self.spinSystemNumber)
            return None
