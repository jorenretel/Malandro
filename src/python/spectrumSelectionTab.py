'''The tab containing a big table with all the spectra.
   The user can select which spectra from the project to
   use with the assignment algorithm. For each spectrum
   a peak list and a labelling scheme can be selected.
'''

from memops.gui.PulldownList import PulldownList
from memops.gui.ScrolledMatrix import ScrolledMatrix


class SpectrumSelectionTab(object):
    '''Class describing the tab containing a big table with
       all the spectra. The user can select which spectra from
       the project to use with the assignment algorithm. For
       each spectrum a peak list and a labelling scheme can be
       selected.
    '''

    def __init__(self, parent, frame):
        '''Init. args: parent: the guiElement that this
                               tab is part of.
                       frame:  the frame this part of the
                               GUI lives in.
        '''

        self.guiParent = parent
        self.frame = frame
        self.project = parent.project
        self.nmrProject = parent.nmrProject
        self.specConfigList = []

        self.displayTable = None
        self.selectedAutoSpec = None
        self.autoLabellingPulldown = None
        self.autoPeakListPulldown = None

        self.body()
        self.updateSpecSelection()
        # self.setupSpectrumSettings()
        self.updateAutoMatrix()

    def body(self):
        '''Sets up the body of this view.'''

        frame = self.frame
        frame.grid_rowconfigure(2, weight=1)
        frame.grid_columnconfigure(0, weight=1)

        headingList = [
            '#', 'Spectrum', 'Peak List', 'use?', 'labelling scheme']

        tipTexts = ['Row number', 'spectrum name',
                    'which peak list to use',
                    'use this spectrum?',
                    'Which labelling scheme belongs to this spectrum?']

        self.autoLabellingPulldown = PulldownList(
            self.guiParent, self.setAutoLabellingScheme)
        self.autoPeakListPulldown = PulldownList(
            self.guiParent, self.setAutoPeakList)

        editWidgets = [None, None, self.autoPeakListPulldown,
                       None, self.autoLabellingPulldown]

        editGetCallbacks = [None, None, self.getAutoPeakLists,
                            self.changeUse, self.getAutoLabellingSchemes]

        editSetCallbacks = [None, None, None, None, None]

        self.displayTable = ScrolledMatrix(frame, headingList=headingList,
                                           callback=self.selectAutoSpec,
                                           editWidgets=editWidgets,
                                           multiSelect=False,
                                           editGetCallbacks=editGetCallbacks,
                                           editSetCallbacks=editSetCallbacks,
                                           tipTexts=tipTexts)
        self.displayTable.grid(row=2, column=0, sticky='nsew')

    def updateSpecSelection(self):
        '''Update the list with spectrum settings.
           This function can be called to initialize
           this list or to update it when some spectra
           are added in the project.
        '''

        # TODO: also implement what has to happen
        # when spectra get deleted from the project.

        presentSpectrumSettings = set([specSetting.ccpnSpectrum for specSetting in self.specConfigList])

        for expt in self.nmrProject.sortedExperiments():
            for spec in expt.sortedDataSources():
                if spec.dataType == 'processed' and spec not in presentSpectrumSettings:

                    newSpectrumSetting = SpectrumSettings()
                    newSpectrumSetting.ccpnSpectrum = spec
                    newSpectrumSetting.peakList = spec.getActivePeakList()
                    self.specConfigList.append(newSpectrumSetting)

    def getLabellingSchemes(self):
        '''Returns all labellinSchemes in the project'''

        return [True, None, ] + self.project.sortedLabelingSchemes()

    def setAutoLabellingScheme(self, scheme):
        '''Select a labelling scheme for the
           selected spectrum.
        '''

        self.selectedAutoSpec.setupLabellingScheme(scheme)
        self.updateAutoMatrix()

    def getAutoLabellingSchemes(self, notifyScheme=None):
        '''Sets up the the pulldown lists with labelling schemes
           used for the pulldown lists.
        '''
        #names = []
        #index = 0
        #
        #scheme = self.selectedAutoSpec.labellingScheme
        #schemes = self.getLabellingSchemes()
        #
        #
        # if schemes:
        #  names = ['Automatic from sample','<None>'] + [sc.name for sc in schemes[2:]]
        #
        #  index = schemes.index(scheme)

        names = ['Automatic from sample']
        schemes = [True]
        index = 0

        self.autoLabellingPulldown.setup(names, schemes, index)

    def setAutoPeakList(self, peakList):
        '''Select a peakList for the selected spectrum.'''

        self.selectedAutoSpec.setupPeakList(peakList)

        self.updateAutoMatrix()

    def selectAutoSpec(self, obj, *args, **kwargs):
        '''Set the selection of a spectrum in the table.
           args:
               obj: (SpectrumSettings) the spectrum that
                    should be selected.
        '''

        self.selectedAutoSpec = obj

    def getAutoPeakLists(self, *args, **kwargs):
        '''Set up the the autoPeakListPulldown for the
           selected spectrum.
        '''

        names = []
        index = 0

        peakList = self.selectedAutoSpec.peakList

        peakLists = self.selectedAutoSpec.ccpnSpectrum.sortedPeakLists()

        if peakLists:
            names = [str(pl.serial) for pl in peakLists]

            index = peakLists.index(peakList)

        self.autoPeakListPulldown.setup(names, peakLists, index)

    def changeUse(self, *args):
        '''Toggle whether a spectrum is used in the
           procedure or not.
        '''

        if self.selectedAutoSpec.used is False:
            self.selectedAutoSpec.changeSpectrumUse(True)
        elif self.selectedAutoSpec.used is True:
            self.selectedAutoSpec.changeSpectrumUse(False)

        self.updateAutoMatrix()

    def updateAutoMatrix(self):
        '''Update the whole table.'''

        textMatrix = []
        colorMatrix = []

        spectra = self.specConfigList

        for i, spectrum in enumerate(spectra):

            dataSource = spectrum.ccpnSpectrum
            expt = dataSource.experiment

            name = '%s:%s' % (expt.name, dataSource.name)

            peakListName = str(spectrum.peakList.serial)

            if spectrum.labellingScheme is None:
                schemeName = 'None'

            elif spectrum.labellingScheme is True:

                schemeName = 'Automatic from sample'

            else:

                schemeName = spectrum.labellingScheme.name

            datum = [i + 1, name, peakListName,
                     spectrum.used and 'Yes' or 'No', schemeName]

            textMatrix.append(datum)

            hexColors = dataSource.analysisSpectrum.posColors
            hexColor = hexColors[int(0.7 * len(hexColors))]

            if spectrum.used:

                colorMatrix.append([hexColor,
                                    hexColor,
                                    hexColor,
                                    hexColor,
                                    hexColor])

            else:

                colorMatrix.append([hexColor, None, None, None, None])

        self.displayTable.update(textMatrix=textMatrix,
                                 objectList=spectra,
                                 colorMatrix=colorMatrix)

    def update(self, *args):
        '''update the view'''

        self.updateSpecSelection()
        self.updateAutoMatrix()


class SpectrumSettings(object):
    '''A small class to temporarely store some
       info about spectra the user configures in the GUI.
       Attributes:
           ccpnSpectrum: (ccp.api.nmr.Nmr.DataSource) the spectrum
                         these configurations belong to.
           labellingScheme: if True, labelling scheme will be
                            automatically fetched from the labelled sample.
                            Else labelling information is not used and
                            every colabelling fraction is 1 by default.
                            Choosing a specific labelling scheme is not
                            possible any longer and is also not
                            nescesarry because a labelled sampe can be
                            configured very easily.
           peakList: (ccp.api.nmr.Nmr.PeakList) that should be used
                     when evaluating this spectrum.
           used: (boolean) True if the spectrum is used by the algorithm.

    '''

    def __init__(self, ccpnSpectrum=None,
                 labellingScheme=True,
                 peakList=None, used=False):
        '''Init'''

        self.ccpnSpectrum = ccpnSpectrum
        self.labellingScheme = labellingScheme
        self.peakList = peakList
        self.used = used

    def setupLabellingScheme(self, labellingScheme):
        '''Set labelling scheme for this spectrum.
           args: labellingScheme: True or None, where
                 True corresponds to automatically fething
                 the labelling from the labelled sample
                 that the spectrum belongs to. None
                 means 'no labelling info'. In that case
                 all colabelling fraction default to 1.
        '''

        self.labellingScheme = labellingScheme

    def setupPeakList(self, peakList):
        '''Set th peakList for this spectrum.
           args: peaksList (ccp.api.nmr.Nmr.PeakList) peak list
                 that should be used when evaluating this specttrum.
        '''
        self.peakList = peakList

    def changeSpectrumUse(self, used):
        '''Set whether spectrum is used or not.
           args: used (boolean) True if the spectrum
                 should be taken into account, False if
                 not.
        '''

        self.used = used

        if not self.ccpnSpectrum.experiment:

            print ('''This spectrum has not been connected '''
                   '''to a experiment type. If you want to '''
                   '''use this spectrum, go to '''
                   '''experiments --> experiment and configure this.''')
