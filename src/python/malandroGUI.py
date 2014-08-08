'''This module contains the toplevel of the graphical user
   interface of the malandro.
'''

def open_reference(argServer):
    """Descrn: Opens the macro.
       Inputs: ArgumentServer
       Output: None
    """

    print 'version...'

    Connector(argServer.parent)


import cPickle
import os
#from memops.general import Implementation
from memops.gui.ButtonList import UtilityButtonList
from memops.gui.Label import Label
from memops.gui.TabbedFrame import TabbedFrame
from memops.gui.MessageReporter import showWarning, showYesNo  # showMulti
from ccpnmr.analysis.popups.BasePopup import BasePopup
#from ccpnmr.analysis.macros.ArgumentServer import ArgumentServer

#from assignmentFunctions import assignSpinSystemstoResidues

from src.cython.malandro import Malandro, myDataModel
import src.python.assignMentTransferTab
reload(src.python.assignMentTransferTab)
from src.python.assignMentTransferTab import AssignMentTransferTab
import src.python.spectrumSelectionTab
reload(src.python.spectrumSelectionTab)
from src.python.spectrumSelectionTab import SpectrumSelectionTab
import src.python.saveLoadTab
reload(src.python.saveLoadTab)
from src.python.saveLoadTab import SaveLoadTab
import src.python.annealingSettingsTab
reload(src.python.annealingSettingsTab)
from src.python.annealingSettingsTab import AnnealingSettingsTab
import src.python.resultsTab
reload(src.python.resultsTab)
from src.python.resultsTab import ResultsTab
#import benchMark
# reload(benchMark)
#from benchMark import createBenchmark


class Connector(object):

    '''This is just a little class that contains all
       settings the user configures in the GUI.
       It acts as a connector between the real logic
       of the monte carlo procedure and the GUI.
    '''

    def __init__(self, guiParent):
        '''Init. args: guiParent'''

        self.chain = None
        self.useAssignments = True
        self.useTentative = True
        self.reTypeSpinSystems = False
        self.typeSpinSystems = True
        self.useDimensionalAssignments = True
        self.selectedSpectra = None
        self.amountOfSteps = 0
        self.amountOfRepeats = 0
        self.minTypeScore = 0.01
        self.leavePeaksOutFraction = 0.0
        self.shiftList = None
        self.nmrProject = None
        self.project = None
        self.minIsoFrac = None
        self.acceptanceConstantList = []
        self.sourceName = None
        self.ranAnnealling = False
        self.ranPreCalculations = False
        self.results = None
        self.GUI = ViewAssignmentPopup(guiParent, self)
        self.auto = Malandro()

        #self.addEnergyPoint = self.GUI.annealingSettingsTab.addEnergyPoint

        self.auto.registerTextObserver(self.updateInfoText)
        self.auto.registerEnergyObserver(
            self.GUI.annealingSettingsTab.addEnergyPoint)

    def update(self):
        '''Pulls all settings out of the GUI and gives
           them to the assignment algorithm written in Cython.
        '''

        self.chain = self.GUI.annealingSettingsTab.chain
        self.useAssignments = self.GUI.annealingSettingsTab.useAssignmentsCheck.get()
        self.useTentative = self.GUI.annealingSettingsTab.useTentativeCheck.get()
        self.reTypeSpinSystems = not self.GUI.annealingSettingsTab.useTypeCheck.get()
        self.typeSpinSystems = self.GUI.annealingSettingsTab.useAlsoUntypedSpinSystemsCheck.get()
        self.useDimensionalAssignments = self.GUI.annealingSettingsTab.useDimensionalAssignmentsCheck.get()
        self.amountOfRepeats = self.GUI.annealingSettingsTab.amountOfRepeats
        self.amountOfSteps = self.GUI.annealingSettingsTab.amountOfSteps
        self.minTypeScore = self.GUI.annealingSettingsTab.minTypeScore
        self.leavePeaksOutFraction = self.GUI.annealingSettingsTab.leavePeaksOutFraction
        self.shiftList = self.GUI.annealingSettingsTab.shiftList
        self.nmrProject = self.GUI.nmrProject
        self.project = self.GUI.project
        self.minIsoFrac = self.GUI.annealingSettingsTab.minIsoFrac

        #self.sourceName = self.GUI.sourceName

        self.selectedSpectra = []

        for spec in self.GUI.spectrumSelectionTab.specConfigList:

            if spec.used:

                self.selectedSpectra.append(spec)

        self.GUI.annealingSettingsTab.updateAcceptanceConstantList()
        self.acceptanceConstantList = self.GUI.annealingSettingsTab.acceptanceConstantList
        self.auto.updateSettings(self)

    def runAllCalculations(self):
        '''Run all calculations. This includes all
           calculations that have to be carried out
           before the annealing starts and the annealing
           itself.
        '''

        if not self.ranPreCalculations:
            self.update()
            self.GUI.annealingSettingsTab.disableIllegalButtonsAfterPrecalculations()
            self.preCalculateDataModel()
        self.startAnnealing()

    def preCalculateDataModel(self):
        '''Run all calculations that have to be carried
           out before the annealing can run. This is mainly
           setting up a data model, simulating spectra and
           matching those simulated spectra with the real
           spectra.
        '''

        print 'updating'

        self.update()

        if self.checkPoint1():

            print 'precalculate data model'
            self.auto.preCalculateDataModel()
            self.ranPreCalculations = True

    def startAnnealing(self):
        '''Run the annealing.'''

        if self.checkPoint2():

            self.auto.startMonteCarlo(amountOfRuns=self.amountOfRepeats,
                                      stepsPerTemperature=self.amountOfSteps,
                                      acceptanceConstants=self.acceptanceConstantList,
                                      fractionOfPeaksLeftOut=self.leavePeaksOutFraction)

            self.ranAnnealling = True

            self.getResults()

            # createBenchmark(self.results)

            # TODO: this is not very nice
            self.GUI.resultsTab.dataModel = self.results
            self.GUI.resultsTab.update()
            self.GUI.assignTab.update()

    def checkPoint1(self):
        '''First of two checkpoints, is run before the precalculations.
           Here the labelling scheme and the chain information is checked.
           Returns True if it is safe to proceed. Returns False otherwise.
        '''

        if self.checkSpectraAndLabelling() and self.checkChain():

            return True

        else:

            return False

    def checkPoint2(self):
        '''Second of two checkpoints, is run before the annealing starts.
           It is checked whether the precalculations already ran and whether
           the list with acceptance constants is valid.
           Returns True if it is safe to proceed. Returns False otherwise.
        '''

        if self.ranPreCalculations:

            if self.GUI.annealingSettingsTab.updateAcceptanceConstantList():
                return True
        else:

            string = '''Without running the pre-calculations,
                        the annealing can not run.
                     '''

            showWarning('Run pre-calculations first', string, parent=self.GUI)

        return False

    def checkSpectraAndLabelling(self):
        '''Check whether some spectra are selected by the user
           and whether a reference experiment has been set
           for this experiment (needed to simulate the spectra)
           and if the user said the labelling scheme should
           come automatically from the labelled sample,
           the labelled mixture should be set.
        '''

        if not self.selectedSpectra:

            string = 'Select one or more spectra.'

            showWarning('No Spectra Selected', string, parent=self.GUI)

            return False

        for spec in self.selectedSpectra:

            name = spec.ccpnSpectrum.experiment.name

            if not spec.ccpnSpectrum.experiment.refExperiment:

                string = name + (''' has no reference experiment, '''
                                 '''go to Experiment --> Experiments --> '''
                                 '''Experiment Types and select '''
                                 '''a experiment type.''')

                showWarning('No Reference Experiment',
                            string, parent=self.GUI)

                return False

            mixtures = spec.ccpnSpectrum.experiment.labeledMixtures

            if spec.labellingScheme is True and not mixtures:

                string = name + (''' is not connected to any labelled'''
                                 ''' sample.In order to determine the'''
                                 ''' labelling scheme automatically'''
                                 ''' from the experiment, this has to'''
                                 ''' be set. Go to'''
                                 ''' Menu-Modecule-Isotope Labelling'''
                                 ''' to do so.''')

                showWarning('No Labelled Sample', string, parent=self.GUI)

                return False

            elif len(mixtures) > 1:

                string = name + (''' is connected to more than one labelled'''
                                 ''' sample. This is physically impossible.'''
                                 ''' Go to Menu-Modecule-Isotope Labelling'''
                                 ''' if you want to change this.''')

                showWarning('Multiple Labelled Samples for Spectrum',
                            string, parent=self.GUI)

                return False

        return True

    def checkChain(self):
        '''Checks whether a chain is configured
           and whether this chain has residues.
        '''

        if not self.chain:

            string = 'Setup a molecular chain in Molecule --> Molecules'
            showWarning('No chain selected', string, parent=self.GUI)
            return False

        if not self.chain.residues:

            string = '''The selected chain does not have residues,
                        set up the residues in
                        Molecule --> Molecules --> Sequences
                     '''
            showWarning('No Residues', string, parent=self.GUI)
            return False

        return True

    def getResults(self):
        '''Gets the results from the assignment algorithm.
        '''

        if self.ranAnnealling:

            self.results = self.auto.getResults()

    def updateInfoText(self, text):
        '''Update the string with information that is
           displayed to the used in the GUI and the terminal.
        '''

        print text

        self.GUI.updateInfoText(text)

    def saveDataToPyc(self, fileName):
        '''Pickle the data model and save it to file to
           be use later. Since this storage method relies
           on pickle, it is not guaranteed that the data
           can be loaded in future versions of this macro.
        '''

        if self.results:

            if not fileName[-4:] == '.pyc':

                fileName = fileName + '.pyc'

            if os.path.exists(fileName):
                if not showYesNo('File exists',
                                 'File "%s" exists, overwrite?' % fileName,
                                 parent=self.GUI):

                    return

            self.updateInfoText('Saving results...')
            writeFile = open(fileName, 'wb')
            cPickle.dump(self.results, writeFile, protocol=2)
            writeFile.close()
            self.updateInfoText('Results succesfully saved.')

        else:

            string = 'There are no results to save at the moment.'
            showWarning('No results to save', string, parent=self.GUI)

    def loadDataFromPyc(self, fileName):
        '''Load previously pickled data.
               args: fileName: name of the file
        '''

        self.updateInfoText('Loading file...')

        # TODO : fix this, getting these values from the GUI is ridiculous
        self.project = self.GUI.project
        self.nmrProject = self.GUI.nmrProject

        results = None

        fileExists = False

        if os.path.exists(fileName):

            fileExists = True

        elif os.path.exists(fileName + '.pyc'):

            fileName = fileName + '.pyc'

            fileExists = True

        if fileExists:

            try:
                with open(fileName, 'rb') as readFile:

                    results = cPickle.load(readFile)

            except:

                string = ('''Something went wrong when trying to '''
                          '''load the previously produced results. '''
                          '''Probably the file you try to '''
                          '''open is not of the correct type.''')

                print string

                showWarning('Can not open file', string, parent=self.GUI)
                raise

            if results:

                if isinstance(results, myDataModel):

                    self.results = results
                    self.updateInfoText('File loaded succesfully.')
                    self.updateInfoText(('''Re-connecting to object in '''
                                         '''the ccpn analysis project...'''))
                    self.results.connectToProject(self.project,
                                                  self.nmrProject)
                    # TODO: this is not very nice
                    self.GUI.resultsTab.dataModel = self.results
                    self.updateInfoText('Done')
                    return

                else:

                    string = (('''You are trying to open a file that '''
                               '''was not created using this macro.'''))


                    showWarning('Can not open file', string, parent=self.GUI)

        else:

            string = 'The file you selected does not exist.'

            showWarning('No Such File', string, parent=self.GUI)

        self.updateInfoText(' ')


class ViewAssignmentPopup(BasePopup):
    '''The main popup that is shown when the macro is loaded.'''

    def __init__(self, parent, controler, *args, **kw):

        self.font = 'Helvetica 10'
        self.sFont = 'Helvetica %d'
        self.project = parent.project
        self.guiParent = parent

        self.sourceName = 'RefDB'

        self.connector = controler
        self.tabbedFrame = None
        self.saveLoadTab = None
        self.assignTab = None
        self.spectrumSelectionTab = None
        self.resultsTab = None
        self.annealingSettingsTab = None
        self.infoLabel = None

        self.waiting = False

        BasePopup.__init__(self, parent, title="Assignment Suggestions", **kw)

    def open(self):
        '''Open the popup. A standard method.'''

        self.updateAfter()
        BasePopup.open(self)

    def body(self, guiFrame):
        '''This method describes the outline of the body of the
           application.
               args: guiFrame: frame the body should live in.
        '''

        self.geometry('900x700')

        guiFrame.grid_columnconfigure(0, weight=1)
        guiFrame.grid_rowconfigure(0, weight=1)

        tabbedFrame = TabbedFrame(guiFrame,
                                  options=['Spectra',
                                           'Annealing',
                                           'Results',
                                           'Bulk Transfer Assignments To Project',
                                           'Save and Load'],
                                  callback=self.toggleTab,
                                  grid=(0, 0))

        self.tabbedFrame = tabbedFrame

        autoFrame, NAFrame, resultsFrame, assignFrame, saveFrame = tabbedFrame.frames

        self.spectrumSelectionTab = SpectrumSelectionTab(self, autoFrame)
        self.assignTab = AssignMentTransferTab(self, assignFrame)
        self.saveLoadTab = SaveLoadTab(self, saveFrame)
        self.annealingSettingsTab = AnnealingSettingsTab(self, NAFrame)
        self.resultsTab = ResultsTab(self, resultsFrame)

        bottomButtons = UtilityButtonList(tabbedFrame.sideFrame,
                                          helpUrl='www.google.com')
        bottomButtons.grid(row=0, column=0, sticky='e')

        self.infoLabel = Label(guiFrame, text=' ')

        self.infoLabel.grid(row=2, column=0, sticky='w')

        self.waiting = False
        self.updateAfter()

        self.administerNotifiers(self.registerNotify)

    def toggleTab(self, index):
        '''Toggle to a specific tab.
               args: index: index of tab.
        '''

        pass

    def administerNotifiers(self, notifyFunc):
        '''Subscribing to get messages when parts of
           the data model of the project changes.
        '''

        notifyFunc(self.spectrumSelectionTab.update,
                   'ccp.nmr.Nmr.Experiment', 'setName')
        notifyFunc(self.spectrumSelectionTab.update,
                   'ccp.nmr.Nmr.DataSource', 'setName')
        notifyFunc(self.spectrumSelectionTab.update,
                   'ccp.nmr.Nmr.DataSource', 'delete')
        notifyFunc(self.spectrumSelectionTab.update,
                   'ccp.nmr.Nmr.DataSource', '__init__')

        notifyFunc(self.annealingSettingsTab.updateChains,
                   'ccp.molecule.MolSystem.Chain', 'delete')
        notifyFunc(self.annealingSettingsTab.updateChains,
                   'ccp.molecule.MolSystem.Chain', '__init__')

    def updateAfter(self, *opt):
        '''update views, but waits untill all functions have returned.'''

        if self.waiting:
            return
        else:
            self.waiting = True
            self.after_idle(self.update)

    def update(self):
        '''update'''

        self.toggleTab(self.tabbedFrame.selected)
        self.waiting = False

    def destroy(self):
        '''Destroy the popup.'''

        self.administerNotifiers(self.unregisterNotify)
        BasePopup.destroy(self)

    def updateInfoText(self, string):
        '''Update the text shown at the bottom of the GUI.
               args: string: the string displayed to the user.
        '''

        self.infoLabel.set(string)
        self.infoLabel.update()
