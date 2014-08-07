'''This module contains the class describing
   the tab where reults can be loaded and saved.
'''

from memops.gui.FileSelect import FileSelect, FileType
from memops.gui.ButtonList import ButtonList
from memops.gui.MessageReporter import showWarning


class SaveLoadTab(object):
    '''The tab that lets the user FileSelect
       a file to open or to save.
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
        self.loadAndSaveButtons = None
        self.fileselectionBox = None
        self.body()

    def body(self):
        '''Setting up the body of this view.'''

        frame = self.frame

        file_types = [FileType('pyc', ['*.pyc'])]
        self.fileselectionBox = FileSelect(
            frame, multiSelect=False, file_types=file_types)
        self.fileselectionBox.grid(
            row=0, column=0, columnspan=6, sticky='nsew')

        texts = ['Load', 'Save']
        commands = [self.loadDataFromPyc, self.saveDataToPyc]
        self.loadAndSaveButtons = ButtonList(
            frame, commands=commands, texts=texts)
        self.loadAndSaveButtons.grid(row=1, column=0, sticky='nsew')

        frame.grid_rowconfigure(0, weight=1)
        frame.grid_columnconfigure(0, weight=1)

    def saveDataToPyc(self):
        '''Save the data to the selected file.'''

        if self.guiParent.connector.results:

            fileName = self.fileselectionBox.getFile()
            self.guiParent.connector.saveDataToPyc(fileName)
            self.fileselectionBox.updateFileList()

        else:

            string = ('''There are no results to save, '''
                      '''the algorithm has not run yet.''')

            showWarning('No Results to Save', string, parent=self.guiParent)

    def loadDataFromPyc(self):
        '''Load the data from the selected file.'''

        fileName = self.fileselectionBox.getFile()
        self.guiParent.connector.loadDataFromPyc(fileName)
        self.guiParent.resultsTab.update()
