'''Just a startup script.'''

import malandro.gui.malandroGUI as GUI
from version import __version__


def open_malandro(argServer):
    """Descrn: Opens the macro.
       Inputs: ArgumentServer
       Output: None
    """
    print 'Malandro Version {}'.format(__version__)
    reload(GUI)
    GUI.Connector(argServer.parent)
