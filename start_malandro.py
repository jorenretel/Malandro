'''Just a startup script.'''

import src.python.malandroGUI as GUI


def open_malandro(argServer):
    """Descrn: Opens the macro.
       Inputs: ArgumentServer
       Output: None
    """
    print 'version...'
    reload(GUI)
    GUI.Connector(argServer.parent)
