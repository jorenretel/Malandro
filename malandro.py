
import src.python.malandroGUI as GUI

from src.cython.malandro import *

def open_malandro(argServer):
  """Descrn: Opens the macro.
     Inputs: ArgumentServer
     Output: None
  """
  print 'version...'
  reload(GUI)
  program = GUI.connector(argServer.parent)