
from src.python.malandroGUI import connector

def open_malandro(argServer):
  """Descrn: Opens the macro.
     Inputs: ArgumentServer
     Output: None
  """
  print 'version...'

  program = connector(argServer.parent)