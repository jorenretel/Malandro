
'''This is a macro for ccpn analysis. The only thing this macro does
   is compile c extentions necesarry for Malandro.
   If you downloaded a 'pre-compiled' version of Malandro you don't
   have to run this macro.
   If you do want to compile the c code (again) you can use this macro.
   The reason this is written in the form of a ccpn macro is that a lot
   of people run a 'pre-compiled' vesrion of CCPN Analysis on their
   computer. This version of Analysis comes with its own python
   interpreter, which might be different than the default one on your
   system. This macro is just an easy way to make sure the code gets
   compiled in a compatible way with the python version your Analysis
   is using.

   To run this macro: open CCPN Analysis, open a project and go to:

       Macro -> Organize Macros

   click on 'Add Macro', navigate to the location of this file and
   selected it. In the bottom half of the dialog, select 'runInstaller',
   and press the 'Load Macro' button. Now the macro should appear in
   the list of macros and you can select and run it. Most probably
   there will be some compiler messages and if the compilation process
   finishes successfully a message will appear:

   'Compiling is done, you can now run the macro.'

   For more installation instruction read the README. 

'''
import sys
import os

def runInstaller(argServer) :

    pathToPython = sys.executable
    print pathToPython
    workingDirectory = os.path.dirname(__file__)

    import subprocess
    print 'Compiling malandro for your version of ccpn analysis...'
    print 'using python version %s' %sys.version
    process = subprocess.Popen([pathToPython, 'setupC.py', 'build_ext', '--inplace'], cwd=workingDirectory, stdout=subprocess.PIPE)
    process.wait()
    if process.returncode == 0:
        print 'Compiling is done, you can now run the macro.'
