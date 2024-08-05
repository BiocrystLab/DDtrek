'''
PyMOL DDtrek plugin. 

The plugin creates a session with protein-ligand structures aligned and truncated around point of interest

Created from PyMOL Qt example plugin: https://github.com/Pymol-Scripts/pymol2-demo-plugin
'''

from __future__ import absolute_import
from __future__ import print_function
import os, sys


# import fails without adjustments to sys.path
plugin_dir = os.path.dirname(__file__)
sys.path.insert(0,  plugin_dir)

from ddtrek import ddtrek


def __init_plugin__(app=None):
    '''
    Add an entry to the PyMOL "Plugin" menu
    '''
    from pymol.plugins import addmenuitemqt
    addmenuitemqt('DDtrek Plugin', run_plugin_gui)


# global reference to avoid garbage collection of our dialog
dialog = None


def run_plugin_gui():
    '''
    Open our custom dialog
    '''
    global dialog

    if dialog is None:
        dialog = make_dialog()

    dialog.show()


def make_dialog():
    # import required libraries
    from pymol import cmd
    from pymol.Qt import QtWidgets
    from pymol.Qt.utils import loadUi

    # Line copied from https://github.com/BobSchiffrin/PyXlinkViewer
    getOpenFileNames = QtWidgets.QFileDialog.getOpenFileNames

    # create a new Window
    dialog = QtWidgets.QDialog()
    


    # populate the Window from our *.ui file which was created with the Qt Designer
    uifile = os.path.join(os.path.dirname(__file__), 'ddtrek.ui')
    form = loadUi(uifile, dialog)


    def open_file():
        """
        Callback for open file function
        """
        curr_dir = os.getcwd()

        try:
            fname = getOpenFileNames(dialog, 'Open DDtrek input file...', curr_dir)[0][0]
        except IndexError:
            # to supress index error if no file was selected
            fname = ''

        if fname:
            form.input_filename.setText(fname)


    def run():
        '''
        TODO:run ddtrek here
        '''
        fname = form.input_filename.text()
        if fname:
            print(f'Opening {fname}...')
            ddtrek(fname)

    # hook up button callbacks
    form.button_open.clicked.connect(run)
    form.button_browse.clicked.connect(open_file)
    form.button_close.clicked.connect(dialog.close)

    return dialog
