'''
PyMOL DDtrek plugin. 

The plugin creates a session with protein-ligand structures aligned and truncated around point of interest

Created from PyMOL Qt example plugin: https://github.com/Pymol-Scripts/pymol2-demo-plugin
'''

from __future__ import absolute_import
from __future__ import print_function
import functools
import os, sys


# import fails without adjustments to sys.path
plugin_dir = os.path.dirname(__file__)
sys.path.insert(0,  plugin_dir)

DEBUG = os.environ.get('DEBUG_DD',False)

from ddtrek import ddtrek, map_extract


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
    from pymol.Qt.utils import getSaveFileNameWithExt


    # create a new Window
    dialog = QtWidgets.QDialog()

    # save/open menu items copied from https://github.com/BobSchiffrin/PyXlinkViewer
    getOpenFileNames = QtWidgets.QFileDialog.getOpenFileNames
    


    # populate the Window from our *.ui file which was created with the Qt Designer
    uifile = os.path.join(os.path.dirname(__file__), 'ddtrek.ui')
    form = loadUi(uifile, dialog)


    def open_file(input_field, tooltip='Open'):
        """
        Callback for open file function
        input field like form.input_filename
        """
        curr_dir = os.getcwd()

        try:
            fname = getOpenFileNames(dialog, tooltip, curr_dir)[0][0]
        except IndexError:
            # to supress index error if no file was selected
            fname = ''

        if fname:
            input_field.setText(fname)

    def save_map(input_field):
        '''
        Specify location of (optional) output map
        '''
        save_str = getSaveFileNameWithExt(dialog,'Save as...')
        if save_str:
            input_field.setText(save_str)
#/Volumes/Data/bigDownloads/emd_3944.map.gz
    def run():
        '''
        run ddtrek using filename in text field
        '''
        fname = form.input_filename_2.text()
        if fname:
            print(f'Opening {fname}...')
            ddtrek(fname)
            dialog.close()

    def extract():
        '''
        extract section of map
        '''
        map_fname = form.input_filename_4.text()
        selection = form.input_selection_1.text()
        mapout_fname = form.input_filename_5.text()
        if map_fname and selection:
            map_extract(map_fname, selection,savedmap=mapout_fname)
            dialog.close()

    if DEBUG:
        form.input_filename_4.setText('/Volumes/Data/bigDownloads/emd_3944.map.gz')
        form.input_filename_5.setText('/Volumes/Data/bigDownloads/emd_3944_test.map')

        
    # hook up button callbacks for ddtrek tab
    form.button_open_2.clicked.connect(run)
    form.button_browse_2.clicked.connect(functools.partial(open_file,form.input_filename_2, 'Open DDtrek file'))
    form.button_close_2.clicked.connect(dialog.close)

    # hook up button callback for map extraction tab
    form.button_extract_1.clicked.connect(extract)
    form.button_browse_4.clicked.connect(functools.partial(open_file, form.input_filename_4, 'Open density map'))
    form.button_browse_5.clicked.connect(functools.partial(save_map, form.input_filename_5))
    form.button_close_4.clicked.connect(dialog.close)


    return dialog
