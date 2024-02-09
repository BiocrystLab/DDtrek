# Description
DDTREK is a script for automatic generation of PyMOL session, which will 
contain a reference structure, aligned structures 
and (optional) electron density maps for ligands.
# REQUIREMENTS:
- PyMOL version >= 2.5.0
- Input document:
    - reference structure (path to file after #REF keyword)
    - list of structures to be added to the session.
    See ddtrek.inp for example

# INSTALLATION:
- open PyMOL
- install Gemmi module from PyMOL command line:
    ```import pip
    
    pip.main(['install','gemmi'])
    ```
- use PyMOL build-in plugin manager to install DDtrek_plugin
