# Description
DDTREK is a script for automatic generation of PyMOL session, which will 
contain a reference structure, aligned structures 
and (optional) electron density maps for ligands.
# REQUIREMENTS:
- PyMOL version >= 2.5.0
- Input document:
    - reference structure (path to file after #REF keyword)
    - list of structures to be added to the session.
    See examples folder for details

# INSTALLATION:
- go to https://github.com/BiocrystLab/DDtrek and download the file `DDtrek_plugin.zip`
    Alternatively, use [direct link](https://raw.githubusercontent.com/BiocrystLab/DDtrek/main/ddtrek.zip)
- open PyMOL 2 (version 2.5 or higher)

 **Note**: DDtrek works in PyMOL 3 but less reliably.

- Click on the top menu "Plugin" and then on "Plugin manager".
- click on the tab "Install New Plugin", click on "Choose file..." and select the file `DDtrek_plugin.zip` that you downloaded earlier.


# TROUBLESHOOTING
- Most commonly occuring error is `Selector-Error: Invalid selection name "current_entry".` with PyMOL failing to create a new object. It is not exacly clear, what causes this error. Solution: close and open PyMOL, and rerun DDtrek
