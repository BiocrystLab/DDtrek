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
`DDtrek` relies on user for correct formatting of input file. In case of trouble, check that input file formatted properly.


Some common errors:

- Error `Selector-Error: Invalid selection name "current_entry".` with PyMOL failing to create a new object. It is not exacly clear, what causes this error. **Solution**: close and open PyMOL, and rerun DDtrek
- Error `gemmi is not definded`: `DDtrek` attempts to install dependency automatically using `pip`. This might fail on Linux machines if user doesn't have access to PyMOL installation, e.g. when PyMOL installed using package manager. **Solution**: Ask system administrator to install `gemmi` package using `sudo pip install gemmi`
