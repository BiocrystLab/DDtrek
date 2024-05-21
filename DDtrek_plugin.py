import os
import random
from pymol import cmd
import gemmi


'''
# DDTREK is a script for automatic generation of PyMOL session, which will 
contain a reference structure, aligned structures 
and (optional) electron density maps for ligands.
REQUIREMENTS:
- PyMOL version >= 2.5.0
- Input document:
    - reference structure (path to file after #REF keyword)
    - list of structures to be added to the session.
    See ddtrek.inp for example

INSTALLATION:
- open PyMOL
- install Gemmi module from PyMOL command line:
    import pip
    pip.main(['install','gemmi'])
- use PyMOL build-in plugin manager to install DDtrek_plugin


Algorithm of DDtrek:
1) get list of structures, maps, chains and ligand names from txt file
2) iterate over structures: extract specified chains and aligns them against reference structure
3) for each ligand create density map object. 

TODO: 
- Draw maps on request by type. If not in the list - print list of maps, delete output folder and exit. Right now the script generates meshes only for 2FOFC maps.
- automate map drawing by detecting required map type
- Automate gemmi installation: import gemmi and run pip install gemmi in case of ImportError
- Implement cut-off distances for coordinates and maps as constants
- Add GUI
- write temporary files ligand.pdb and masked.ccp4 to the folder specified by TMP environment variable

'''

DEBUG = True # True will produce extra output text in terminal


def load_mtz_map_fragment(mtzfilename:str, mapobjname:str, margin=3) -> None:
    '''
    DESCRIPTION
    :mtzfilename: path to mtz map file
    :mapobjname: arbitrary name of map entry
    :margin: map cutoff distance around the ligand in angstroms
    The subroutine extracts map fragment around ligand atoms using the following protocol
    NOTE: requires rw access to the folder with ddtrek file

    Algorithm:
    1. Generate 2Fo-Fc map from the specified MTZ file :mtzfilename:
    2. Save  fragment of the map around the ligand into temporary file
    3. load map into session
    '''
    mtz = gemmi.read_mtz_file(mtzfilename)
    m = gemmi.Ccp4Map()
    # ATM only 2Fo-Fc maps with column names 2FOFCWT or FTW(Refmac5)  are recognized
    try:
        m.grid = mtz.transform_f_phi_to_map('2FOFCWT','PH2FOFCWT', sample_rate=3)# column labels for mtz with map coefficients in Phenix and Buster
    except:
        m.grid = mtz.transform_f_phi_to_map('FWT','PHWT', sample_rate=3)#refmac5 default names for map coefficients
    m.update_ccp4_header()
    ligand = gemmi.read_structure('ligand.pdb') # if map file was specified, this file is generated in main body of DDtrek 
    # command below seems to generate Gb files of map fragments and overflows memory in PyMOL < 2.5
    m.set_extent(ligand.calculate_fractional_box(margin=margin)) #cut map fragment with margin around the ligand
    m.write_ccp4_map('masked.ccp4')
    cmd.load('masked.ccp4', mapobjname)
    # remove temporary files
    os.remove('masked.ccp4')
    os.remove('ligand.pdb')
    return

def ddtrek(input_filename: str) -> None:
    '''
    DESCRIPTION

    Aligns structures from input file to the reference structure. Optionally generates electron density maps
    
    See publication: DDtrek (2024)

    USAGE

    ddtrek [ path to input file ]

    where
    :input_filename: path (absolute or relative) to DDtrek.in file
    '''

    #Change working directory to location of input file
    abspath = os.path.abspath(input_filename)
    dname = os.path.dirname(abspath)
    fname = os.path.basename(input_filename)
    os.chdir(dname)
    assert (os.access(dname, os.W_OK)), "Folder with DDtrek input file should be writable"
    ## get the list of preloaded structures in currently opened PyMOL window
    preloaded_structures = cmd.get_object_list()

    # USER input file with user-specified structures and ligand selectors
    pdb_mtz_list = open(fname,'r').readlines()

    ### Iterate over entries in input list and align them against reference structure

    for entry in pdb_mtz_list:
        '''
        For every line in the input file:
        0. process control cards (line starts with #)
            or 
           read line with pdb, [optional mtz with map coefficients], chain, resi, entry name, [ optional alignment protein chain]. 
            Where:
            - Chain, resi define ligand in pdb
            - optional chain control which protein chain from extracted structure is aligned to reference chain
        1.load pdb
        2. extract into new object: ligand specified by chain, resi + its close environment: protein chains + symmetry related entries + water molecules 
        3. align new object with the reference structure 
        4. group new object into the group defined by #G card
        5. Optinally, extract map fragment and generate mesh map
        '''
        #PATH to currently processed pdb and mtz files are stored here
        pdb, mtz = None,None
        # optional  protein chain for alignment (useful when ASU contains more than one chain and alignment is not working as intended)
        align_chain = None

        # PARSE INPUT LINE
        ## Remove trailing whitespaces
        entry = entry.rstrip()
        
        # mandatory param line with #G card serves as a prefix for entries. 
        # Convenient for grouping in PyMOL,e.g.: `group pref,pref_*` will make `pref` group for entries following `#G pref` line
        if entry.startswith('#G'): # for example #G Glu32
            group_name = entry.split()[1]
            continue # go to next line in input file
        # Load reference structure after #REF key if reference is missing from PyMOL session
        if entry.startswith('#REF') and ('reference' not in preloaded_structures):
            ref_pdb = entry.split()[1]
            cmd.load(ref_pdb,'tmp_reference')
            # if reference PDB has >1 chains or chain is specified in DDtrek.in
            # identify a single protein chain for use as a reference for alignment:
            if len(entry.split()) == 3:
                reference_chain = entry.split()[2]
            else:
                reference_chain = cmd.get_chains('tmp_reference and polymer')[0]
            cmd.create('reference',f'tmp_reference and polymer and chain {reference_chain}')
            cmd.delete('tmp_reference')
            continue

        if entry.startswith('# ') or len(entry) == 0: # skip comments and empty lines
            continue
        
        # Check number of elements per line. Mtz and alignment chains are optional
        # pdb, (optional) mtz, ligand chain, ligand resi, entry_name, (optional) alignment chain
        # chain and resi define ligand

        # Process correct lines
        #split every line into corresponding elements
        
        pdb = entry.split()[0]
        # discard non-PDB lines and malformed lines
        if not pdb.lower().endswith('pdb') or (pdb.lower().endswith('pdb') and len(entry.split())) < 4:
            print('Malformed line with PDB %s . Skipping...' % pdb )
            continue

        # load map mtz and check extension _ only mtz maps are allowed
        # TODO: may be rewrite this section. It looks terrible BUT it works...
        mtz = entry.split()[1]
        if not mtz.lower().endswith('mtz'):
            # if second element is not mtz
            mtz = None
            # finally, parse entry name,  selectors for the ligand (resi, chain) and optional alignment chain
            # It was easier to copy paste rather than create a function...
            # looks bulky but easier to read
            if len(entry.split()) == 4:
                entry_name = entry.split()[-1] # last element should be the entry name
                resi = entry.split()[-2]
                chain = entry.split()[-3]
            elif len(entry.split()) == 5: # if optional protein chain was specified
                align_chain = entry.split()[-1]
                entry_name = entry.split()[-2] # last element should be entry name
                resi = entry.split()[-3]
                chain = entry.split()[-4]
        else: # else if mtz map is specified, consider cases with and w/o alignment chain
            if len(entry.split()) == 5:
                entry_name = entry.split()[-1] # last element should be entry name
                resi = entry.split()[-2]
                chain = entry.split()[-3]
            elif len(entry.split()) == 6: # selection fo the case with optional protein chain
                align_chain = entry.split()[-1]
                entry_name = entry.split()[-2] # last element should be entry name
                resi = entry.split()[-3]
                chain = entry.split()[-4]

 

        # Skip preloaded structure by checking name of the current entry against all entries loaded into PyMOL session
        if (any([i.find(entry_name)>=0 for i in preloaded_structures])): # be careful with 'any'
            print('Entry %s is already loaded! Skipping...' % (entry_name))
            continue

        #GENERATE ENTRY FROM PDB STRUCTURE
        ### ALL MAGIC is HERE
        # load coordinates 
        cmd.load(pdb,'current_entry')
        
        # Define PyMOL formatted ligand selection
        ligand_residue = '(current_entry and chain %s and resi %s)' % ( chain, resi)

        # If mtz is specified, save ligand coordinates for map fragment truncation
        if mtz:
            cmd.save('ligand.pdb', '%s' % (ligand_residue), format='pdb')
            if DEBUG:
                print("Extracting ligand:",ligand_residue)

        # CREATE ENTRY in PYMOL

        ### Generate symmetry mates 
        ### expand selection to nearby protein chains - useful when binding pocket is near crystallographic interface
        ### Protein chains will have unique segid
        cmd.symexp(prefix='symmetry',object='current_entry',selection=ligand_residue, cutoff=4.5, segi=1)
        # remove  atoms in symmetry mates beyond 6 A from the ligand in order to avoid unnecessary calculations
        cmd.remove("%s beyond 6 of %s" % ('symmetry*',ligand_residue))

        
        # Now select the chain that contains the ligand, ligand itself and everything in 4.5 A radius around the ligand of interest. 
        # Select ligand, all closely situated protein chains and water molecules
        extracted_selection = "current_entry and (%s or ((bychain %s around 4.5) and polymer) or (%s around 4.5 and resn HOH))" % (ligand_residue,ligand_residue,ligand_residue)
        # now add symmetry mates to it
        #extracted_selection = "(%s) or symmetry*" % tmp_selection


        if DEBUG:
            print('selection:%s' % extracted_selection)

        #Extract coordinates into the new entry
        cmd.create(entry_name, selection=extracted_selection,
                    source_state=0, target_state=0,discrete=0)


        #ALIGN newly created structure to reference using all CA atoms
        if align_chain:
            alignment_selection = "%s and polymer and chain %s and name CA" % (entry_name, align_chain)
        else:
            alignment_selection = "%s and polymer and name CA" % (entry_name)
        if DEBUG:
            print('Alignment selection:%s\n' % alignment_selection)
        cmd.align(alignment_selection, 'reference and name CA') # align using C-alpha atoms to reference

        #Cleanup. Remove loaded PDB and generated symmetry mates
        cmd.delete('current_entry') # remove original PDB
        cmd.delete('symmetry*') # remove symmetry mates
        
        #Truncate atoms in new entry. Limit by 4.5 A radius around the ligand
        pocket_selection = "(byres(%s within 4.5 of (%s and chain %s and resi %s)))" % (entry_name,entry_name,chain, resi)
        cmd.remove("%s and not %s" % (entry_name, pocket_selection))

        # STYLING PDB
        # Change colors by random and represenation of the protein
        # select a random color from pymol set of colors and color carbons only
        structure_color = random.choice(cmd.get_color_indices())[0]
        #TODO fix string formatting to match single stype
        print("Structure " + entry_name + " is colored by " + structure_color)
        cmd.color(structure_color,entry_name + " and elem C")
        cmd.hide('everything', entry_name)
        cmd.show('lines', "%s and polymer" % entry_name)
        cmd.show('sticks', "%s and hetatm" % entry_name)# likely the ligand is the only hetatm so rather generous selection
        cmd.show('nb_spheres', "%s and resn hoh" % entry_name) #waters by spheres

        # Finally, add new entry to the group
        cmd.group(group_name, entry_name, 'add')

        #END PDB manipulation
        
        #MAPS MESH GENERATION

        # check if correct filepath loaded into mtz variable. Otherwise, skip map generation and go to next entry
        if mtz is None:
            continue
        mtz_name = "%s_map" % (entry_name) # like Z111112321_mtz
        ligand_selection = '%s and chain %s and resi %s' % (entry_name, chain, resi) # ligand selection in newly created entry
        if DEBUG:
            print("Map generation name:%s\n" % ligand_selection)
        map_name = entry_name + '_mesh' 
        load_mtz_map_fragment(mtz, mtz_name)
        cmd.matrix_copy(entry_name, mtz_name)
        # pymol isomesh generation map at 1 sigma around ligand in 1.8 angstrom radius
        cmd.isomesh(map_name, mtz_name, 1, ligand_selection, carve=1.8)

        # MAP REPRESENTATION 
        # color will match color of carbon atoms
        cmd.color(structure_color, map_name)

        #GROUPING: 
        # group add all related entries to the group as specified by #G
        cmd.group(group_name, map_name, 'add')
    
    # At the end of run put mtz files with map fragment into separate folder
    cmd.group ('map_objects',"*map", 'add')
    # Finally:   
    # Draw surface for the reference and make it transparent
    cmd.show('surface', 'reference')
    cmd.color('gray60', 'reference')
    cmd.set('transparency', '0.5')
    cmd.refresh()
    # disable mtz group — we don't need to show it anyway
    cmd.disable('map_objects')
    cmd.zoom()
    #cmd.save('session.pse')

cmd.extend("ddtrek", ddtrek)
