import os
import random
import tempfile
from pymol import cmd

try:
    import gemmi
except ImportError:
    import pip
    print('Installing missing Gemmi library...')
    pip.main(['install','gemmi'])


'''
Biocrystallography, KU Leuven 2024

DDTREK is a script for generation of PyMOL session, which will 
contain a reference structure, aligned structures 
and (optional) electron density maps for the ligands.
REQUIREMENTS:
- PyMOL version >= 2.5.0
- Input document with:
    - reference structure (path to file after #REF keyword)
    - list of structures to be added to the session.
    See ddtrek.inp for example



Algorithm of DDtrek:
1) get list of structures, maps, chains and ligand names from txt file
2) iterate over structures: extract specified chains and aligns them against reference structure
3) for each ligand create density map object. 

TODO: 
- Draw maps on request by type. If not in the list - print list of maps, delete output folder and exit. Right now the script generates meshes only for 2FOFC maps.
- automate map drawing by detecting required map type
- Implement cut-off distances for coordinates and maps as constants
- write temporary files ligand.pdb and masked.ccp4 to the folder specified by TMP environment variable

'''

DEBUG = os.environ.get('DEBUG_DD',False)

def map_extract(mapobj, selection, margin=9, savedmap='') -> None:
    '''
    DESCRIPTION
    For a given :selection: and atoms within :margin: extracts its density map 
    from the file specified in :mapobj:
    
    Optinally, extracted map fragment could be saved in file specified by :savedmap:
    #TODO: add tab in UI to allow map extraction from there
    '''
    tmp_dir = tempfile.gettempdir()
    os.chdir(tmp_dir)
    cmd.save('ligand.pdb', f'{selection} or {selection} around {margin}', format='pdb')
    if savedmap:
        print(f'Saving {savedmap}')
        load_mtz_map_fragment(mapobj,'extracted_map',savedmap=savedmap)
    else:
        load_mtz_map_fragment(mapobj,'extracted_map')

def load_cryoem_map_fragment(mapfile:str, savedmap='',mapout='masked.ccp4', ligand='ligand.pdb') -> None:
    '''
    DESCRIPTION
    :mapfile: path to cryoem map file
    :mapout: path to output file
    :ligand: path to selection containing pdb
    
    '''
    # 0. Load map into CcpMap
    map_obj = gemmi.read_ccp4_map(mapfile,setup=True)
    mapgrid = map_obj.grid
    # 1. Define new unit cell confined within a box around the ligand
    
    lig_obj = gemmi.read_pdb(ligand)
    box = lig_obj.calculate_box(margin=3)#box with + 3A from all atoms?
    boxsize = box.get_size()
    box_xdim = boxsize.x
    box_ydim = boxsize.y
    box_zdim = boxsize.z
    gridspacing = mapgrid.spacing
    start = box.minimum
    startx = start.x
    starty = start.y
    startz = start.z
    startpoint = [int(startx / gridspacing[0]), int(starty / gridspacing[1]), int(startz / gridspacing[2])]
    boxgridsize = [int(box_xdim / gridspacing[0]), int(box_ydim / gridspacing[1]), int(box_zdim / gridspacing[2])]
    # 2. copy section of input map into new object
    mapfragment = mapgrid.get_subarray(start=startpoint, shape=boxgridsize)
    if DEBUG:
        print("Extracted map fragment props:xyz of starting point and box size")
        print(list([*startpoint,*boxgridsize]))
    ccp4map = gemmi.Ccp4Map()
    ccp4map.grid = gemmi.FloatGrid(mapfragment)
    # 3. adjust headers by specifying
    ccp4map.update_ccp4_header()
    ccp4map.grid.spacegroup=gemmi.SpaceGroup('P1')

    '''
    set size (1-3) and position(5-7) of map fragment inside ASU
    1      NC              # of Columns    (fastest changing in map)
    2      NR              # of Rows
    3      NS              # of Sections   (slowest changing in map)
    4      2 for REAL values of density map
    5      NCSTART         Number of first COLUMN  in map
    6      NRSTART         Number of first ROW     in map
    7      NSSTART         Number of first SECTION in map
    '''
    for i,value in enumerate([*boxgridsize, 2, *startpoint]):
        # x,y,z of start point
        #
        ccp4map.set_header_i32(i+1, value)

    '''
    copy original unit cell params
    8      NX              Number of intervals along X
    9      NY              Number of intervals along Y
    10      NZ              Number of intervals along Z

    FROM: https://www.ccp4.ac.uk/html/maplib.html
    '''
    #
    for i in [8,9,10]:
        ccp4map.set_header_i32(i,map_obj.header_i32(i))

    # copy unit cell params
    '''
    11      X length        Cell Dimensions (Angstroms)
    12      Y length                     "
    13      Z length                     "
    14      Alpha           Cell Angles     (Degrees)
    15      Beta                         "
    16      Gamma                        "
    '''
    for i in range(11,17):
        ccp4map.set_header_float(i,map_obj.header_float(i))


    ccp4map.update_ccp4_header()
    ccp4map.write_ccp4_map(mapout)
    if savedmap:
        ccp4map.write_ccp4_map(savedmap)

def load_mtz_map_fragment(mtzfilename:str, mapobjname:str, margin=3, savedmap='') -> None:
    '''
    DESCRIPTION
    :mtzfilename: path to mtz map file
    :mapobjname: arbitrary name of map entry
    :margin: map cutoff distance around the ligand in angstroms
    The subroutine extracts map fragment around ligand atoms using the following protocol
    Algorithm:
    1. Generate 2Fo-Fc map from the specified MTZ file :mtzfilename:
    2. Save  fragment of the map around the ligand into temporary file
    3. load map into session

    NOTE: this function requires rw access to the folder with ddtrek file
    '''
    ligand = gemmi.read_structure('ligand.pdb') # if map file was specified, this file is generated in main body of DDtrek 

    if mtzfilename.endswith('mtz'):
        mtz = gemmi.read_mtz_file(mtzfilename)
        m = gemmi.Ccp4Map()
        # ATM only 2Fo-Fc maps with column names 2FOFCWT or FTW(Refmac5)  are recognized
        # TODO: add option to select type of map to draw
        # Option to add:
        # omit map "mFo-DFc_omit"  and phi "PHImFo-DFc_omit"
        if 'polder' in mtzfilename:
            #if filename contains word Polder, then extract omit map
            m.grid = mtz.transform_f_phi_to_map('mFo-DFc_omit', 'PHImFo-DFc_omit', sample_rate=3)
        else:
            try:
                m.grid = mtz.transform_f_phi_to_map('2FOFCWT','PH2FOFCWT', sample_rate=3)# column labels for mtz with map coefficients in Phenix and Buster
            except:
                m.grid = mtz.transform_f_phi_to_map('FWT','PHWT', sample_rate=3)#refmac5 default names for map coefficients
        m.update_ccp4_header()
        # command below seems to generate Gb files of map fragments and overflows memory in PyMOL < 2.5
        m.set_extent(ligand.calculate_fractional_box(margin=margin)) #cut map fragment with margin around the ligand
        m.write_ccp4_map('masked.ccp4')
        if savedmap:
            m.write_ccp4_map(savedmap)

    elif mtzfilename.endswith(('ccp4','map','map.gz')):
        # load cryoEM map
        load_cryoem_map_fragment(mtzfilename,savedmap)
    cmd.load('masked.ccp4', mapobjname)
    # remove temporary files
    os.remove('masked.ccp4')
    os.remove('ligand.pdb')
    return

def ddtrek(input_filename: str, coordinate_cutoff = 7, map_cutoff = 7) -> None:
    '''
    DESCRIPTION

    Aligns structures from input file to the reference structure. Optionally generates electron density maps
    
    See publication: DDtrek (2024)

    USAGE

    ddtrek [ path to input file ]

    Variables
    :input_filename: path (absolute or relative) to DDtrek.in file
    :coordinate_cutoff: distance from ligand in angstrom defines truncation radius
    :map_cutoff: cutoff distance for map fragment. Passed to load_mtz_map_fragment as :margin:
    :DEBUG: produces extra output into PyMOL terminal
    '''

    #Change working directory to location of input file. Folder should be writable
    abspath = os.path.abspath(input_filename)
    dname = os.path.dirname(abspath)
    fname = os.path.basename(input_filename)
    os.chdir(dname)
    assert (os.access(dname, os.W_OK)), "Folder with DDtrek input file should be writable"
    ## get the list of preloaded structures in currently opened PyMOL window
    preloaded_structures = cmd.get_object_list()
    # set default group name to avoid errors if user forgot to specify the name by #G
    group_name = 'default'

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
            print('Malformed line with PDB %s . Check pdb filename and presence of all elements. Skipping...' % pdb )
            continue

        # load map mtz and check extension _ only mtz maps are allowed
        # TODO: may be rewrite this section. It looks terrible BUT it works...
        mtz = entry.split()[1]
        if not mtz.lower().endswith(('mtz','map','ccp4')):
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
        cmd.symexp(prefix='symmetry',object='current_entry',selection=ligand_residue, cutoff=coordinate_cutoff, segi=1)
        # remove  atoms in symmetry mates beyond 6 A from the ligand in order to avoid unnecessary calculations
        cmd.remove("%s beyond 6 of %s" % ('symmetry*',ligand_residue))

        
        # Now select the chain that contains the ligand, ligand itself and everything in coordinate_cutoff radius around the ligand of interest. 
        # Select ligand, all closely situated protein chains and water molecules
        # TODO reformat string into new python way f'' and add cutoff distances as variables here
        tmp_selection = "current_entry and (%s or ((bychain %s around 7) and polymer) or (%s around 7 and resn HOH))" % (ligand_residue,ligand_residue,ligand_residue)
        # now add symmetry mates to it
        extracted_selection = "(%s) or symmetry*" % tmp_selection


        if DEBUG:
            print('selection:%s' % extracted_selection)

        #Extract coordinates into the new entry
        cmd.create(entry_name, selection=extracted_selection,
                    source_state=0, target_state=0,discrete=0)

        # Prior to alignment, check that reference is loaded. If not, create reference from 1st polymer chain of the currently loaded entry
        if 'reference' not in cmd.get_object_list():
            reference_chain = cmd.get_chains('current_entry and polymer')[0]
            cmd.create('reference',f'current_entry and polymer and chain {reference_chain}')
        #ALIGN newly created structure to reference using all CA atoms
        if align_chain:
            alignment_selection = "%s and polymer and chain %s and name CA" % (entry_name, align_chain)
        else:
            alignment_selection = "%s and polymer and name CA" % (entry_name)
        if DEBUG:
            print('Alignment selection:%s\n' % alignment_selection)
        cmd.super(alignment_selection, 'reference and name CA') # align using C-alpha atoms to reference

        #Cleanup. Remove temporary PDB and its generated symmetry mates
        cmd.delete('current_entry') # remove original PDB
        cmd.delete('symmetry*') # remove symmetry mates
        
        #Truncate atoms in new entry. Limit by 7 A radius around the ligand
        pocket_selection = "(byres(%s within 7 of (%s and chain %s and resi %s)))" % (entry_name,entry_name,chain, resi)
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
        if 'omit' in map_name:
            #contour omit maps at 3 sigma
            cmd.isomesh(map_name, mtz_name, 3, ligand_selection, carve=1.8)
        else:
            #2Fo-Fc maps are contoured at 1 sigma level
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
cmd.extend("map_extract", map_extract)