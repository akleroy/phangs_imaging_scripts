# Pure python parts of the PHANGS pipeline, used in some of our calls
# that run CASA scripts, do outside analysis, etc.

import os
import numpy as np
import scipy.ndimage as ndimage
import glob

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# Interface to the text file keys that steer the process
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

def read_ms_key(fname='../scripts/ms_file_key.txt'):
    """
    Read the measurement set key into a big dictionary. This maps
    locations of reduced files to galaxy, project, and data set name.
    """
    infile = open(fname, 'r')

    ms_key = {}

    while True:
        line  = infile.readline()    
        if len(line) == 0:
            break
        if line[0] == '#':
            continue
        words = line.split()
        if len(words) < 4:
            continue

        this_gal = words[0]
        this_proj = words[1]
        this_ms = words[2]
        this_file = words[3]

        if (this_gal in ms_key.keys()) == False:
            ms_key[this_gal] = {}
        if ms_key[this_gal].has_key(this_proj) == False:
            ms_key[this_gal][this_proj] = {}
        ms_key[this_gal][this_proj][this_ms] = this_file
        
    infile.close()
    
    return ms_key

def get_uvdata_key(gal=None,
                   just_proj=None,
                   just_ms=None,
                   just_array=None,
                   quiet=False):
    """
    Figure out the root directory for calibrating the uv data for a
    galaxy. Reads from the ms_file_key and works backwards from the
    calibrated data in that directory.
    """

    if gal == None:
        if quiet == False:
            print("Please specify a galaxy.")
        return None

    ms_key = read_ms_key()
    
    # Get the measurement set for this specific galaxy

    if ms_key.has_key(gal) == False:
        if quiet == False:
            print("Galaxy "+gal+" not found in the measurement set key.")
        return None
    gal_specific_key = ms_key[gal]

    uvdata_dict = {}

    # Loop over projects in the measurement set key

    for this_proj in gal_specific_key.keys():

        # If a project is specificied skip all but the relevant project

        if just_proj != None:
            if type(just_proj) == type([]):
                if just_proj.count(this_proj) == 0:
                    continue
            else:
                if this_proj != just_proj:
                    continue

        proj_specific_key = gal_specific_key[this_proj]

        # Loop over all measurement sets in the project key

        for this_ms in proj_specific_key.keys():

            if just_ms != None:
                if type(just_ms) == type([]):
                    if just_ms.count(this_ms) == 0:
                        continue
                    else:
                        if this_ms != just_ms:
                            continue
 
            if just_array != None:
                if this_ms.count(just_array) == 0:
                    continue
                
            if this_ms.count('7m') == 1:
                this_array = '7m'
            elif this_ms.count('12m') == 1:
                this_array = '12m'
            else:
                this_array = '???'

            # We now have the calibrated file name

            this_calibrated_file = proj_specific_key[this_ms]

            components = this_calibrated_file.split('/calibrated/')
            if len(components) != 2:
                print("")
                print("WARNING! Something is wrong with file "+this_calibrated_file)
                print("We assume that there is one and only one /calibrated/ in the directory.")
                print("Fix this or whatever else is going wrong and rerun. Skipping for now.")
                print("")
                continue

            this_dir = components[0]+'/'
            this_uid = components[1]

            if uvdata_dict.has_key(this_dir):
                current_uids = uvdata_dict[this_dir]['uid']
                if current_uids.count(this_uid) == 0:
                    current_uids.append(this_uid)
                    uvdata_dict[this_dir]['uid'] = current_uids
            else:
                uvdata_dict[this_dir] = {'gal':gal,
                                         'array':this_array,
                                         'uid':[this_uid]}

    return uvdata_dict

def read_dir_key(fname='../scripts/dir_key.txt'):
    """
    Read the directory key, which gives us a general way to sort out
    which MS files go in which directory. This is relevant mainly for
    cases where multiple science goals target a single galaxy. In
    those cases, we want the whole galaxy in one directory. This gives
    a way to map, e.g., ngc3627north to ngc3627/
    """
    infile = open(fname, 'r')
    
    dir_key = {}

    while True:
        line  = infile.readline()    
        if len(line) == 0:
            break
        if line[0] == '#':
            continue
        words = line.split()
        if len(words) < 2:
            continue

        this_ms = words[0]
        this_dir = words[1]    
        dir_key[this_ms] = this_dir
        
    infile.close()
    
    return dir_key

def dir_for_gal(gal=None,
                quiet=False):
    """
    Return the working directory given a galaxy name. See above.
    """

    if gal == None:
        if quiet == False:
            print("Please specify a galaxy.")
        return

    dir_key = read_dir_key()
    if gal in dir_key.keys():
        this_dir = '../'+dir_key[gal]+'/'
    else:
        this_dir = '../'+gal+'/'

    return this_dir

def list_gal_names():
    """
    List the full set of galaxy names known from the ms_file_key.
    """
    ms_key = read_ms_key()
    gal_names = ms_key.keys()
    gal_names.sort()
    return gal_names

def list_wholegal_names():
    """
    List only full galaxy names (no parts, e.g., _1 or _2).
    """
    ms_key = read_ms_key()
    dir_key = read_dir_key()
    part_names = ms_key.keys()
    wholegal_names = []
    for this_part in part_names:
        if this_part in dir_key.keys():
            this_name = dir_key[this_part]
        else:
            this_name = this_part
        if wholegal_names.count(this_name) > 0:
            continue
        wholegal_names.append(this_name)
    wholegal_names.sort()
    return wholegal_names    

def mosaic_key():
    """
    List the parts that contribute to mosaic.
    """
    ms_key = read_ms_key()
    dir_key = read_dir_key()
    part_names = ms_key.keys()

    mosaic_key = {}
    for this_part in part_names:
        if this_part in dir_key.keys():
            this_name = dir_key[this_part]
            if this_name in mosaic_key.keys():
                current_list = mosaic_key[this_name]
                if this_part in current_list:
                    continue
                else:
                    current_list.append(this_part)
                mosaic_key[this_name] = current_list
            else:
                mosaic_key[this_name] = [this_part]
        else:
            continue

    return mosaic_key

def read_mosaic_key(fname='../scripts/mosaic_definitions.txt'):
    """
    Read the file containing the centers and velocities for each
    mosaic. Note that for cases where the galaxy is observed several
    times, the RA and DEC refer to the intended center of the mosaic,
    NOT the center of the galaxy. In other cases, we tend use the NED
    center.
    """
    infile = open(fname, 'r')

    mosaic_key = {}

    while True:
        line  = infile.readline()    
        if len(line) == 0:
            break
        if line[0] == '#':
            continue
        words = line.split()

        if len(words) < 5:
            continue

        this_gal = words[0]
        this_ra = words[1]
        this_dec = words[2]
        this_vsys = words[3]
        this_vwidth = words[4]
        
        mosaic_key[this_gal] = {}
        mosaic_key[this_gal]['rastring'] = this_ra
        mosaic_key[this_gal]['decstring'] = this_dec
        mosaic_key[this_gal]['vsys'] = float(this_vsys)
        mosaic_key[this_gal]['vwidth'] = float(this_vwidth)

    infile.close()
    
    return mosaic_key

def read_multipart_key(fname='../scripts/multipart_fields.txt'):
    """
    Read the file containing the definitions of the mosaics
    for the multi-part fields.
    """
    infile = open(fname, 'r')
    
    multipart_key = {}

    while True:
        line  = infile.readline()    
        if len(line) == 0:
            break
        if line[0] == '#':
            continue
        words = line.split()

        if len(words) < 5:
            continue

        this_gal = words[0]
        this_ra_ctr = words[1]
        this_dec_ctr = words[2]
        this_deltara = words[3]
        this_deltadec = words[4]
        
        multipart_key[this_gal] = {}
        multipart_key[this_gal]['ra_ctr_deg'] = float(this_ra_ctr)
        multipart_key[this_gal]['dec_ctr_deg'] = float(this_dec_ctr)
        multipart_key[this_gal]['delta_ra_as'] = float(this_deltara)
        multipart_key[this_gal]['delta_dec_as'] = float(this_deltadec)

    infile.close()
    
    return multipart_key

def read_override_mosaic_params(
    fname='../scripts/override_mosaic_params.txt'
    ):
    """
    Read hand set overrides for cell and image sizes.
    """

    infile = open(fname, 'r')

    override_dict = {}
    while True:
        line = infile.readline()    
        if len(line) == 0:
            break
        if line[0] == '#':
            continue
        words = line.split()
        if len(words) < 3:
            continue
        vis_override = words[0]
        param_override = words[1]
        value_override = words[2]
        if override_dict.has_key(vis_override) == False:
            override_dict[vis_override] = {}
        override_dict[vis_override][param_override] = value_override
            
    infile.close()

    return override_dict

def read_override_imaging_params(
    fname='../scripts/override_imaging_params.txt'
    ):
    """
    Read hand set overrides for imaging parameters.
    """

    infile = open(fname, 'r')

    override_dict = {}
    while True:
        line = infile.readline()    
        if len(line) == 0:
            break
        if line[0] == '#':
            continue
        words = line.split()
        if len(words) < 3:
            continue
        vis_override = words[0]
        param_override = words[1]
        value_override = words[2]
        if override_dict.has_key(vis_override) == False:
            override_dict[vis_override] = {}
        override_dict[vis_override][param_override] = value_override
            
    infile.close()

    return override_dict

def read_singledish_key(
    fname='../scripts/singledish_key.txt'
    ):
    """
    Read the single dish key.
    """
    infile = open(fname, 'r')

    singledish_key = {}
    while True:
        line = infile.readline()    
        if len(line) == 0:
            break
        if line[0] == '#':
            continue
        words = line.split()
        if len(words) < 3:
            continue

        this_gal = words[0]
        this_file = words[1]
        this_product = words[2]
        
        if (this_gal in singledish_key.keys()) == False:
            singledish_key[this_gal] = {}
        singledish_key[this_gal][this_product] = this_file
            
    infile.close()

    return singledish_key
