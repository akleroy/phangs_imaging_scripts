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

        if ms_key.has_key(this_gal) == False:
            ms_key[this_gal] = {}
        if ms_key[this_gal].has_key(this_proj) == False:
            ms_key[this_gal][this_proj] = {}
        ms_key[this_gal][this_proj][this_ms] = this_file
        
    infile.close()
    
    return ms_key

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

def dir_for_gal(gal=None):
    """
    Return the working directory given a galaxy name. See above.
    """

    if gal == None:
        if quiet == False:
            print "Please specify a galaxy."
        return

    dir_key = read_dir_key()
    if dir_key.has_key(gal):
        this_dir = '../'+dir_key[gal]+'/'
    else:
        this_dir = '../'+gal+'/'

    return this_dir

def list_gal_names():
    """
    List the full set of galaxy names known from the ms_file_key
    """
    ms_key = read_ms_key()
    gal_names = ms_key.keys()
    gal_names.sort()
    return gal_names

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
