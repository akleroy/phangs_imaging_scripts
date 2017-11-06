# NB: CASA doesn't always include the pwd in the module search path. I
# had to modify my init.py file to get this to import.

import os
import numpy as np
import glob

# Other PHANGS scripts
import line_list

# CASA imports
import analysisUtils as au
from split import split
from statwt import statwt
from mstransform import mstransform
from concat import concat
from uvcontsub import uvcontsub
from flagdata import flagdata
from imstat import imstat
from tclean import tclean
from exportfits import exportfits

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

def read_override_image_params(
    fname='../scripts/override_image_params.txt'):
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

def read_imaging_key(fname='../scripts/imaging_key.txt'):
    """
    Read the key used to guide PHANGS imaging.
    """
    infile = open(fname, 'r')
    
    imaging_key = {}

    while True:
        line  = infile.readline()    
        if len(line) == 0:
            break
        if line[0] == '#':
            continue
        words = line.split()
        if len(words) < 7:
            continue

        this_input_vis = words[0]
        this_cube_root = words[1]
        this_pb_limit = words[2]
        this_multiscale_snr_thresh = words[3]
        this_smallscalebias = words[4]
        this_clean_mask = words[5]
        this_singlescale_snr_thresh = words[6]

        imaging_key[this_input_vis] = {}
        imaging_key[this_input_vis]['cube_root'] = this_cube_root
        imaging_key[this_input_vis]['pb_limit'] = this_pb_limit
        imaging_key[this_input_vis]['multiscale_thresh'] = this_multiscale_snr_thresh
        imaging_key[this_input_vis]['smallscalebias'] = this_smallscalebias
        imaging_key[this_input_vis]['clean_mask'] = this_clean_mask
        imaging_key[this_input_vis]['singlescale_thresh'] = this_singlescale_snr_thresh

    infile.close()

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# Routines to move data around.
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

# All of these know about the PHANGS keys. They're called as part of
# the pipeline to set up the imaging.

def copy_data(gal=None,
              just_proj=None,
              just_ms=None,
              just_array=None,
              do_split=True,
              do_statwt=True,
              quiet=False):
    """
    Copies data from its original location, which is specified in a
    text file ms_key.txt. Then splits out only the science target.
    """

    if gal == None:
        if quiet == False:
            print "Please specify a galaxy."
        return

    ms_key = read_ms_key()

    if ms_key.has_key(gal) == False:
        if quiet == False:
            print "Galaxy "+gal+" not found in the measurement set key."
        return
    gal_specific_key = ms_key[gal]

    # Change to the right directory

    this_dir = dir_for_gal(gal)
    os.chdir(this_dir)

    if quiet == False:
        print "--------------------------------------------------------"
        print "START: Copying data from original location."
        print "--------------------------------------------------------"

        print "Galaxy: ", gal
        if just_array != None: print "Project: ", just_proj
        if just_array != None: print "Measurements Set: ", just_ms
        if just_array != None: print "Array: ", just_array

    # Loop over files in the measurement set key

    for this_proj in gal_specific_key.keys():
        if just_proj != None:
            if type(just_proj) == type([]):
                if just_proj.count(this_proj) == 0:
                    continue
            else:
                if this_proj != just_proj:
                    continue

        proj_specific_key = gal_specific_key[this_proj]
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
           
            # Set up a copy command, overwriting previous versions

            in_file = proj_specific_key[this_ms]            

            if do_split:
                copied_file = gal+'_'+this_proj+'_'+this_ms+'_copied.ms'
            else:
                copied_file = gal+'_'+this_proj+'_'+this_ms+'.ms'

            os.system('rm -rf '+copied_file)
            os.system('rm -rf '+copied_file+'.flagversions')

            command = 'cp -r -H '+in_file+' '+copied_file
            print command
            var = os.system(command)    
            print var

            # Call split and statwt if desired (default is yes)

            if do_split:

                if quiet == False:
                    print "Splitting out science target data."

                out_file = gal+'_'+this_proj+'_'+this_ms+'.ms'

                os.system('rm -rf '+out_file)
                os.system('rm -rf '+out_file+'.flagversions')

                split(vis=copied_file
                      , intent ='OBSERVE_TARGET#ON_SOURCE'
                      , datacolumn='DATA'
                      , outputvis=out_file)        

                os.system('rm -rf '+copied_file)
                os.system('rm -rf '+copied_file+'.flagversions')

            if do_statwt:

                if quiet == False:
                    print "Using statwt to re-weight the data."

                statwt(vis=out_file,
                       datacolumn='DATA')

    if quiet ==False:
        print "--------------------------------------------------------"
        print "END: Copying data from original location."
        print "--------------------------------------------------------"

def concat_line_for_gal(
    gal=None,
    just_proj=None,
    just_ms=None,
    just_array=None,
    line='co21',
    tag='',
    do_statwt=True,
    do_chan0=True,
    quiet=False):
    """
    Combine all measurement sets for one line and one galaxy.
    """

    # Identify the data sets to combine

    if gal == None:
        if quiet == False:
            print "Please specify a galaxy."
        return

    ms_key = read_ms_key()

    if ms_key.has_key(gal) == False:
        if quiet == False:
            print "Galaxy "+gal+" not found in the measurement set key."
        return
    gal_specific_key = ms_key[gal]

    files_to_concat = []

    for this_proj in gal_specific_key.keys():
        if just_proj != None:
            if type(just_proj) == type([]):
                if just_proj.count(this_proj) == 0:
                    continue
            else:
                if this_proj != just_proj:
                    continue

        proj_specific_key = gal_specific_key[this_proj]
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

            this_in_file = gal+'_'+this_proj+'_'+this_ms+'_'+line+'.ms'    
            if os.path.isdir(this_in_file) == False:
                continue
            files_to_concat.append(this_in_file)

    if len(files_to_concat) == 0:
        print "No files to concatenate found. Returning."
        return

    # Concatenate all of the relevant files

    if tag != '':
        out_file =  gal+'_'+tag+'_'+line+'.ms'
    else:
        out_file =  gal+'_'+line+'.ms'

    os.system('rm -rf '+out_file)
    os.system('rm -rf '+out_file+'.flagversions')

    concat(vis=files_to_concat,
           concatvis=out_file)

    # Re-weight the data empirically

    if do_statwt:
        statwt(vis=out_file,
               datacolumn='DATA')

    # Collapse to form a "channel 0" measurement set

    if do_chan0 == False:
        return
    
    if tag != '':
        chan0_vis = gal+'_'+tag+'_'+line+'_chan0.ms'
    else:
        chan0_vis = gal+'_'+line+'_chan0.ms'

    os.system('rm -rf '+chan0_vis)
    os.system('rm -rf '+chan0_vis+'.flagversions')
    split(vis=out_file
          , datacolumn='DATA'
          , spw=''
          , outputvis=chan0_vis
          , width=10000)

def concat_cont_for_gal(
    gal=None,
    just_proj=None,
    just_ms=None,
    just_array=None,
    tag='',
    ):
    """
    Concatenate continuum data sets.
    """
    pass

    # Identify the data sets to combine

    if gal == None:
        if quiet == False:
            print "Please specify a galaxy."
        return

    ms_key = read_ms_key()

    if ms_key.has_key(gal) == False:
        if quiet == False:
            print "Galaxy "+gal+" not found in the measurement set key."
        return
    gal_specific_key = ms_key[gal]

    files_to_concat = []

    for this_proj in gal_specific_key.keys():
        if just_proj != None:
            if type(just_proj) == type([]):
                if just_proj.count(this_proj) == 0:
                    continue
            else:
                if this_proj != just_proj:
                    continue

        proj_specific_key = gal_specific_key[this_proj]
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

            this_in_file = gal+'_'+this_proj+'_'+this_ms+'_cont.ms'
            if os.path.isdir(this_in_file) == False:
                continue
            files_to_concat.append(this_in_file)

    if len(files_to_concat) == 0:
        print "No files to concatenate found. Returning."
        return

    # Concatenate all of the relevant files

    if tag != '':
        out_file =  gal+'_'+tag+'_cont.ms'
    else:
        out_file =  gal+'_cont.ms'

    os.system('rm -rf '+out_file)
    os.system('rm -rf '+out_file+'.flagversions')

    concat(vis=files_to_concat,
           concatvis=out_file)

def extract_phangs_lines(   
    gal=None,
    just_array=None,
    ext='',
    quiet=False
    ):
    """
    Extract all phangs lines and continuum for a galaxy.
    """

    # Could add sio54, which is generally covered in PHANGS but almost
    # always likely to be a nondetection.

    if quiet == False:
        print "--------------------------------------------------------"
        print "START: Extracting spectral lines from data set."
        print "--------------------------------------------------------"

    chan_width = {}
    chan_width['co21'] = 2.5
    chan_width['c18o21'] = 6.0

    for line in ['co21', 'c18o21']:

        extract_line_for_galaxy(   
            gal=gal,
            just_array=just_array,
            line=line,
            ext=ext,
            chan_width=chan_width[line],    
            quiet=quiet
            )

        if just_array != '12m':
            concat_line_for_gal(
                gal=gal,
                just_array='7m',
                tag='7m',
                line=line,
                do_statwt=True,
                do_chan0=True)

        if just_array != '7m':
            concat_line_for_gal(
                gal=gal,
                just_array='12m',
                tag='12m',
                line=line,
                do_statwt=True,
                do_chan0=True)

        has_7m = len(glob.glob(gal+'*7m*'+line+'*')) > 0
        has_12m = len(glob.glob(gal+'*12m*'+line+'*')) > 0
        if has_12m == False or has_7m == False:
            continue

        if just_array == None:

            concat_line_for_gal(
                gal=gal,
                just_array = None,
                tag='12m+7m',
                line=line,
                do_statwt=True,
                do_chan0=True)
            

    if quiet == False:
        print "--------------------------------------------------------"
        print "END: Extracting spectral lines from data set."
        print "--------------------------------------------------------"

def extract_phangs_continuum(   
    gal=None,
    just_array=None,
    ext='',
    quiet=False
    ):
    """
    Extract all phangs lines and continuum for a galaxy.
    """

    if quiet == False:
        print "--------------------------------------------------------"
        print "START: Extracting continuum from data set."
        print "--------------------------------------------------------"

    lines_to_flag = line_list.lines_co+line_list.lines_13co+line_list.lines_c18o

    extract_continuum_for_galaxy(   
        gal=gal,
        just_array=just_array,
        lines_to_flag=lines_to_flag,
        ext=ext,
        do_statwt=True,
        do_collapse=True,
        quiet=quiet
        )

    if just_array != '12m':
        concat_cont_for_gal(
            gal=gal,
            just_array = '7m',
            tag = '7m')

    if just_array != '7m':
        concat_cont_for_gal(
            gal=gal,
            just_array = '12m',
            tag = '12m')

    has_7m = len(glob.glob(gal+'*7m*cont*')) > 0
    has_12m = len(glob.glob(gal+'*12m*cont*')) > 0

    if just_array == None and has_7m and has_12m:
        concat_cont_for_gal(
            gal=gal,
            just_array = None,
            tag = '12m+7m')

    if quiet == False:
        print "--------------------------------------------------------"
        print "END: Extracting continuum from data set."
        print "--------------------------------------------------------"

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# Routines to extract lines from a measurement set
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

def list_lines_in_ms(
    in_file= None,
    vsys=0.0,
    gal=None,
    ):    
    """
    List the lines likely to be present in a measurement set. This can
    be a general purpose utility.
    """

    sol_kms = 2.99e5

    # pull the parameters from the galaxy in the mosaic file
    if gal != None:
        mosaic_parms = read_mosaic_key()
        if mosaic_parms.has_key(gal):
            vsys = mosaic_parms[gal]['vsys']
            vwidth = mosaic_parms[gal]['vwidth']

    # Set up the input file

    if os.path.isdir(in_file) == False:
        if quiet == False:
            print "Input file not found."
        return

    lines_in_ms = []
    for line in line_list.line_list.keys():
        restfreq_ghz = line_list.line_list[line]

        # work out the frequency of the line and the line wings

        target_freq_ghz = restfreq_ghz*(1.-vsys/sol_kms)

        this_spw_list = au.getScienceSpwsForFrequency(in_file, target_freq_ghz*1e9)    
        if len(this_spw_list) == 0:
            continue
        
        lines_in_ms.append(line)

    return lines_in_ms

def chanwidth_for_line(
    in_file=None,
    line='co21',
    gal=None,
    vsys=0.0,
    vwidth=500.,
    quiet=False):
    """
    Return the coarsest channel width among spectral windows that
    overlap a line. This can be a general purpose utility.
    """

    # pull the parameters from the galaxy in the mosaic file
    if gal != None:
        mosaic_parms = read_mosaic_key()
        if mosaic_parms.has_key(gal):
            vsys = mosaic_parms[gal]['vsys']
            vwidth = mosaic_parms[gal]['vwidth']

    sol_kms = 2.99e5

    # Set up the input file

    if os.path.isdir(in_file) == False:
        if quiet == False:
            print "Input file not found."
        return

    # Look up the line

    if line_list.line_list.has_key(line) == False:
        if quiet == False:
            print "Line not found. Give lower case abbreviate found in line_list.py"
        return
    restfreq_ghz = line_list.line_list[line]

    # Work out which spectral windows contain the line contain

    target_freq_ghz = restfreq_ghz*(1.-vsys/sol_kms)
    target_freq_high = restfreq_ghz*(1.-(vsys-0.5*vwidth)/sol_kms)
    target_freq_low = restfreq_ghz*(1.-(vsys+0.5*vwidth)/sol_kms)

    spw_list_string = ''    
    first = True
    spw_list = []

    for target_freq in [target_freq_high, target_freq_ghz, target_freq_low]:
        this_spw_list = au.getScienceSpwsForFrequency(in_file, target_freq*1e9)    
        for spw in this_spw_list:
            if spw_list.count(spw) != 0:
                continue
            spw_list.append(spw)
            if not first:
                spw_list_string += ','
            else:
                first = False
            spw_list_string += str(spw)

    if len(spw_list) == 0:
        if quiet == False:
            print "No spectral windows contain this line at this redshift."
        return

    # Figure out how much averaging is needed to reach the target resolution
    chan_width_hz = au.getChanWidths(in_file, spw_list_string)

    # Convert to km/s and return
    chan_width_kms = abs(chan_width_hz / (restfreq_ghz*1e9)*sol_kms)

    return chan_width_kms

def extract_line(in_file=None,
                 out_file=None,
                 line='co21',
                 gal=None,
                 vsys=0.0,
                 vwidth=500.,
                 chan_width=2.5,
                 quiet=False):
    """
    Extract a spectral line from a measurement set and regrid onto a
    new velocity grid with the desired spacing. This doesn't
    necessarily need the PHANGS keys in place and may be a general
    purpose utility.
    """

    sol_kms = 2.99e5

    # pull the parameters from the galaxy in the mosaic file
    if gal != None:
        mosaic_parms = read_mosaic_key()
        if mosaic_parms.has_key(gal):
            vsys = mosaic_parms[gal]['vsys']
            vwidth = mosaic_parms[gal]['vwidth']

    # Set up the input file

    if os.path.isdir(in_file) == False:
        if quiet == False:
            print "Input file not found."
        return

    # Look up the line

    if line_list.line_list.has_key(line) == False:
        if quiet == False:
            print "Line not found. Give lower case abbreviate found in line_list.py"
        return
    restfreq_ghz = line_list.line_list[line]

    # Work out which spectral windows contain the line contain

    target_freq_ghz = restfreq_ghz*(1.-vsys/sol_kms)
    target_freq_high = restfreq_ghz*(1.-(vsys-0.5*vwidth)/sol_kms)
    target_freq_low = restfreq_ghz*(1.-(vsys+0.5*vwidth)/sol_kms)

    spw_list_string = ''    
    first = True
    spw_list = []

    for target_freq in [target_freq_high, target_freq_ghz, target_freq_low]:
        this_spw_list = au.getScienceSpwsForFrequency(in_file, target_freq*1e9)    
        for spw in this_spw_list:
            if spw_list.count(spw) != 0:
                continue
            spw_list.append(spw)
            if not first:
                spw_list_string += ','
            else:
                first = False
            spw_list_string += str(spw)

    if len(spw_list) == 0:
        if quiet == False:
            print "No spectral windows contain this line at this redshift."
        return

    # Figure out the starting velocity and number of channels.

    start_vel_kms = (vsys - vwidth/2.0)
    nchan = int(np.floor(vwidth / chan_width)+1)

    # Convert to strings. The precision is hardcoded with
    # extragalactic ALMA work in mind.

    restfreq_string = "{:10.6f}".format(restfreq_ghz)+'GHz'
    chan_dv_string =  "{:5.2f}".format(chan_width)+'km/s'
    start_vel_string =  "{:6.1f}".format(start_vel_kms)+'km/s'

    # Figure out how much averaging is needed to reach the target resolution

    chan_width_hz = au.getChanWidths(in_file, spw_list_string)
    current_chan_width_kms = abs(chan_width_hz / (restfreq_ghz*1e9)*sol_kms)    
    target_width_hz = chan_width/sol_kms*restfreq_ghz*1e9
    rebin_factor = min(target_width_hz / chan_width_hz)

    if current_chan_width_kms > chan_width:
        print "Requested channel width is smaller than the starting width. Returning."
        return
    
    if rebin_factor < 2:
        chanbin = 1
    else:
        chanbin = int(np.floor(rebin_factor/2.))

    if chanbin > 1:
        chanaverage = True
    else:
        chanaverage = False

    # Report the call

    if quiet == False:        
        print "FILE:", in_file
        print "LINE TAG: ", line
        print "REST FREQUENCY: ", restfreq_string
        print "SPECTRAL WINDOWS: ", spw_list
        print "STARTING CHANNEL WIDTH: ", current_chan_width_kms
        print "TARGET CHANNEL WIDTH: ", chan_dv_string
        print "SUPPLIED SOURCE VELOCITY: ", str(vsys)
        print "DESIRED VELOCITY WIDTH: ", str(vwidth)
        print "START VELOCITY: ", start_vel_string
        print "NUMBER OF CHANNELS: ", str(nchan)
        print "CHANNELS TO BIN TOGETHER FIRST: ", chanbin

    # Call mstransform
    
    os.system('rm -rf '+out_file+'.temp')
    os.system('rm -rf '+out_file+'.temp.flagversions')

    mstransform(vis=in_file,
                outputvis=out_file+'.temp',
                spw=spw_list_string,
                datacolumn='DATA',
                chanaverage=chanaverage,
                chanbin=chanbin,
                hanning=False,
                interpolation='cubic',
                )

    os.system('rm -rf '+out_file)
    os.system('rm -rf '+out_file+'.flagversions')

    mstransform(vis=out_file+'.temp',
                outputvis=out_file,
                datacolumn='DATA',
                combinespws=True,
                regridms=True,
                mode='velocity',
                interpolation='cubic',
                start=start_vel_string,
                nchan=nchan,
                width=chan_dv_string,
                restfreq=restfreq_string,
                outframe='lsrk',
                veltype='radio',
                )

    os.system('rm -rf '+out_file+'.temp')
    os.system('rm -rf '+out_file+'.temp.flagversions')

    return

def extract_line_for_galaxy(   
    gal=None,
    just_proj=None,
    just_ms=None,
    just_array=None,
    line='co21',
    vsys=0.0,
    vwidth=500.,
    chan_width=2.5,    
    ext='',
    quiet=False
    ):
    """
    Extract a given line for all data sets for a galaxy. This knows
    about the PHANGS measurement set keys and is specific to our
    projects.
    """
    
    if gal == None:
        if quiet == False:
            print "Please specify a galaxy."
        return

    ms_key = read_ms_key()

    if ms_key.has_key(gal) == False:
        if quiet == False:
            print "Galaxy "+gal+" not found in the measurement set key."
        return
    gal_specific_key = ms_key[gal]

    # Look up the galaxy specific parameters

    mosaic_parms = read_mosaic_key()
    if mosaic_parms.has_key(gal):
        vsys = mosaic_parms[gal]['vsys']
        vwidth = mosaic_parms[gal]['vwidth']

    # Change to the right directory

    this_dir = dir_for_gal(gal)
    os.chdir(this_dir)

    # Loop over all projects and measurement sets

    for this_proj in gal_specific_key.keys():

        if just_proj != None:
            if type(just_proj) == type([]):
                if just_proj.count(this_proj) == 0:
                    continue
            else:
                if this_proj != just_proj:
                    continue

        proj_specific_key = gal_specific_key[this_proj]
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
            
            in_file = gal+'_'+this_proj+'_'+this_ms+ext+'.ms'
            out_file = gal+'_'+this_proj+'_'+this_ms+'_'+line+'.ms'    

            lines_in_ms = list_lines_in_ms(in_file, gal=gal)
            if lines_in_ms.count(line) == 0:
                print "Line not found in measurement set."
                return

            extract_line(in_file=in_file,
                         out_file=out_file,
                         line=line,
                         gal=gal,
                         chan_width=chan_width,
                         quiet=quiet)            

    return
    
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# Routines to extract continuum from a measurement set
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

def contsub(
    in_file=None,
    lines_to_flag=None,
    gal=None,
    vsys=0.0,
    vwidth=500.,
    quiet=False    
    ):
    """
    Carry out uv continuum subtraction on a measurement set. First
    figures out channels corresponding to spectral lines for a suite
    of bright lines.
    """

    sol_kms = 2.99e5

    # Set up the input file

    if os.path.isdir(in_file) == False:
        if quiet == False:
            print "Input file not found."
        return

    # pull the parameters from the galaxy in the mosaic file

    if gal != None:
        mosaic_parms = read_mosaic_key()
        if mosaic_parms.has_key(gal):
            vsys = mosaic_parms[gal]['vsys']
            vwidth = mosaic_parms[gal]['vwidth']

    # set the list of lines to flag

    if lines_to_flag == None:
        lines_to_flag = line_list.lines_co + line_list.lines_13co + line_list.lines_c18o

    vm = au.ValueMapping(in_file)

    spw_flagging_string = ''
    first = True
    for spw in vm.spwInfo.keys():
        this_spw_string = str(spw)+':0'
        if first:
            spw_flagging_string += this_spw_string
            first = False
        else:
            spw_flagging_string += ','+this_spw_string            

    for line in lines_to_flag:
        rest_linefreq_ghz = line_list.line_list[line]

        shifted_linefreq_hz = rest_linefreq_ghz*(1.-vsys/sol_kms)*1e9
        hi_linefreq_hz = rest_linefreq_ghz*(1.-(vsys-vwidth/2.0)/sol_kms)*1e9
        lo_linefreq_hz = rest_linefreq_ghz*(1.-(vsys+vwidth/2.0)/sol_kms)*1e9

        spw_list = au.getScienceSpwsForFrequency(this_infile,
                                                 shifted_linefreq_hz)
        if spw_list == []:
            continue

        print "Found overlap for "+line
        for this_spw in spw_list:
            freq_ra = vm.spwInfo[this_spw]['chanFreqs']
            chan_ra = np.arange(len(freq_ra))
            to_flag = (freq_ra >= lo_linefreq_hz)*(freq_ra <= hi_linefreq_hz)
            to_flag[np.argmin(np.abs(freq_ra - shifted_linefreq_hz))]
            low_chan = np.min(chan_ra[to_flag])
            hi_chan = np.max(chan_ra[to_flag])                
            this_spw_string = str(this_spw)+':'+str(low_chan)+'~'+str(hi_chan)
            if first:
                spw_flagging_string += this_spw_string
                first = False
            else:
                spw_flagging_string += ','+this_spw_string

    print "... proposed channels to avoid "+spw_flagging_string

    os.system('rm -rf '+in_file+'.contsub')
    uvcontsub(vis=in_file,
              fitspw=spw_flagging_string,
              excludechans=True)

    return

def extract_continuum(
    in_file=None,
    out_file=None,
    lines_to_flag=None,
    gal=None,
    vsys=0.0,
    vwidth=500.,
    do_statwt=True,
    do_collapse=True,
    quiet=False):
    """
    Extract a continuum measurement set, flagging any specified lines,
    reweighting using statwt, and then collapsing to a single "channel
    0" measurement.
    """

    sol_kms = 2.99e5

    # Set up the input file

    if os.path.isdir(in_file) == False:
        if quiet == False:
            print "Input file not found."
        return

    # pull the parameters from the galaxy in the mosaic file

    if gal != None:
        mosaic_parms = read_mosaic_key()
        if mosaic_parms.has_key(gal):
            vsys = mosaic_parms[gal]['vsys']
            vwidth = mosaic_parms[gal]['vwidth']

    # set the list of lines to flag

    if lines_to_flag == None:
        lines_to_flag = line_list.lines_co + line_list.lines_13co + line_list.lines_c18o

    # Make a continuum copy of the data

    os.system('rm -rf '+out_file)
    os.system('rm -rf '+out_file+'.flagversions')

    command = 'cp -r -H '+in_file+' '+out_file
    print command
    var = os.system(command)
    print var
    
    # Figure out the line channels and flag them

    vm = au.ValueMapping(out_file)

    spw_flagging_string = ''
    first = True
    for spw in vm.spwInfo.keys():
        this_spw_string = str(spw)+':0'
        if first:
            spw_flagging_string += this_spw_string
            first = False
        else:
            spw_flagging_string += ','+this_spw_string            

    for line in lines_to_flag:
        rest_linefreq_ghz = line_list.line_list[line]

        shifted_linefreq_hz = rest_linefreq_ghz*(1.-vsys/sol_kms)*1e9
        hi_linefreq_hz = rest_linefreq_ghz*(1.-(vsys-vwidth/2.0)/sol_kms)*1e9
        lo_linefreq_hz = rest_linefreq_ghz*(1.-(vsys+vwidth/2.0)/sol_kms)*1e9

        spw_list = au.getScienceSpwsForFrequency(out_file,
                                                 shifted_linefreq_hz)
        if spw_list == []:
            continue

        print "Found overlap for "+line
        for this_spw in spw_list:
            freq_ra = vm.spwInfo[this_spw]['chanFreqs']
            chan_ra = np.arange(len(freq_ra))
            to_flag = (freq_ra >= lo_linefreq_hz)*(freq_ra <= hi_linefreq_hz)
            to_flag[np.argmin(np.abs(freq_ra - shifted_linefreq_hz))]
            low_chan = np.min(chan_ra[to_flag])
            hi_chan = np.max(chan_ra[to_flag])                
            this_spw_string = str(this_spw)+':'+str(low_chan)+'~'+str(hi_chan)
            if first:
                spw_flagging_string += this_spw_string
                first = False
            else:
                spw_flagging_string += ','+this_spw_string
        
    print "... proposed flagging "+spw_flagging_string

    if spw_flagging_string != '':
        flagdata(vis=out_file,
                 spw=spw_flagging_string,
                 )
        
    if do_statwt:
        print "... deriving emprical weights using STATWT."
        statwt(vis=out_file,
               datacolumn='DATA')

    if do_collapse:
        print "... Collapsing the continuum to a single channel."

        os.system('rm -rf '+out_file+'.temp_copy')
        os.system('rm -rf '+out_file+'.temp_copy.flagversions')

        command = 'mv '+out_file+' '+out_file+'.temp_copy'
        print command
        var = os.system(command)
        print var

        command = 'mv '+out_file+'.flagversions '+out_file+'.temp_copy.flagversions'
        print command
        var = os.system(command)
        print var

        split(vis=out_file+'.temp_copy',
              outputvis=out_file,
              width=10000,
              datacolumn='DATA',
              keepflags=False)        

        os.system('rm -rf '+out_file+'.temp_copy')
        os.system('rm -rf '+out_file+'.temp_copy.flagversions')
        
    return    

def extract_continuum_for_galaxy(   
    gal=None,
    just_proj=None,
    just_ms=None,
    just_array=None,
    lines_to_flag=None,
    ext='',
    do_statwt=True,
    do_collapse=True,
    quiet=False
    ):
    """
    Extract continuum for all data sets for a galaxy. This knows about
    the PHANGS measurement set keys and is specific to our projects.
    """
    
    if gal == None:
        if quiet == False:
            print "Please specify a galaxy."
        return

    ms_key = read_ms_key()

    if ms_key.has_key(gal) == False:
        if quiet == False:
            print "Galaxy "+gal+" not found in the measurement set key."
        return
    gal_specific_key = ms_key[gal]

    # Look up the galaxy specific parameters

    mosaic_parms = read_mosaic_key()
    if mosaic_parms.has_key(gal):
        vsys = mosaic_parms[gal]['vsys']
        vwidth = mosaic_parms[gal]['vwidth']

    # Change to the right directory

    this_dir = dir_for_gal(gal)
    os.chdir(this_dir)

    # Loop over all projects and measurement sets

    for this_proj in gal_specific_key.keys():

        if just_proj != None:
            if type(just_proj) == type([]):
                if just_proj.count(this_proj) == 0:
                    continue
            else:
                if this_proj != just_proj:
                    continue

        proj_specific_key = gal_specific_key[this_proj]
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
            
            in_file = gal+'_'+this_proj+'_'+this_ms+ext+'.ms'
            out_file = gal+'_'+this_proj+'_'+this_ms+'_cont.ms'

            extract_continuum(
                in_file=in_file,
                out_file=out_file,
                lines_to_flag=lines_to_flag,
                gal=gal,
                do_statwt=do_statwt,
                do_collapse=do_collapse)

    return

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# Routines to characterize measurement sets
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

def pick_phangs_cell_and_imsize(
    in_file=None,
    oversamp=5
    ):
    """
    Wraps estimate_cell_and_imsize and also allows our custom
    overrides.
    """

    cell_size_string, x_size_string, y_size_string = \
        estimate_cell_and_imsize(in_file, oversamp)

    override_dict = read_override_image_params()

    # Check for overrides
    if override_dict.has_key(in_file):
        if override_dict.has_key('cell_size'):
            cell_size_string = override_dict[this_vis]['cell_size']
        if override_dict.has_key('x_size'):
            x_size_string = override_dict[this_vis]['x_size']
        if override_dict.has_key('y_size'):
            y_size_string = override_dict[this_vis]['y_size']    

    return cell_size_string, x_size_string, y_size_string

def estimate_cell_and_imsize(
    in_file=None,    
    oversamp=5
    ):
    """
    Pick a cell and image size for a measurement set. Requests an
    oversampling factor, which is by default 5. Will pick a good size
    for the FFT and will try to pick a round number for the cell size.
    """

    if os.path.isdir(in_file) == False:
        print "File not found."
        return
    
    valid_sizes = []
    for ii in range(10):
        for kk in range(3):
            for jj in range(3):
                valid_sizes.append(2**(ii+1)*5**(jj)*3**(kk))
    valid_sizes.sort()
    valid_sizes = np.array(valid_sizes)

    # Cell size implied by baseline distribution from analysis
    # utilities.

    au_cellsize, au_imsize, au_centralField = \
        au.pickCellSize(in_file, imsize=True, npix=oversamp)
    xextent = au_cellsize*au_imsize[0]*1.2
    yextent = au_cellsize*au_imsize[1]*1.2

    # Make the cell size a nice round number

    if au_cellsize < 0.1:
        cell_size = au_cellsize
    if au_cellsize >= 0.1 and au_cellsize < 0.5:
        cell_size = np.floor(au_cellsize/0.05)*0.05
    if au_cellsize >= 0.5 and au_cellsize < 1.0:
        cell_size = np.floor(au_cellsize/0.1)*0.1
    if au_cellsize >= 1.0 and au_cellsize < 2.0:
        cell_size = np.floor(au_cellsize/0.25)*0.25
    if au_cellsize >= 2.0 and au_cellsize < 5.0:
        cell_size = np.floor(au_cellsize/0.5)*0.5
    if au_cellsize >= 5.0:
        cell_size = np.floor(au_cellsize/1.0)*0.5

    # Now make the image size a good number for the FFT

    need_cells_x = xextent / cell_size
    need_cells_y = yextent / cell_size

    cells_x = np.min(valid_sizes[valid_sizes > need_cells_x])
    cells_y = np.min(valid_sizes[valid_sizes > need_cells_y])

    image_size = [int(cells_x), int(cells_y)]
    cell_size_string = str(cell_size)+'arcsec'

    x_size_string = str(image_size[0])
    y_size_string = str(image_size[1])

    return cell_size_string, x_size_string, y_size_string

# TBD: Add the baseline data extractor to make plots (extract_uv_plots.py)

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# Routines to characterize and manipulate cubes
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

def stat_clean_cube(cube_file=None):
    """
    Calculate statistics for an image cube.
    """
    if cube_file == None:
        print "No cube file specified. Returning"
        return
    imstat_dict = imstat(cube_file)
    
    return imstat_dict

def save_copy_of_cube(
    input_root=None,
    output_root=None):
    """
    Copy a cube to a new name. Used to make a backup copy. Overwrites
    the previous cube of that name.
    """

    wipe_cube(output_root)
    
    os.system('cp -r '+input_root+'.image '+output_root+'.image')
    os.system('cp -r '+input_root+'.model '+output_root+'.model')
    os.system('cp -r '+input_root+'.mask '+output_root+'.mask')
    os.system('cp -r '+input_root+'.pb '+output_root+'.pb')
    os.system('cp -r '+input_root+'.psf '+output_root+'.psf')
    os.system('cp -r '+input_root+'.residual '+output_root+'.residual')
    os.system('cp -r '+input_root+'.psf '+output_root+'.weight')
    os.system('cp -r '+input_root+'.residual '+output_root+'.sumwt')

def wipe_cube(
    cube_root=None):
    """
    Wipe files associated with a cube.
    """
    if cube_root == None:
        return
    os.system('rm -rf '+cube_root+'.image')
    os.system('rm -rf '+cube_root+'.model')
    os.system('rm -rf '+cube_root+'.mask')
    os.system('rm -rf '+cube_root+'.pb')
    os.system('rm -rf '+cube_root+'.psf')
    os.system('rm -rf '+cube_root+'.residual')
    os.system('rm -rf '+cube_root+'.weight')
    os.system('rm -rf '+cube_root+'.sumwt')

def replace_cube_with_copy(
    to_root=None,
    from_root=None):
    """
    Replace a cube with a copy.
    """

    wipe_cube(to_root)

    os.system('cp -r '+from_root+'.image '+to_root+'.image')
    os.system('cp -r '+from_root+'.model '+to_root+'.model')
    os.system('cp -r '+from_root+'.mask '+to_root+'.mask')
    os.system('cp -r '+from_root+'.pb '+to_root+'.pb')
    os.system('cp -r '+from_root+'.psf '+to_root+'.psf')
    os.system('cp -r '+from_root+'.residual '+to_root+'.residual')
    os.system('cp -r '+from_root+'.psf '+to_root+'.weight')
    os.system('cp -r '+from_root+'.residual '+to_root+'.sumwt')


# TBD: Add mask reprojection and combination (makeMask)


# TBD: Add export to fits routines (postProcessCubes)

def export_to_fits(
    cube_root=None):
    """
    Export the various products associated with a CASA cube to FITS.
    """
    
    exportfits(imagename=cube_root+'.image',
               fitsimage=cube_root+'.fits',
               velocity=True, overwrite=True, dropstokes=True, 
               dropdeg=True, bitpix=bitpix)

    exportfits(imagename=cube_root+'_dirty.image',
               fitsimage=cube_root+'_dirty.fits',
               velocity=True, overwrite=True, dropstokes=True, 
               dropdeg=True, bitpix=bitpix)

    exportfits(imagename=cube_root+'.model',
               fitsimage=cube_root+'_model.fits',
               velocity=True, overwrite=True, dropstokes=True, 
               dropdeg=True, bitpix=bitpix)

    exportfits(imagename=cube_root+'.residual',
               fitsimage=cube_root+'_residual.fits',
               velocity=True, overwrite=True, dropstokes=True, 
               dropdeg=True, bitpix=bitpix)

    exportfits(imagename=cube_root+'.mask',
               fitsimage=cube_root+'_mask.fits',
               velocity=True, overwrite=True, dropstokes=True, 
               dropdeg=True, bitpix=bitpix)
    
    exportfits(imagename=cube_root+'.pb',
               fitsimage=cube_root+'_pb.fits',
               velocity=True, overwrite=True, dropstokes=True, 
               dropdeg=True, bitpix=bitpix)

    return


# TBD: Add a routine to actually write the feathering scripts? (feather_script_12m and feather_script_7m)

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# Routines to image the data
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

def call_clean(
    vis = None,
    image_root = None,
    phase_center = "",
    image_size = None,
    cell_size = None,
    restfreq_ghz = -1.0,
    calcres=True,
    calcpsf=True,
    specmode = 'cube',
    deconvolver = 'hogbom',
    threshold = '0.0mJy/beam',
    scales = [0],
    smallscalebias = 0.9,    
    briggs_weight = 0.5,
    niter = 0,
    cycle_niter = 200,
    minpsffraction = 0.5,
    pb_limit = 0.25,
    uv_taper_string = '',
    restoringbeam = '',
    usemask='user',
    mask=0.0,
    interactive = False,
    reset = False,
    logfile = None,
    ):
    """
    Tclean wrapper with sane defaults. Wrapped by our other calls.
    """
    
    if vis == None:
        print "No visibility. Returning."
        return    

    if os.path.isdir(vis) == False:
        print "Visibility file not found. Returning."
        return

    if restfreq_ghz < 0:
        restfreq_str = ''
    else:
        restfreq_str = str(restfreq_ghz)+'GHz'

    if logfile != None:
        oldlogfile = casalog.logfile()
        casalog.setlogfile(logfile)

    if reset:
        wipe_cube(image_root)

    tclean(vis=vis,
           imagename=image_root,
           # Spatial axes
           phasecenter=phase_center,
           cell=cell_size,
           imsize=image_size,
           gridder='mosaic',
           # Spectral axis
           specmode=specmode,
           restfreq=restfreq_str,
           outframe='lsrk',
           veltype='radio',
           # Workflow
           calcres=calcres,
           calcpsf=calcpsf,
           # Deconvolver
           deconvolver=deconvolver,
           scales=scales,
           smallscalebias=smallscalebias,
           pblimit=pb_limit,
           normtype='flatnoise',
           # Restoring beam
           restoringbeam=restoringbeam,
           # U-V plane gridding
           weighting='briggs',
           robust=briggs_weight,
           uvtaper=uv_taper_string,
           # Stopping criterion
           niter=niter,
           threshold=threshold,
           cycleniter=cycle_niter,
           cyclefactor=3.0,
           minpsffraction=minpsffraction,
           # Mask
           usemask=usemask,
           mask=mask,
           pbmask=pb_limit,
           # UI
           interactive=False,
           )

    if logfile != None:
        casalog.setlogfile(oldlogfile)



# TBD: Add dirty map creation (imageMultiscale and imageMultiscale2)

# TBD: Add multiscale imaging (imageMultiscale and imageMultiscale2)

# TBD: Add single scale imaging (imageMultiscale and imageMultiscale2)

# TBD: Add the "recipe" level routines used in the actual imaging (phangsImagingPipeline and phangsImagingPipeline2)

