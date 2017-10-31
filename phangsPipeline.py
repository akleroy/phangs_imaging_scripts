# NB: CASA doesn't always include the pwd in the module search path.

import os
import numpy as np

# Other PHANGS scripts
import line_list

# CASA imports
import analysisUtils as au
from split import split
from statwt import statwt
from mstransform import mstransform
from concat import concat
from uvcontsub import uvcontsub

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# Interface to the text file keys that steer the process
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

def read_ms_key(fname='../scripts/ms_file_key.txt'):
    """
    Read the measurement set key into a big dictionary.
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
    Read the directory key, which gives us a general way to sort out which MS files go in which directory.
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
    Return the working directory given a galaxy name.
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
    Read the file containing the centers and velocities for each source.
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

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# Routines to move data around
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

def copy_data(gal=None,
              just_proj=None,
              just_ms=None,
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

            if do_split == False:
                continue

            out_file = gal+'_'+this_proj+'_'+this_ms+'.ms'

            os.system('rm -rf '+out_file)
            os.system('rm -rf '+out_file+'.flagversions')

            split(vis=copied_file
                  , intent ='OBSERVE_TARGET#ON_SOURCE'
                  , datacolumn='DATA'
                  , outputvis=out_file)        

            os.system('rm -rf '+copied_file)
            os.system('rm -rf '+copied_file+'.flagversions')

            if do_statwt == False:
                continue

            statwt(vis=out_file,
                   datacolumn='DATA')

    if quiet ==False:
        print "--------------------------------------------------------"
        print "END: Copying data from original location."
        print "--------------------------------------------------------"

def concat_ms_for_line(gal=None,
                       just_proj=None,
                       just_ms=None,
                       line='co21',
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

    # TBD

    # Concatenate all of the relevant files

    out_file =  gal+'_'+tag+'_'+line+'.ms'

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
    
    chan0_vis = gal+'_'+tag+'_'+line+'_chan0.ms'

    os.system('rm -rf '+chan0_vis)
    os.system('rm -rf '+chan0_vis+'.flagversions')
    split(vis=out_file
          , datacolumn='DATA'
          , spw=''
          , outputvis=chan0_vis
          , width=10000)

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# Routines to extract lines from a measurement set
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

def list_lines_in_ms(
    in_file= None,
    vsys=0.0,
    gal=None,
    ):    
    """
    List the lines likely to be present in a measurement set.
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
    overlap a line.
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
    new velocity grid with the desired spacing.
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
    target_width_hz = chan_width/sol_kms*restfreq_ghz*1e9
    rebin_factor = min(target_width_hz / chan_width_hz)
    
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
    line='co21',
    vsys=0.0,
    vwidth=500.,
    chan_width=2.5,    
    ext='',
    quiet=False
    ):
    """
    Extract a given line for all data sets for a galaxy.
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

def extract_phangs_lines(   
    gal=None,
    just_proj=None,
    just_ms=None,
    quiet=False
    ):
    """
    Extract all phangs lines and continuum for a galaxy.
    """

    extract_line_for_galaxy(   
        gal=gal,
        line='co21',
        chan_width=2.5,    
        quiet=quiet
        )

    extract_line_for_galaxy(   
        gal=gal,
        line='c18o21',
        chan_width=6.0, 
        quiet=quiet
        )
    
    # could add sio5to4 once I have some more disk space

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
        rest_linefreq_ghz = line_list[line]

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
    quiet=False):
    """
    Extract a continuum measurement set, flagging any specified lines.
    """
    pass

def concat_continuum(
    ):
    """
    Concatenate continuum data sets.
    """
    pass
    
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# Routines to image the data
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
