"""
Standalone routines to analyze and manipulate visibilities.
"""

#region Imports and definitions

import os, sys, re, shutil
import numpy as np
import pyfits # CASA has pyfits, not astropy
import glob

import logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

# Analysis utilities
import analysisUtils as au

# CASA stuff
import casaStuff

# Pipeline versionining
from pipelineVersion import version as pipeVer

#endregion

#region Routines for basic characterization

#endregion

#region Routines to analyze and extract lines in measurement sets

def split_science_targets(
    in_ms, 
    out_ms, 
    do_split = True, 
    do_statwt = False, 
    use_symlink = True, 
    overwrite = False, 
    quiet = False, 
    ):
    """
    Split science targets from the input ALMA measurement set to a new measurement set.
    
    Setting do_split = False will just make a copy of the original measurement set without splitting science targets data.
    
    Setting use_symlink = False will make a copy of the original measurement set so as to avoid any touching of the original data. 
    
    Setting overwrite = True will overwrite existing output data. 
    """
    
    # 
    # This funcion was part of the copy_data() function in older version phangsPipeline.py.
    # 
    
    # 
    # check input ms data dir
    if os.path.isdir(in_ms):
        in_file = in_ms
    else:
        logger.error('Error! The input uv data measurement set "'+in_ms+'"does not exist!')
        raise Exception('Error! The input uv data measurement set "'+in_ms+'"does not exist!')
    # 
    # check output suffix
    if not re.match(r'^(.*)\.ms$', out_ms, re.IGNORECASE):
        out_file = out_ms + '.ms'
    else:
        out_file = out_ms
    # 
    # check existing copied data in the imaging directory
    if os.path.isdir(out_file):
        if not overwrite:
            logger.warning('Found existing copied data '+out_file+', will not re-copy it.')
            return
        else:
            shutil.rmtree(out_file)
            if os.path.isdir(out_file+'.flagversions'):
                shutil.rmtree(out_file+'.flagversions')
    # 
    # prepare the copied data folder name (copied_file is a data folder).
    # If we are going to do some additional processing, make this an intermediate file ("_copied"). 
    if do_split:
        copied_file = re.sub(r'^(.*)\.ms$', r'\1_copied.ms', out_file, re.IGNORECASE)
    else:
        copied_file = out_file
        use_symlink = False
    # 
    # Copy. We could place a symbolic link here using ln -s
    # instead, but instead I think the right move is to make
    # the intermediate files and then clean them up. This
    # avoids "touching" the original data at all.
    if use_symlink:
        command = 'ln -fsT '+in_file+' '+copied_file
        logger.info(command)
        var = os.system(command)
        logger.info(var)
        
        command = 'ln -fsT '+in_file+'.flagversions'+' '+copied_file+'.flagversions'
        logger.info(command)
        var = os.system(command)
        logger.info(var)
        
    else:
        command = 'cp -Lr '+in_file+' '+copied_file
        logger.info(command)
        var = os.system(command)  
        logger.info(var)
        
        command = 'cp -Lr '+in_file+'.flagversions'+' '+copied_file+'.flagversions'
        logger.info(command)
        var = os.system(command)
        logger.info(var)
    # 
    # check copied_file, make sure copying was done
    if not os.path.isdir(copied_file):
        logger.error('Failed to copy the uv data to '+os.path.abspath(out_file)+'! Please check your file system writing permission!')
        raise Exception('Failed to copy the uv data to the imaging directory! Please check your file system writing permission!')
    # 
    # call CASA split
    if do_split:
        # 
        if not quiet:
            logger.info("Splitting out science target data.")
        # 
        # If present, we use the corrected column. If not,
        # then we use the data column.
        #mytb = au.createCasaTool(tbtool)
        #mytb.open(copied_file)
        #colnames = mytb.colnames()
        casaStuff.tb.open(copied_file, nomodify = True)
        colnames = casaStuff.tb.colnames()
        if 'CORRECTED_DATA' in colnames:
            logger.info("Data has a CORRECTED column. Will use that.")
            use_column = 'CORRECTED'
        else:
            logger.info("Data lacks a CORRECTED column. Will use DATA column.")
            use_column = 'DATA'
        casaStuff.tb.close()
        # 
        # 
        casaStuff.split(vis = copied_file, 
                        intent = 'OBSERVE_TARGET#ON_SOURCE', 
                        datacolumn = use_column, 
                        outputvis = out_file)
        # 
        # 
        os.system('rm -rf '+copied_file)
        os.system('rm -rf '+copied_file+'.flagversions')
    # 
    # call CASA statwt
    if do_statwt:
        # 
        if not quiet:
            logger.info("Using statwt to re-weight the data.")
        # 
        casaStuff.statwt(vis = out_file, 
                         datacolumn = 'DATA')
    # 
    # end of split_science_targets()



def list_lines_in_ms(
    in_file= None,
    vsys=0.0,
    gal=None,
    quiet=False,
    ):    
    """
    List the lines likely to be present in a measurement set. This can
    be a general purpose utility.
    """

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
                 chan_fine=0.5,
                 rebin_factor=5,
                 do_statwt=False,
                 edge_for_statwt=-1,
                 quiet=False):
    """
    Extract a spectral line from a measurement set and regrid onto a
    new velocity grid with the desired spacing. This doesn't
    necessarily need the PHANGS keys in place and may be a general
    purpose utility. There are some minor subtleties here related to
    regridding and rebinning.
    """

    if quiet == False:
        print "--------------------------------------"
        print "EXTRACT_LINE begins:"

    # pull the parameters from the galaxy in the mosaic file. This is
    # PHANGS-specific. Just ignore the gal keyword to use the routine
    # for non-PHANGS applications.

    if gal != None:
        mosaic_parms = read_mosaic_key()
        if mosaic_parms.has_key(gal):
            vsys = mosaic_parms[gal]['vsys']
            vwidth = mosaic_parms[gal]['vwidth']

    # Set up the input file

    if os.path.isdir(in_file) == False:
        if quiet == False:
            print "... input file not found."
        return

    # Look up the line

    if line_list.line_list.has_key(line) == False:
        if quiet == False:
            print "... line not found. Give lower case abbreviate found in line_list.py"
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
            print "... no spectral windows contain this line at this redshift."
        return

    if quiet == False:
        print "... spectral windows to consider: "+spw_list_string

    # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
    # STEP 1. Shift the zero point AND change the channel width (slightly).
    # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

    start_vel_kms = (vsys - vwidth/2.0)
    chan_width_hz = au.getChanWidths(in_file, spw_list_string)
    current_chan_width_kms = abs(chan_width_hz / (restfreq_ghz*1e9)*sol_kms)        
    if chan_fine == -1:
        nchan_for_recenter = int(np.max(np.ceil(vwidth / current_chan_width_kms)))
    else:
        nchan_for_recenter = int(np.max(np.ceil(vwidth / chan_fine)))

    # Cast to text with specified precision.
    restfreq_string = "{:12.8f}".format(restfreq_ghz)+'GHz'
    start_vel_string =  "{:12.8f}".format(start_vel_kms)+'km/s'
    chanwidth_string =  "{:12.8f}".format(chan_fine)+'km/s'

    if quiet == False:
        print "... shifting the fine grid (before any regridding)"
        print "... rest frequency: "+restfreq_string
        print "... new starting velocity: "+start_vel_string
        print "... original velocity width: "+str(current_chan_width_kms)
        print "... target velocity width: "+str(chan_fine)
        print "... number of channels at this stage: "+str(nchan_for_recenter)

    os.system('rm -rf '+out_file+'.temp')
    os.system('rm -rf '+out_file+'.temp.flagversions')
    if chan_fine == -1:
        mstransform(vis=in_file,
                    outputvis=out_file+'.temp',
                    spw=spw_list_string,
                    datacolumn='DATA',
                    combinespws=False,
                    regridms=True,
                    mode='velocity',
                    #interpolation='linear',
                    interpolation='cubic',
                    start=start_vel_string,
                    nchan=nchan_for_recenter,
                    restfreq=restfreq_string,
                    outframe='lsrk',
                    veltype='radio',
                    )
    else:
        mstransform(vis=in_file,
                    outputvis=out_file+'.temp',
                    spw=spw_list_string,
                    datacolumn='DATA',
                    combinespws=False,
                    regridms=True,
                    mode='velocity',
                    #interpolation='linear',
                    interpolation='cubic',
                    start=start_vel_string,
                    nchan=nchan_for_recenter,
                    restfreq=restfreq_string,
                    width=chanwidth_string,
                    outframe='lsrk',
                    veltype='radio',
                    )

    current_file = out_file+'.temp'

    # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
    # STEP 2. Change the channel width by integer binning. 
    # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

    if quiet == False:
        print "... channel averaging"
        print "... rebinning factor: "+str(rebin_factor)

    if rebin_factor > 1:
        os.system('rm -rf '+out_file+'.temp2')
        os.system('rm -rf '+out_file+'.temp2.flagversions')
        mstransform(vis=current_file,
                    outputvis=out_file+'.temp2',
                    spw='',
                    datacolumn='DATA',
                    regridms=False,
                    chanaverage=True,
                    chanbin=rebin_factor,
                    )
        current_file = out_file+'.temp2'

    # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
    # STEP 3. Combine the SPWs
    # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

    if quiet == False:
        print "... combining spectral windows"

    os.system('rm -rf '+out_file)
    os.system('rm -rf '+out_file+'.flagversions') 
    mstransform(vis=current_file,
                outputvis=out_file,
                spw='',
                datacolumn='DATA',
                regridms=False,
                chanaverage=False,
                combinespws=True
                )    

    if quiet == False:
        print "... deleting old files"
        
    # Clean up
    os.system('rm -rf '+out_file+'.temp')
    os.system('rm -rf '+out_file+'.temp.flagversions')
    os.system('rm -rf '+out_file+'.temp2')
    os.system('rm -rf '+out_file+'.temp2.flagversions')

    # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
    # STEP 4. Re-weight the data using statwt
    # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

    # N.B. need the sliding time bin to make statwt work.

    if do_statwt:

        if edge_for_statwt == -1:
            exclude_str = ''
        else:
            nchan_final = int(np.floor(nchan_for_recenter / rebin_factor)+1)
            exclude_str = '*:'+str(edge_for_statwt-1)+'~'+\
                str(nchan_final-(edge_for_statwt-2))

        print "... running statwt with exclusion: "+exclude_str

        # This needs to revert to oldstatwt, it seems not to work in the new form

        test = statwt(vis=out_file,
                      timebin='0.001s',
                      slidetimebin=False,
                      chanbin='spw',
                      statalg='classic',
                      datacolumn='data',
                      excludechans=exclude_str,
                      )

    print "--------------------------------------"

    return

#endregion

#region Continuum-related functions

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
            print("Input file not found: "+in_file)
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
        
    # Here - this comman needs to be examined and refined in CASA
    # 5.6.1 to see if it can be sped up. Right now things are
    # devastatingly slow.
    if do_statwt:
        print "... deriving empirical weights using STATWT."
        statwt(vis=out_file,
               timebin='0.001s',
               slidetimebin=False,
               chanbin='spw',
               statalg='classic',
               datacolumn='data',
               )

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

#endregion

#region Measurement set characterization, e.g., noise


def noise_spectrum(
    vis=None,
    stat_name="medabsdevmed",
    start_chan=None,
    stop_chan=None):
    """
    Calculates the u-v based noise spectrum and returns it as an array.
    """
    if vis == None:
        return None
    
   # Note the number of channels in SPW 0

    vm = au.ValueMapping(vis)
    nchan = vm.spwInfo[0]['numChannels']
    spec = np.zeros(nchan)
    for ii in range(nchan):
        if start_chan != None:
            if ii < start_chan:
                continue
        if stop_chan != None:
            if ii > stop_chan:
                continue
        print "Channel "+str(ii)+" / "+str(nchan)
        result = visstat(vis=vis,
                         axis='amp',
                         spw='0:'+str(ii),
                         )
        if result == None:
            print "Skipping channel."
            continue
        spec[ii] = result[result.keys()[0]][stat_name]
        
    return spec

#endregion
