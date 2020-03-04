"""
Standalone routines to analyze and manipulate visibilities.
"""

# 20200226: introduced os.mkdir(outfile+'.touch') os.rmdir(outfile+'.touch') to make sure we can handle sudden system break.

#region Imports and definitions

import os, sys, re, shutil, inspect, copy
import numpy as np
#import pyfits # CASA has pyfits, not astropy
#import glob

import logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

# Analysis utilities
import analysisUtils as au

# CASA stuff
import casaStuff

# Spectral lines
import line_list

# Pipeline versionining
from pipelineVersion import version as pipeVer

#endregion

#region Routines for basic characterization

#endregion

#region Routines to analyze and extract lines in measurement sets

# Physical constants
sol_kms = 2.9979246e5

#######################################
# Split and combine measurement sets. #
#######################################
    
def split_science_targets(
    infile = None, 
    outfile = None, 
    do_split = True,
    split_field = '', 
    split_intent = 'OBSERVE_TARGET#ON_SOURCE', 
    split_spw = '',
    do_statwt = False, 
    use_symlink = True, 
    overwrite = False, 
    ):
    """
    Split science targets from the input ALMA measurement set to form
    a new, science-only measurement set. Optionally reweight the data
    using statwt. If only a copy is requested, use a symlink.
    
    Args:

    infile (str): The input measurement set data.
    
    outfile (str): The output measurement set data.

    field (str): The field used for selection.
        
    do_split (bool): Set to False to only make a copy of the original
    measurement set without splitting science targets data. The
    default is True, splitting science targets data.
        
    use_symlink (bool): Set to True to use a symlink instead of making
    an actual copy of the original measurement set. Only works with
    do_split=False
        
    overwrite (bool): Set to True to overwrite existing output
    data. The default is False, not overwriting anything.
    
    Inputs:

    infile: ALMA measurement set data folder.
    
    Outputs:

    outfile: ALMA measurement set data folder.
    
    """
        
    # Check inputs

    if infile is None:
        logging.error("Please specify infile.")
        raise Exception("Please specify infile.")
    
    if outfile is None:
        logging.error("Please specify outfile.")
        raise Exception("Please specify outfile.")
    
    if not os.path.isdir(infile):
        logger.error('Error! The input uv data measurement set "'+infile+'"does not exist!')
        raise Exception('Error! The input uv data measurement set "'+infile+'"does not exist!')    
     
    if do_statwt and not do_split:
        logging.warning("Stat weight option (do_statwt) only enabled with do_split.")

    # Check for presence of existing outfile and abort if it is found
    # without overwrite permission.

    if os.path.isdir(outfile) and not os.path.isdir(outfile+'.touch'):
        if not overwrite:
            logger.warning('Found existing data "'+outfile+'", will not overwrite.')
            return()

    # Delete existing output data.

    for suffix in ['', '.flagversions', '.touch']:
        if os.path.islink(outfile+suffix):
            os.unlink(outfile+suffix)
            logger.debug('os.unlink "'+outfile+'"')
            
        if os.path.isdir(outfile+suffix):
            shutil.rmtree(outfile+suffix)

    ###########################################################################
    # COPY: If split is turned off, we are just making a copy.
    ###########################################################################

    if do_split == False:

        logger.info("")
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%")
        logger.info("I will just copy or link the data.")
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%")
        logger.info("")

        if use_symlink:
            
            # Make links

            if os.path.isdir(infile):
                os.symlink(infile, outfile)
                logger.debug('os.symlink "'+infile+'", "'+outfile+'"')

            if os.path.isdir(infile+'.flagversions'):
                os.symlink(infile+'.flagversions', outfile+'.flagversions')
                logger.debug('os.symlink "'+infile+'.flagversions'+'", "'+outfile+'.flagversions'+'"')

            # Check

            if not os.path.islink(outfile):
                logger.error('Failed to link the uv data to '+os.path.abspath(outfile)+'!')
                logger.error('Please check your file system writing permission or system breaks.')
                raise Exception('Failed to link the uv data to the imaging directory.')

            return()

        else:

            # Check existing output data

            has_existing_outfile = False
            if os.path.isdir(outfile) and not os.path.isdir(outfile+'.touch'):
                if not overwrite:
                    has_existing_outfile = True

            # delete existing copied data if not overwriting

            if not has_existing_outfile:
                for suffix in ['', '.flagversions', '.touch']:
                    if os.path.isdir(outfile+suffix):
                        shutil.rmtree(outfile+suffix)
                        logger.debug('shutil.rmtree "'+outfile+suffix+'"')
                        
            # copy the data (.touch directory is a temporary flagpost)

            os.mkdir(outfile+'.touch')

            if os.path.isdir(infile):
                shutil.copytree(infile, outfile)
                logger.debug('shutil.copytree "'+infile+'", "'+outfile+'"')

            if os.path.isdir(infile+'.flagversions'):
                shutil.copytree(infile+'.flagversions', outfile+'.flagversions')
                logger.debug('shutil.copytree "'+infile+'.flagversions'+'", "'+outfile+'.flagversions'+'"')

            os.rmdir(outfile+'.touch')            
             
            # check copied_file, make sure copying was done

            if not os.path.isdir(outfile) or \
                    os.path.isdir(outfile+'.touch'):
                logger.error('Failed to copy the uv data to '+os.path.abspath(oufile)+'!')
                logger.error('Please check your file system writing permission or system breaks.')
                raise Exception('Failed to copy the uv data to the imaging directory.')

            return()

    ###########################################################################
    # SPLIT: If split is turned on, split and extract the science target
    ###########################################################################

    if do_split:

        logger.info("")
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%")
        logger.info("I will split out the data.")
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%")
        logger.info("")
        
        logger.info('Splitting from '+infile+' to '+outfile)

        # Verify the column to use. If present, we use the corrected
        # column. If not, then we use the data column.

        casaStuff.tb.open(infile, nomodify = True)
        colnames = casaStuff.tb.colnames()
        if 'CORRECTED_DATA' in colnames:
            logger.info("Data has a CORRECTED column. Will use that.")
            use_column = 'CORRECTED'
        else:
            logger.info("Data lacks a CORRECTED column. Will use DATA column.")
            use_column = 'DATA'
        casaStuff.tb.close()

        logger.debug('field='+split_intent)
        logger.debug('intent='+split_field)

        os.mkdir(outfile+'.touch')
        casaStuff.split(vis = infile, 
                        intent = split_intent, 
                        field = split_field,
                        spw = split_spw,
                        datacolumn = use_column, 
                        outputvis = outfile)
        os.rmdir(outfile+'.touch')

        # Re-weight the data if desired.

        if do_statwt:
            logger.info("Using statwt to re-weight the data.")
            logger.debug('casa statwt vis="'+outfile+'"')
            os.mkdir(outfile+'.touch')
            casaStuff.statwt(vis = outfile, 
                             datacolumn = 'DATA')
            os.rmdir(outfile+'.touch')
            
        return()

    return()

def concat_ms(
    infile_list=None,  
    outfile=None,  
    freqtol='',
    dirtol='',
    overwrite = False, 
    ):
    """
    Concatenate a list of measurement sets into one measurement set. A
    thin wrapper to concat.
    
    Args:
        infile_list (list or str): The input list of measurement sets. 
        outfile (str): The output measurement set data with suffix ".ms".
    
    Inputs:
        infile: ALMA measurement set data folder.
    
    Outputs:
        outfile: ALMA measurement set data folder.
    
    """
    # Check inputs

    if infile_list is None:
        logging.error("Please specify infile_list.")
        raise Exception("Please specify infile_list.")
    
    if outfile is None:
        logging.error("Please specify outfile.")
        raise Exception("Please specify outfile.")

    # make sure the input infile_list is a list
    if np.isscalar(infile_list):
        infile_list = [infile_list]
    
    # check file existence
    for this_infile in infile_list:
        if not os.path.isdir(this_infile):
            logger.error('Error! The input measurement set "'+this_infile+'" not found')
            raise Exception('Error! The input measurement set "'+this_infile+'" not found')
    
    # PROPOSE TO DEPRECATE THIS (e.g., .contsub violates this rule)
    # check output suffix
    #if re.match(r'^(.*)\.ms$', outfile, re.IGNORECASE):
    #    out_name = re.sub(r'^(.*)\.ms$', r'\1', outfile, re.IGNORECASE)
    #    outfile = out_name + '.ms'
    #else:
    #    out_name = outfile
    #    outfile = out_name + '.ms'
    
    # Quit if output data are present and overwrite is off.
    if os.path.isdir(outfile) and not os.path.isdir(outfile+'.touch'):
        if not overwrite:
            logger.warning('Found existing output data "'+outfile+'", will not overwrite it.')
            return()

    # if overwrite or no file present, then delete existing output data.
    for suffix in ['', '.flagversions', '.touch']:
        if os.path.isdir(outfile+suffix):
            shutil.rmtree(outfile+suffix)

    # Concatenate all of the relevant files
    os.mkdir(outfile+'.touch')
    casaStuff.concat(vis = infile_list, 
                     concatvis = outfile,
                     freqtol=freqtol, dirtol=dirtol)
                     #<TODO># what about freqtol? set as an input? dirtol?
    os.rmdir(outfile+'.touch')

    return()

##########################
# Continuum subtraction. #
##########################

def contsub(
    infile = None, 
    outfile = None,
    lines_to_exclude = None, 
    vsys = None, 
    vwidth = None, 
    overwrite = False, 
    ):
    """
    Carry out uv continuum subtraction on a measurement set. First
    figures out channels corresponding to spectral lines for a
    provided suite of bright lines.
    """
    
    if infile is None:
        logging.error("Please specify infile.")
        raise Exception("Please specify infile.")
    
    if outfile is None:
        outfile = infile+'.contsub'

    if not os.path.isdir(infile):
        logger.error('The input uv data measurement set "'+infile+'"does not exist.')
        raise Exception('The input uv data measurement set "'+infile+'"does not exist.')

    # check existing output data in the imaging directory
    if os.path.isdir(infile+'.contsub') and not os.path.isdir(infile+'.contsub'+'.touch'):
        if not overwrite:
            logger.warning('Found existing output data "'+infile+'.contsub'+'", will not overwrite it.')
            return
    if os.path.isdir(infile+'.contsub'):
        shutil.rmtree(infile+'.contsub')
    if os.path.isdir(infile+'.contsub'+'.touch'):
        shutil.rmtree(infile+'.contsub'+'.touch')

    # Figure out which channels to exclude from the fit.

    # find_spw_channels_for_lines_to_flag
    spw_flagging_string = find_spw_channels_for_lines_to_flag(infile = infile, 
                                                              lines_to_flag = lines_to_flag, 
                                                              vsys = vsys, 
                                                              vwidth = vwidth)
    # 
    # uvcontsub, this outputs infile+'.contsub'
    os.mkdir(infile+'.contsub'+'.touch')
    casaStuff.uvcontsub(vis = infile,
                        fitspw = spw_flagging_string,
                        excludechans = True,
                        combine='spw',
                        )
                        # Exception: Error in uvcontsub: combine must include 'spw' when the fit is being applied to spws outside fitspw.

    os.rmdir(infile+'.contsub'+'.touch')
    # 
    return

##########################################################
# Interface between spectral lines and spectral windows. #
##########################################################

def find_spws_for_line(
    infile = None, 
    line = None, 
    vsys = 0.0, 
    vwidth = 0.0, 
    exit_on_error = True, 
    ):
    """
    List the spectral windows in the input ms data that contains the
    input line, given the line velocity (vsys) and line width
    (vwidth), which are in units of km/s. Defaults to rest frequency
    (with vsys = 0.0 and vwidth = 0.0).
    """

    # Check inputs
    
    if infile is None:
        logging.error("Please specify infile.")
        raise Exception("Please specify infile.")

    if line is None:
        logging.error("Please specify an input line.")
        raise Exception("Please specify an input line.")
        
    # Verify file existence

    if not os.path.isdir(infile):
        logger.error('Error! The input uv data measurement set "'+infile+'"does not exist!')
        raise Exception('Error! The input uv data measurement set "'+infile+'"does not exist!')
 
    # Get the line name and rest-frame frequency in the line_list
    # module for the input line

    line_name, restfreq_ghz = line_list.get_line_name_and_frequency(line, exit_on_error=exit_on_error)

    # Work out the frequencies at the line edes.

    target_freq_ghz = restfreq_ghz*(1.-vsys/sol_kms)
    target_freq_high = restfreq_ghz*(1.-(vsys-0.5*vwidth)/sol_kms)
    target_freq_low = restfreq_ghz*(1.-(vsys+0.5*vwidth)/sol_kms)

    # Work out which spectral windows contain the line contain

    spw_list = []
    
    # Loop over line edges and center frequencies and find the match
    # to the spectral window (spw) for each. In theory, this could be
    # improved to deal with the case where an intermediate spectral
    # window is TOTALLY encompassed between the low/high and center.

    for target_freq in [target_freq_high, target_freq_ghz, target_freq_low]:

        # run analysisUtils.getScienceSpwsForFrequency() to get the corresponding spectral window (spw)
        this_spw_list = au.getScienceSpwsForFrequency(infile, target_freq*1e9)    
        spw_list.extend(this_spw_list)

    # If we don't find the line in this data set, issue a warning and
    # return.

    if len(spw_list) == 0:
        logger.warning('No spectral windows contain the input line.')
        spw_list_string = None
    else:
        
        # sort and remove duplicates
        spw_list = sorted(list(set(spw_list)))

        # make spw_list_string appropriate for use in selection
        spw_list_string = ','.join(np.array(spw_list).astype(str))
     
    return(spw_list_string)

def find_spw_channels_for_lines_to_flag(
    infile, 
    lines_to_flag, 
    vsys = None, 
    vwidth = None, 
    ):
    """
    List the spectral window and channels corresponding to the input lines in the input ms data. 
    Galaxy system velocity (vsys) and velocity width (vwidth) in units of km/s are needed. 
    """
    
    # 
    # check input ms data dir
    if not os.path.isdir(infile):
        logger.error('Error! The input uv data measurement set "'+infile+'"does not exist!')
        raise Exception('Error! The input uv data measurement set "'+infile+'"does not exist!')
    # 
    # check vsys
    if vsys is None:
        logger.error('Error! Please input vsys for the galaxy systematic velocity in units of km/s.')
        raise Exception('Error! Please input vsys.')
    # 
    # check vwidth
    if vwidth is None:
        logger.error('Error! Please input vwidth for the galaxy systematic velocity in units of km/s.')
        raise Exception('Error! Please input vwidth.')
    # 
    # set the list of lines to flag
    if lines_to_flag is None:
        lines_to_flag = line_list.line_families['co'] + line_list.line_families['13co'] + line_list.line_families['c18o']
    else:
        lines_to_flag_copied = copy.copy(lines_to_flag)
        lines_to_flag = []
        for line_to_flag_copied in lines_to_flag_copied:
            matched_line_names = line_list.get_line_names_in_line_family(line_to_flag_copied, exit_on_error = False)
            if len(matched_line_names) > 0:
                lines_to_flag.extend(matched_line_names)
            else:
                matched_line_name, matched_line_freq = line_list.get_line_name_and_frequency(line_to_flag_copied, exit_on_error = False)
                if matched_line_name is not None:
                    lines_to_flag.append(matched_line_name)
    # 
    if len(lines_to_flag) == 0:
        logger.debug('No line to flag for the continuum .')
    else:
        logger.debug('Lines to flag for the continuum: '+str(lines_to_flag))
    # 
    # Figure out the line channels and flag them
    vm = au.ValueMapping(infile)
    # 
    spw_flagging_string = ''
    first = True
    for spw in vm.spwInfo.keys():
        this_spw_string = str(spw)+':0' #<TODO># flag the first channel for all spws
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
        
        spw_list = au.getScienceSpwsForFrequency(infile,
                                                 shifted_linefreq_hz)
        if len(spw_list) == 0:
            continue
        
        for this_spw in spw_list:
            freq_ra = vm.spwInfo[this_spw]['chanFreqs']
            chan_ra = np.arange(len(freq_ra))
            to_flag = (freq_ra >= lo_linefreq_hz)*(freq_ra <= hi_linefreq_hz)
            to_flag[np.argmin(np.abs(freq_ra - shifted_linefreq_hz))]
            low_chan = np.min(chan_ra[to_flag])
            hi_chan = np.max(chan_ra[to_flag])                
            this_spw_string = str(this_spw)+':'+str(low_chan)+'~'+str(hi_chan)
            logger.info("... found line "+line+" in spw "+this_spw_string)
            if first:
                spw_flagging_string += this_spw_string
                first = False
            else:
                spw_flagging_string += ','+this_spw_string
    
    logger.info("... proposed line channels to flag "+spw_flagging_string)
    
    return spw_flagging_string

def compute_chanwidth_for_line(
    infile, 
    line, 
    vsys = None, 
    vwidth = None, 
    ): 
    """Calculates the coarsest channel width among all spectral windows in the input measurement set that contain the input line.
    
    Args:
        infile (str): The input measurement set data with suffix ".ms".
        line (str): Line name. 
        output_spw (bool): Set to True to output not only found line names but also corresponding spectral window (spw) number. 
    
    Returns:
        chan_width_kms (float): 
    
    """
    
    # 
    # This funcion is modified from the chanwidth_for_line() function in older version phangsPipeline.py.
    # 
    
    # 
    # Get the line name and line center rest-frame frequency in the line_list module for the input line
    line_name, restfreq_ghz = line_list.get_line_name_and_frequency(line, exit_on_error=True) # exit_on_error = True
    # 
    # Find spws for line
    spw_list_string = find_spws_for_line(infile, line, vsys = vsys, vwidth = vwidth)
    # 
    # Figure out how much averaging is needed to reach the target resolution
    chan_width_hz = au.getChanWidths(infile, spw_list_string)
    # 
    # Convert to km/s and return
    chan_width_kms = abs(chan_width_hz / (restfreq_ghz*1e9)*sol_kms) #<TODO># here we use the rest-frame frequency to compute velocities - Agree we need to fix this.
    # 
    # Return
    return chan_width_kms

#################################################################
# Extract single-line or continuum data from a measurement set. #
#################################################################

def extract_line(
    infile, 
    outfile, 
    line = 'co21', 
    vsys = None, 
    vwidth = None, 
    chan_fine = 0.0, 
    rebin_factor = 5, 
    do_regrid_first = True, 
    do_regrid_only = False, 
    do_statwt = False, 
    edge_for_statwt = -1, 
    overwrite = False, 
    ):
    """Extract a spectral line uv data from a measurement set with optimized regridding. 
    
    Extract a spectral line from a measurement set and regrid onto a
    new velocity grid with the desired spacing. This doesn't
    necessarily need the PHANGS keys in place and may be a general
    purpose utility. There are some minor subtleties here related to
    regridding and rebinning.
    
    This functions calls CASA mstransform three times:
        1. `mstransform(regridms=True, width=common_chanwidth*one_plus_eps)`
           Regrid the velocity axis to the common coarest channel width
           
        2. `mstransform(regridms=False, chanbin=rebin_factor)`
           Bin the velocity axis from common coarest channel width to the target channel width 
           by an integer rebin_factor (this can be 1 if target channel width is not large enough)
           
        3. `mstransform(regridms=False, chanaverage=False, combinespws=True)`
           Combine spws
    
    Alternatively, we can also do this with a different order (`regrid_method = 2`):
        1. `mstransform(regridms=False, chanbin=rebin_factor)`
           Rebin the velocity axis by an integer factor so as to get as close to 
           the target channel width as possible.
           
        2. `mstransform(regridms=True, width=target_chanwidth)`
           Regrid the velocity axis to the exact target channel width
        
        3. `mstransform(regridms=False, chanaverage=False, combinespws=True)`
           Combine spws
    
    Args:
        infile (str): The input measurement set data with suffix ".ms".
        outfile (str): The output measurement set data with suffix ".ms".
        line (str): Line name. 
        chan_fine (float): Channel width in units of km/s to regrid to. 
                           If it is -1 then we will keep the original channel width.
                           If do_regrid_only is True then this will be the exact output channel width.
        rebin_factor (int): channel rebin factor, must be an integer.  
        do_regrid_first (bool): If True, then first do regridding then do rebinning, and the final channel 
                                width will be `chan_fine * rebin_factor`. If False, then first do rebinning
                                then do regridding, and the final channel width will be `chan_fine`. 
                                The default is True. 
        do_regrid_only (bool): If True then only regridding will be done. The default is False, that is, 
                               we do both regridding and rebinning. 
        do_statwt (bool): 
        edge_for_statwt (int): 
    
    Inputs:
        infile: ALMA measurement set data folder.
    
    Outputs:
        outfile: ALMA measurement set data folder.
    
    """
    
    # 
    # This funcion is modified from the extract_line() function in older version phangsPipeline.py.
    # 
    
    # 
    # check input ms data dir
    if not os.path.isdir(infile):
        logger.error('Error! The input uv data measurement set "'+infile+'"does not exist!')
        raise Exception('Error! The input uv data measurement set "'+infile+'"does not exist!')
    # 
    # check output suffix
    if re.match(r'^(.*)\.ms$', outfile, re.IGNORECASE):
        out_name = re.sub(r'^(.*)\.ms$', r'\1', outfile, re.IGNORECASE)
        outfile = out_name + '.ms'
    else:
        out_name = outfile
        outfile = out_name + '.ms'
    # 
    # check existing output data
    if os.path.isdir(outfile) and not os.path.isdir(outfile+'.touch'):
        if not overwrite:
            logger.warning('Found existing output data "'+outfile+'", will not overwrite it.')
            return
    # if overwrite, then delete existing output data.
    for suffix in ['', '.flagversions', '.temp', '.temp.flagversions', '.temp2', '.temp2.flagversions', '.touch', '.temp.touch', '.temp2.touch']:
        if os.path.isdir(outfile+suffix):
            shutil.rmtree(outfile+suffix)
    # 
    # check vsys
    if vsys is None:
        logger.error('Error! Please input vsys for the galaxy systematic velocity in units of km/s.')
        raise Exception('Error! Please input vsys.')
    # 
    # check vwidth
    if vwidth is None:
        logger.error('Error! Please input vwidth for the galaxy systematic velocity in units of km/s.')
        raise Exception('Error! Please input vwidth.')
    # 
    # find spws for line
    spw_list_string = find_spws_for_line(infile, line, vsys = vsys, vwidth = vwidth)
    # 
    # exit if no line found
    if spw_list_string == '':
        # there has already a warning message inside find_spws_for_line()
        return
    # 
    # get the line name and line center rest-frame frequency in the line_list module for the input line
    line_name, restfreq_ghz = line_list.get_line_name_and_frequency(line, exit_on_error=True) # exit_on_error = True
    # 
    # print starting message
    logger.info("EXTRACT_LINE begins:")
    logger.info("... line: "+line)
    logger.info("... spectral windows to consider: "+spw_list_string)
    # 
    start_vel_kms = (vsys - vwidth/2.0)
    chan_width_hz = au.getChanWidths(infile, spw_list_string)
    if not np.isscalar(chan_width_hz): chan_width_hz = chan_width_hz[0]
    current_chan_width_kms = abs(chan_width_hz / (restfreq_ghz*1e9)*sol_kms)
    restfreq_string = "{:12.8f}".format(restfreq_ghz)+'GHz'
    start_vel_string =  "{:12.8f}".format(start_vel_kms)+'km/s'
    current_chanwidth_string = "{:12.8f}".format(current_chan_width_kms)+'km/s'
    if chan_fine > 0:
        nchan_for_recenter = int(np.max(np.ceil(vwidth / chan_fine)))
        chanwidth_string =  "{:12.8f}".format(chan_fine)+'km/s'
    else:
        nchan_for_recenter = int(np.max(np.ceil(vwidth / current_chan_width_kms)))
        chanwidth_string =  "{:12.8f}".format(current_chan_width_kms)+'km/s'
    # 
    logger.info("... rest frequency: "+restfreq_string)
    logger.info("... new starting velocity: "+start_vel_string)
    logger.info("... original velocity width: "+str(current_chan_width_kms))
    logger.info("... target velocity width: "+str(chan_fine))
    logger.info("... number of channels at this stage: "+str(nchan_for_recenter))
    # 
    # determine regridding/rebinning/combinespw order
    # we can do regridding-rebinning-combinespw, 
    # or rebinning-regridding-combinespw, 
    # or only regridding-combinespw. 
    mstransform_call_list = []
    mstransform_call_message_list = []
    # 
    # build default mstransform params for regridding
    mstransform_params_for_regrid = {'combinespws': False, 'regridms': True, 'mode': 'velocity', 'interpolation': 'cubic', 
                                     'spw': spw_list_string, 'datacolumn': 'DATA', 
                                     'start': start_vel_string.strip(), 'restfreq': restfreq_string.strip(), 
                                     'outframe': 'lsrk', 'veltype': 'radio', 
                                     'nchan': nchan_for_recenter, 'width': chanwidth_string.strip() }
    mstransform_message_for_regrid = 'mstransform to regrid to line velocity and channel width of '+chanwidth_string.strip()+' from the original channel width of '+current_chanwidth_string.strip()+'.'
    # 
    if chan_fine <= 0.0:
        del mstransform_params_for_regrid['width'] # if user has not input a valid chan_fine, then we will keep the original channel width.
        mstransform_message_for_regrid = 'mstransform to regrid to line velocity while keeping the original channel width of '+current_chanwidth_string.strip()+'.'
    # 
    # build default mstransform params for rebinning
    mstransform_params_for_rebin = {'combinespws': False, 'regridms': False, 
                                    'spw': '', 'datacolumn': 'DATA', 
                                    'chanaverage': True, 'chanbin': rebin_factor }
    mstransform_message_for_rebin = 'mstransform to rebin by a factor of '+str(rebin_factor)+'.'
    # 
    # build default mstransform params for combinespw
    mstransform_params_for_combinespw = {'combinespws': True, 'regridms': False, 
                                         'spw': '', 'datacolumn': 'DATA', 
                                         'chanaverage': False,
                                         'keepflags': False }
    mstransform_message_for_combinespw = 'mstransform to combine spws.'
    # 
    # first mstransform call
    if do_regrid_only or do_regrid_first:
        # do regrid
        mstransform_params = copy.copy(mstransform_params_for_regrid)
        mstransform_message = mstransform_message_for_regrid
        mstransform_call_list.append(mstransform_params)
        mstransform_call_message_list.append(mstransform_message)
    # 
    # second mstransform call (or the first if not do_regrid_only and not do_regrid_first)
    if not do_regrid_only and rebin_factor > 1:
        # do rebin if user has input a valid rebin_factor > 1
        # rebin will be done to all spws
        mstransform_params = copy.copy(mstransform_params_for_rebin)
        mstransform_message = mstransform_message_for_rebin
        mstransform_call_list.append(mstransform_params)
        mstransform_call_message_list.append(mstransform_message)
    # 
    # third mstransform call (or the second if not do_regrid_only and not do_regrid_first) (same as the first mstransform call)
    if not do_regrid_only and not do_regrid_first:
        # do regrid after rebin (if there is no rebin earlier because rebin_factor<=1, this will be the only mstransform call except for the last combinespw call.)
        mstransform_params = copy.copy(mstransform_params_for_regrid)
        mstransform_message = mstransform_message_for_regrid
        mstransform_call_list.append(mstransform_params)
        mstransform_call_message_list.append(mstransform_message)
    # 
    # last mstransform call, combine spw
    mstransform_params = copy.copy(mstransform_params_for_combinespw)
    mstransform_message = mstransform_message_for_combinespw
    mstransform_call_list.append(mstransform_params)
    mstransform_call_message_list.append(mstransform_message)
    # 
    # loop mstransform call list
    logger.info('... we will have '+str(len(mstransform_call_list))+' mstransform calls')
    for k, mstransform_params in enumerate(mstransform_call_list):
        if k == 0:
            mstransform_params['vis'] = infile
            mstransform_params['outputvis'] = outfile+'.temp%d'%(k+1)
        elif k == len(mstransform_call_list)-1:
            mstransform_params['vis'] = outfile+'.temp%d'%(k)
            mstransform_params['outputvis'] = outfile
        else:
            mstransform_params['vis'] = outfile+'.temp%d'%(k)
            mstransform_params['outputvis'] = outfile+'.temp%d'%(k+1)
        # 
        logger.info("... "+mstransform_call_message_list[k])
        logger.debug("... "+'mstransform('+', '.join("{!s}={!r}".format(t, mstransform_params[t]) for t in mstransform_params.keys())+')')
        os.mkdir(mstransform_params['outputvis']+'.touch')
        casaStuff.mstransform(**mstransform_params)
        os.rmdir(mstransform_params['outputvis']+'.touch')
    # 
    # Clean up
    if os.path.isdir(outfile):
        logger.info("... deleting temporary files")
        for k in range(len(mstransform_call_list)):
            for suffix in ['.temp%d'%(k), '.temp%d.flagversions'%(k), '.temp%d.touch'%(k)]:
                if os.path.isdir(outfile+suffix):
                    shutil.rmtree(outfile+suffix)
    
    # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
    # STEP 4. Re-weight the data using statwt
    # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
    
    # N.B. need the sliding time bin to make statwt work.
    
    if do_statwt:
        # 
        if edge_for_statwt == -1:
            exclude_str = ''
        else:
            nchan_final = int(np.floor(nchan_for_recenter / rebin_factor)+1)
            exclude_str = '*:'+str(edge_for_statwt-1)+'~'+\
                               str(nchan_final-(edge_for_statwt-2))
            logger.info("... running statwt with exclusion: "+exclude_str)
        # 
        if 'fitspw' in inspect.getargspec(casaStuff.statwt)[0]:
            # CASA version somewhat >= 5.5.0
            statwt_params = {'vis': outfile, 'timebin': '0.001s', 'slidetimebin': False, 'chanbin': 'spw', 
                             'statalg': 'classic', 'datacolumn': 'data', 
                             'fitspw': exclude_str, 'excludechans': True}
        else:
            # CASA version <= 5.4.1
            statwt_params = {'vis': outfile, 'timebin': '0.001s', 'slidetimebin': False, 'chanbin': 'spw', 
                             'statalg': 'classic', 'datacolumn': 'data', 
                             'excludechans': exclude_str}
        # 
        os.mkdir(outfile+'.touch')
        test = casaStuff.statwt(**statwt_params)
        os.rmdir(outfile+'.touch')
    
    logger.info("---")
    
    return


def extract_continuum(
    infile, 
    outfile, 
    lines_to_flag = None, 
    vsys = None, 
    vwidth = None, 
    do_statwt = False, 
    do_collapse = True, 
    overwrite = False, 
    ):
    """Extract continuum uv data from a measurement set and collapse into a single channel. 
    
    Extract a continuum measurement set, flagging any specified lines,
    reweighting using statwt, and then collapsing to a single "channel
    0" measurement.
    
    Args:
        infile (str): The input measurement set data with suffix ".ms".
        outfile (str): The output measurement set data with suffix ".ms".
        lines_to_flag (list): A list of line names to flag. Lines names must be in our line_list module. If it is None, then we use all 12co, 13co and c18o lines.
        do_statwt (bool): 
        do_collapse (bool): Always True to produce the single-channel continuum data.
    
    Inputs:
        infile: ALMA measurement set data folder.
    
    Outputs:
        outfile: ALMA measurement set data folder.
    
    """
    
    # 
    # check input ms data dir
    if not os.path.isdir(infile):
        logger.error('Error! The input uv data measurement set "'+infile+'"does not exist!')
        raise Exception('Error! The input uv data measurement set "'+infile+'"does not exist!')
    # 
    # check output suffix
    if re.match(r'^(.*)\.ms$', outfile, re.IGNORECASE):
        out_name = re.sub(r'^(.*)\.ms$', r'\1', outfile, re.IGNORECASE)
        outfile = out_name + '.ms'
    else:
        out_name = outfile
        outfile = out_name + '.ms'
    # 
    # check existing output data
    if os.path.isdir(outfile) and not os.path.isdir(outfile+'.touch'):
        if not overwrite:
            logger.warning('Found existing output data "'+outfile+'", will not overwrite it.')
            return
    # if overwrite, then delete existing output data.
    for suffix in ['', '.flagversions', '.temp', '.temp.flagversions', '.temp_copy', '.temp_copy.flagversions', '.touch']:
        if os.path.isdir(outfile+suffix):
            shutil.rmtree(outfile+suffix)
    # 
    # find_spw_channels_for_lines_to_flag
    spw_flagging_string = find_spw_channels_for_lines_to_flag(infile = infile, 
                                                              lines_to_flag = lines_to_flag, 
                                                              vsys = vsys, 
                                                              vwidth = vwidth)
    # 
    # make a continuum copy of the data
    os.mkdir(outfile+'.touch')
    shutil.copytree(infile, outfile)
    os.rmdir(outfile+'.touch')
    # 
    # flagdata
    if spw_flagging_string != '':
        os.mkdir(outfile+'.touch')
        casaStuff.flagdata(vis=outfile,
                           spw=spw_flagging_string,
                           )
        os.rmdir(outfile+'.touch')
    # 
    # statwt
    # Here - this comman needs to be examined and refined in CASA
    # 5.6.1 to see if it can be sped up. Right now things are
    # devastatingly slow.
    if do_statwt:
        logger.info("... deriving empirical weights using STATWT.")
        os.mkdir(outfile+'.touch')
        casaStuff.statwt(vis=outfile,
                         timebin='0.001s',
                         slidetimebin=False,
                         chanbin='spw',
                         statalg='classic',
                         datacolumn='data',
                         )
        os.rmdir(outfile+'.touch')
    # 
    # collapse
    if do_collapse:
        logger.info("... Collapsing the continuum to a single channel.")
        
        if os.path.isdir(outfile):
            shutil.move(outfile, outfile+'.temp_copy')
        if os.path.isdir(outfile+'.flagversions'):
            shutil.move(outfile+'.flagversions', outfile+'.temp_copy'+'.flagversions')
        
        os.mkdir(outfile+'.touch')
        casaStuff.split(vis=outfile+'.temp_copy',
                        outputvis=outfile,
                        width=10000,
                        datacolumn='DATA',
                        keepflags=False)
                        #<TODO><20200210># num_chan or width
        os.rmdir(outfile+'.touch')
        # 
        # clean up
        if os.path.isdir(outfile+'.temp_copy'):
            shutil.rmtree(outfile+'.temp_copy')
        if os.path.isdir(outfile+'.temp_copy.flagversions'):
            shutil.rmtree(outfile+'.temp_copy.flagversions')
    # 
    return

##################
# Analysis tasks #
##################

def noise_spectrum(
    vis=None,
    stat_name="medabsdevmed",
    start_chan=None,
    stop_chan=None):
    """
    Calculates the u-v based noise spectrum and returns it as an array. 
    """
    
    # This function is not used for now.
    
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
        logger.debug("Channel "+str(ii)+" / "+str(nchan))
        result = visstat(vis=vis,
                         axis='amp',
                         spw='0:'+str(ii),
                         )
        if result == None:
            logger.debug("Skipping channel.")
            continue
        spec[ii] = result[result.keys()[0]][stat_name]
        
    return spec








