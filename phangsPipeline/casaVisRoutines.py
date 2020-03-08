"""
Standalone routines to analyze and manipulate visibilities.
"""

# 20200226: introduced os.mkdir(outfile+'.touch') os.rmdir(outfile+'.touch') to make sure we can handle sudden system break.

#region Imports and definitions

import os, sys, re, shutil, inspect, copy
import numpy as np
from scipy.ndimage import label
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
import utilsLines as lines

# Pipeline versionining
from pipelineVersion import version as pipeVer

#endregion

#region Routines for basic characterization

#endregion

#region Routines to analyze and extract lines in measurement sets

# Physical constants
sol_kms = 2.9979246e5

##########################################
# Split, copy, combine measurement sets. #
##########################################
    
def copy_ms(
    infile = None, 
    outfile = None, 
    use_symlink = True, 
    overwrite = False, 
    ):
    """
    Copy a measurement set, optionally using symlink instead of
    actually copying.
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

    return()

def split_science_targets(
    infile = None, 
    outfile = None, 
    field = '', 
    intent = 'OBSERVE_TARGET*', 
    spw = '',
    timebin = '0s',
    do_statwt = False, 
    overwrite = False, 
    ):
    """
    Split science targets from the input ALMA measurement set to form
    a new, science-only measurement set. Optionally reweight the data
    using statwt. 

    Relatively thin wrapper to split that smooths out some things like
    handling of flagversions and which data tables to use.
    
    Args:

    infile (str): The input measurement set data.
    
    outfile (str): The output measurement set data.

    field, spw, intent (str): The field, spw, intent used for selection.
                
    timebin: The time bin applied.

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
        
    logger.debug('intent='+intent)
    logger.debug('field='+field)

    os.mkdir(outfile+'.touch')
    casaStuff.split(vis = infile, 
                    intent = intent, 
                    field = field,
                    spw = spw,
                    datacolumn = use_column, 
                    outputvis = outfile, 
                    keepflags = False,
                    timebin = timebin)
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
    thin wrapper to concat. Thin wrapper to concat. Might build out in
    the future.
    
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
    ranges_to_exclude = [],
    solint = 'int',
    fitorder = 0,
    combine = '',
    overwrite = False, 
    ):
    """
    Carry out uv continuum subtraction on a measurement set. First
    figures out channels corresponding to spectral lines for a
    provided suite of bright lines.
    """
        
    # Error and file existence checking

    if infile is None:
        logging.error("Please specify infile.")
        raise Exception("Please specify infile.")
    
    if outfile is None:
        outfile = infile+'.contsub'

    if not os.path.isdir(infile):
        logger.error('The input uv data measurement set "'+infile+'"does not exist.')
        return()

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

    # find_spw_channels_for_lines

    spw_flagging_string = spw_string_for_freq_ranges(
        infile = infile, 
        freq_ranges_ghz = ranges_to_exclude,
        fail_on_empty = True,
        )

    if spw_flagging_string is None:
        logger.error("All data are masked in at least one spectral window. Returning.")
        return()

    os.mkdir(infile+'.contsub'+'.touch')

    # uvcontsub, this outputs infile+'.contsub'

    casaStuff.uvcontsub(
        vis = infile,
        fitspw = spw_flagging_string,
        excludechans = True,
        combine=combine,
        fitorder=fitorder,
        solint=solint,
        want_cont=False
        )

    os.rmdir(infile+'.contsub'+'.touch')

    # Could manipulate outfile names here.

    return()

##########################################################
# Interface between spectral lines and spectral windows. #
##########################################################

def find_spws_for_line(
    infile = None, 
    line = None, restfreq_ghz = None,
    vsys_kms=None, vwidth_kms=None, vlow_kms=None, vhigh_kms=None,
    max_chanwidth_kms = None,
    exit_on_error = True, 
    as_list = False,
    ):
    """
    List the spectral windows in the input ms data that contains the
    input line, given the line velocity (vsys_kms) and line width
    (vwidth_kms), which are in units of km/s. Defaults to rest frequency
    (with vsys_kms = 0.0 and vwidth_kms = 0.0).
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

    if restfreq_ghz is None:
        if line is None:
            logging.error("Specify a line name or provide a rest frequency in GHz.")
            raise Exception("No rest frequency specified.")

        restfreq_ghz = (lines.get_line_name_and_frequency(line, exit_on_error=True))[1]

    # Work out the frequencies at the line edes.

    line_low_ghz, line_high_ghz = \
        lines.get_ghz_range_for_line(restfreq_ghz=restfreq_ghz, 
                                     vsys_kms=vsys_kms, vwidth_kms=vwidth_kms
                                     , vlow_kms=vlow_kms, vhigh_kms=vhigh_kms)

    # If channel width restrictions are in place, calculate the
    # implied channel width requirement in GHz.

    if max_chanwidth_kms is not None:
        line_freq_ghz = (line_low_ghz+line_high_ghz)*0.5
        max_chanwidth_ghz = line_freq_ghz*max_chanwidth_kms/sol_kms
    else:
        max_chanwidth_ghz = None

    # Work out which spectral windows contain the line by looping over
    # SPWs one at a time.

    spw_list = []
    
    vm = au.ValueMapping(infile)

    for this_spw in vm.spwInfo.keys():
        
        spw_high_ghz = np.max(vm.spwInfo[this_spw]['edgeChannels'])/1e9
        spw_low_ghz = np.min(vm.spwInfo[this_spw]['edgeChannels'])/1e9

        if spw_high_ghz < line_low_ghz:
            continue

        if spw_low_ghz > line_high_ghz:
            continue
        
        if max_chanwidth_ghz is not None:
            spw_chanwidth_ghz = abs(vm.spwInfo[this_spw]['chanWidth'])/1e9
            if spw_chanwidth_ghz > max_chanwidth_ghz:
                continue

        spw_list.append(this_spw)

    # If we don't find the line in this data set, issue a warning and
    # return.

    if len(spw_list) == 0:
        logger.warning('No spectral windows contain the input line.')
        spw_list_string = None # can't be '', that selects all
    else:
        
        # sort and remove duplicates
        spw_list = sorted(list(set(spw_list)))

        # make spw_list_string appropriate for use in selection
        spw_list_string = ','.join(np.array(spw_list).astype(str))
     
    if as_list:
        return(spw_list)
    else:
        return(spw_list_string)

def spw_string_for_freq_ranges(
    infile = None, 
    freq_ranges_ghz = [],
    just_spw = [],
    complement = False,
    fail_on_empty = False,
    ):

    """
    Given an input measurement set, return the spectral 
    List the spectral window and channels corresponding to the input
    lines in the input ms data.  Galaxy system velocity (vsys) and
    velocity width (vwidth) in units of km/s are needed.
    """
    
    # Check file existence

    if infile is None:
        logging.error("Please specify an input file.")
        raise Exception("Please specify an input file.")

    if not os.path.isdir(infile):
        logger.error('The input measurement set "'+infile+'"does not exist.')
        raise Exception('The input measurement set "'+infile+'"does not exist.')

    # Make sure that we have a list 
    if type(freq_ranges_ghz) != type([]):
        freq_ranges_ghz = [freq_ranges_ghz]

    if type(just_spw) != type([]):
        just_spw = [just_spw]

    vm = au.ValueMapping(infile)

    # Loop over spectral windows
    spw_flagging_string = ''
    first_string = True
    for this_spw in vm.spwInfo:
        
        if len(just_spw) > 0:
            if this_spw not in just_spw:
                continue

        freq_axis = vm.spwInfo[this_spw]['chanFreqs']
        half_chan = abs(freq_axis[1]-freq_axis[0])
        chan_axis = np.arange(len(freq_axis))
        mask_axis = np.zeros_like(chan_axis,dtype='bool')

        for this_freq_range in freq_ranges_ghz:
            
            low_freq_hz = this_freq_range[0]*1e9
            high_freq_hz = this_freq_range[1]*1e9

            ind = ((freq_axis-half_chan) >= low_freq_hz)*((freq_axis+half_chan) <= high_freq_hz)
            mask_axis[ind] = True
            
        if complement:
            mask_axis = np.invert(mask_axis) 

        if fail_on_empty:
            if np.sum(np.invert(mask_axis)) == 0:
                return(None)
        
        regions = (label(mask_axis))[0]
        max_reg = np.max(regions)
        for ii in range(1,max_reg+1):
            this_mask = (regions == ii)
            low_chan = np.min(chan_axis[this_mask])
            high_chan = np.max(chan_axis[this_mask])
            this_spw_string = str(this_spw)+':'+str(low_chan)+'~'+str(high_chan)
            if first_string:
                spw_flagging_string += this_spw_string
                first_string = False
            else:
                spw_flagging_string += ','+this_spw_string
        
    logger.info("... returning SPW selection string:")
    logger.info(spw_flagging_string)
    
    return(spw_flagging_string)

def compute_common_chanwidth(
    infile_list = None, 
    line = None, 
    vsys_kms = None, 
    vwidth_kms = None, 
    vlow_kms = None, 
    vhigh_kms = None, 
    ): 
    """
    Calculates the coarsest channel width among all spectral windows
    in the input measurement set that contain the input line.
    
    Args:
    
    Returns:
    
    """
        
    if infile_list is None:
        logging.error("Please specify one or more input files via infile_list.")
        Exception("Please specify one or more input files via infile_list.")

    if np.isscalar(infile_list):
        infile_list = [infile_list]

    # Get the line name and line center rest-frame frequency in the line_list module for the input line
    line_name, restfreq_ghz = lines.get_line_name_and_frequency(line, exit_on_error=True)

    # Work out the frequencies at the line edes and central frequency
    line_low_ghz, line_high_ghz = lines.get_ghz_range_for_line(line=line_name, vsys_kms=vsys_kms, vwidth_kms=vwidth_kms,
                                                               vlow_kms=vlow_kms, vhigh_kms=vhigh_kms)
    
    line_freq_ghz = (line_high_ghz+line_low_ghz)/2.0

    coarsest_channel = None
    for this_infile in infile_list:
        # Find spws for line
        spw_list_string = find_spws_for_line(this_infile, line, vsys_kms = vsys_kms, vwidth_kms = vwidth_kms)
        
        chan_widths_hz = au.getChanWidths(this_infile, spw_list_string)

        # Convert to km/s and return
        for this_chan_width_hz in chan_widths_hz:
            chan_width_kms = abs(this_chan_width_hz / (line_freq_ghz*1e9)*sol_kms)
            if coarsest_channel is None:
                coarsest_channel = chan_width_kms
            else:
                if chan_width_kms > coarsest_channel:
                    coarsest_channel = chan_width_kms

    return(coarsest_channel)

#################################################################
# Extract single-line or continuum data from a measurement set. #
#################################################################
    
def reweight_data():
    
    if edge_for_statwt == -1:
        exclude_str = ''
    else:
        nchan_final = int(np.floor(nchan_for_recenter / rebin_factor)+1)
        exclude_str = '*:'+str(edge_for_statwt-1)+'~'+\
            str(nchan_final-(edge_for_statwt-2))
        logger.info("... running statwt with exclusion: "+exclude_str)

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

def batch_extract_line(
    infile_list = [],
    
    ):
    

def extract_line(
    infile = None, 
    outfile = None, 
    line = 'co21', 
    restfreq_ghz = None,
    method = 'regrid_then_rebin',
    target_chanw_kms = None,
    vlow_kms = None,
    vhigh_kms = None,
    vsys_kms = None,
    vwidth_kms = None, 
    nchan = None,
    binfactor = None,
    overwrite = False, 
    ):

    # Check the method
    
    valid_methods = ['regrid_then_rebin','rebin_then_regrid','just_regrid','just_rebin']
    if method.lower().strip() not in valid_methods:
        logger.error("Not a valid line extraction medod - "+str(method))
        raise Exception("Please specify a valid line extraction method.")

    # Check input

    if infile is None:
        logging.error("Please specify an input file.")
        raise Exception("Please specify an input file.")

    if outfile is None:
        logging.error("Please specify an output file.")
        raise Exception("Please specify an output file.")

    if not os.path.isdir(infile):
        logger.error('The input measurement set "'+infile+'"does not exist.')
        raise Exception('The input measurement set "'+infile+'"does not exist.')

    # Check existence of output data and abort if found and overwrite is off

    if os.path.isdir(outfile) and  not os.path.isdir(outfile+'.touch'):            
        if not overwrite:
            logger.warning('Found existing output data "'+outfile+'", will not overwrite it.')
            return()

    # Else, clear all previous files and temporary files

    # TBD suffixes/syntax need updating
    for suffix in ['', '.flagversions', '.temp', '.temp.flagversions', 
                   '.temp2', '.temp2.flagversions', '.touch', '.temp.touch', '.temp2.touch']:
        if os.path.isdir(outfile+suffix):
            shutil.rmtree(outfile+suffix)

    # Get the line name and rest-frame frequency in the line_list
    # module for the input line.

    if restfreq_ghz is None:
        if line is None:
            logging.error("Specify a line name or provide a rest frequency in GHz.")
            raise Exception("No rest frequency specified.")

        restfreq_ghz = (lines.get_line_name_and_frequency(line, exit_on_error=True))[1]

    # Handle velocity windows, etc.
            
    # TBD

    # Identify SPWs - note whether we have multiple windows

    spw_list = find_spws_for_line(
        infile = line, restfreq_ghz = restfreq_ghz,
        vsys_kms=vsys_kms, vwidth_kms=vwidth_kms, vlow_kms=vlow_kms, vhigh_kms=vhigh_kms,
        exit_on_error = True, 
        as_list = True,
        )
    if spw_list is None:
        logging.error("No SPWs for selected line and velocity range.")
        return()
    multiple_spws = len(spw_list) > 1
    spw = spw_list.join(',')

    # ............................................
    # Initialize the calls
    # ............................................
        
    if method == 'just_regrid' or method == 'regrid_then_rebin' or \
            method == 'rebin_then_regrid':

        regrid_params, regrid_msg =  build_mstransform_call(
            infile=infile, outfile=outfile, restfreq_ghz=restfreq_ghz, spw=spw,
            vstart_kms=vstart_kms, nchan=nchan, method='regrid')
                               
    if method == 'just_rebin' or method == 'regrid_then_rebin' or \
            method == 'rebin_then_regrid':

        # Verify that we have the bin factor

        rebin_params, rebin_msg =  build_mstransform_call(
            infile=infile, outfile=outfile, restfreq_ghz=restfreq_ghz, spw=spw,
            binfactor=binfactor, method='rebin')

    if multiple_spws:
        combine_params, combine_msg =  build_mstransform_call(
            infile=infile, outfile=outfile, restfreq_ghz=restfreq_ghz, spw=spw,
            method='combine')

    # ............................................
    # string the calls together in the desired order
    # ............................................

    params_list = []
    msg_list = []
    if method == 'just_regrid':
        params_list.append(regrid_params)
        msg_list.append(regrid_msg)

    if method == 'just_rebin':
        params_list.append(rebin_params)
        msg_list.append(rebin_msg)

    if method == 'rebin_then_regrid':
        params_list.append(rebin_params)
        msg_list.append(rebin_msg)

        params_list.append(regrid_params)
        msg_list.append(regrid_msg)

    if method == 'regrid_then_rebin':
        params_list.append(regrid_params)
        msg_list.append(regrid_msg)

        params_list.append(rebin_params)
        msg_list.append(rebin_msg)

    if multiple_spws:
        params_list.append(combine_params)
        msg_lis.append(combine_msg)

    # ............................................
    # Execute the list of mstransform calls
    # ............................................

    n_calls = len(param_list)
    logger.info('... we will have '+str(n_calls)+' mstransform calls')

    for kk in range(n_calls):
        this_params = params_list[kk]
        this_msg = msg_list[kk]
        if kk == 0:
            this_params['vis'] = infile
            this_params['outputvis'] = outfile+'.temp%d'%(kk+1)
        elif kk == n_calls-1:
            this_params['vis'] = outfile+'.temp%d'%(kk)
            this_params['outputvis'] = outfile
        else:
            this_params['vis'] = outfile+'.temp%d'%(kk)
            this_params['outputvis'] = outfile+'.temp%d'%(kk+1)
         
        logger.info("... "+this_msf)
        logger.debug("... "+'mstransform('+', '.join("{!s}={!r}".format(t, this_params[t]) for t in this_params.keys())+')')

        os.mkdir(mstransform_params['outputvis']+'.touch')
        casaStuff.mstransform(**this_params)
        os.rmdir(mstransform_params['outputvis']+'.touch')

    # ............................................
    # Clean up leftover files
    # ............................................

    # TBD revisit

    if os.path.isdir(outfile):
        logger.info("... deleting temporary files")
        for k in range(len(mstransform_call_list)):
            for suffix in ['.temp%d'%(k), '.temp%d.flagversions'%(k), '.temp%d.touch'%(k)]:
                if os.path.isdir(outfile+suffix):
                    shutil.rmtree(outfile+suffix)
    
    logger.info("---")

def build_mstransform_call(
    infile = None, 
    outfile = None, 
    restfreq_ghz = None,
    spw = None,
    vstart_kms = None,
    vwidth_kms = None,
    datacolumn = None,
    method = 'regrid',
    target_chan_kms = None,
    nchan = None,
    binfactor = None,
    overwrite = False, 
    ):
    """    
    Extract a spectral line from a measurement set and regrid onto a
    new velocity grid with the desired spacing. There are some minor
    subtleties here related to regridding and rebinning.
    
    """

    # ............................................
    # Error checking and setup
    # ............................................

    # Check that the requested method is understood

    valid_methods = ['rebin','regrid','combine']
    if method.lower().strip() not in valid_methods:
        logger.error("Not a valid line extraction medod - "+str(method))
        raise Exception("Please specify a valid line extraction method.")
    
    # Check input

    if infile is None:
        logging.error("Please specify an input file.")
        raise Exception("Please specify an input file.")

    if outfile is None:
        logging.error("Please specify an output file.")
        raise Exception("Please specify an output file.")

    # If not supplied by the user, find which SPWs should be included
    # in the processing.
    
    if spw is None:
        if restfreq_ghz is not None:
            spw = find_spws_for_line(
                infile=infile, restfreq_ghz = restfreq_ghz,
                vlow_kms=vlow_kms, vhigh_kms=vhigh_kms)

            # Exit if no SPWs contain the line.
            if spw is None:
                # there has already a warning message inside find_spws_for_line()
                return()
        else:
            logger.info("Defaulting to all SPW selections.")
            spw = ''

    # Determine the column to use

    if datacolumn is None:
        casaStuff.tb.open(infile, nomodify = True)
        colnames = casaStuff.tb.colnames()
        if 'CORRECTED_DATA' in colnames:
            logger.info("Data has a CORRECTED column. Will use that.")
            datacolumn = 'CORRECTED'
        else:
            logger.info("Data lacks a CORRECTED column. Will use DATA column.")
            datacolumn = 'DATA'
        casaStuff.tb.close()

    # ............................................
    # Common parameters
    # ............................................

    params = {'vis': infile, 'outputvis': outfile, 'datacolumn':datacolumn,
              'spw': spw, 
              }

    # ............................................
    # Regridding
    # ............................................
    
    if method == 'regrid':

        # Check that we are provided a rest frequency or line name
        
        if restfreq_ghz is None:
            logger.error("Please specify a rest frequency in GHz.")
            raise Exception("No rest frequency specified.")
        restfreq_string = ("{:12.8f}".format(restfreq_ghz)+'GHz').strip()

        # Check that we have a velocity start and width

        if vstart_kms is None:
            logger.error("Please specify a starting velocity in km/s.")
            raise Exception("No starting velocity specified.")        
        start_vel_string =  ("{:12.8f}".format(vlow_kms)+'km/s').strip()
        
        # Check that we have a velocity width

        if vwidth_kms is None and nchan is None:
            logger.error("Please specify a velocity width in km/s or number of channels.")
            raise Exception("No starting velocity specified.")

        # Figure out the channel spacing, catching a few possible errors
        
        max_chan_hz = np.max(np.abs(au.getChanWidths(infile, spw)))
        current_chan_kms = max_chan_hz/restfreq_ghz*sol_kms

        skip_width = False
        if target_chan_kms is None:
            target_chan_kms = current_chan_kms
            skip_width = True
        elif current_chan_kms > target_chan_kms:
            target_chan_kms = current_chan_kms
            skip_width = True

        chanwidth_string =  "{:12.8f}".format(target_chan_kms)+'km/s'

        # Figure the number of channels if not supplied

        if nchan is None:
            nchan = int(np.max(np.ceil(vwidth_kms / target_chan_kms)))
        
        params = {'combinespws': False, 'regridms': True, 'chanaverage': False,
                  'mode': 'velocity', 'interpolation': 'cubic', 
                  'outframe': 'lsrk', 'veltype': 'radio', 'restfreq': restfreq_string, 
                  'start': start_vel_string, 'nchan': nchan, 'width': chanwidth_string }

        if skip_width:
            del params['width']

        message = '... regrid channel width '+chanwidth_string+' and nchan '+str(nchan)

    # ............................................
    # Rebin
    # ............................................
    
    if method == 'rebin':
        
        params.update({'combinespws': False, 'regridms': False, 'chanaverage' : True, 
                       'chanbin': binfactor }

        message = '... rebin by a factor of '+str(binfactor)

    # ............................................
    # Combine SPWs
    # ............................................
    
    if method = 'combine':
     
        params.update({'combinespws': True, 'regridms': False, 'chanaverage' : False, 
                       'keepflags': False })
        
        message = '... combine attempting to merge spectral windows.'

    return(params, message)


def extract_continuum(
    infile, 
    outfile, 
    lines_to_flag = None, 
    vsys_kms = None, 
    vwidth_kms = None, 
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
    # find_spw_channels_for_lines
    spw_flagging_string = find_spw_channels_for_lines(infile = infile, 
                                                      lines_to_flag = lines_to_flag, 
                                                      vsys_kms = vsys_kms, 
                                                      vwidth_kms = vwidth_kms)
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








