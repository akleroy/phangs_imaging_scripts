"""
Standalone routines to analyze and manipulate visibilities.
"""

# 20200226: introduced os.mkdir(outfile+'.touch') os.rmdir(outfile+'.touch') to make sure we can handle sudden system break.

#region Imports and definitions

import os, sys, re, shutil, inspect, copy
import numpy as np
from scipy.ndimage import label
#import pyfits # CASA has pyfits, not astropy
import glob

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
sol_kms = 2.99792458e5

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

        if not os.path.isdir(outfile+'.touch'):
            os.mkdir(outfile+'.touch')

        if os.path.isdir(infile):
            shutil.copytree(infile, outfile)
            logger.debug('shutil.copytree "'+infile+'", "'+outfile+'"')

        if os.path.isdir(infile+'.flagversions'):
            shutil.copytree(infile+'.flagversions', outfile+'.flagversions')
            logger.debug('shutil.copytree "'+infile+'.flagversions'+'", "'+outfile+'.flagversions'+'"')

        if os.path.isdir(outfile+'.touch'):
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
        
    logger.info('... intent: '+intent)
    logger.info('... field: '+field)
    logger.info('... spw: '+spw)
    
    if not os.path.isdir(outfile+'.touch'):
        os.mkdir(outfile+'.touch') # mark the beginning of our processing
        
    split_params = {'vis':infile, 'intent':intent, 'field':field, 'spw':spw, 
                    'datacolumn':use_column, 'outputvis':outfile, 
                    'keepflags':False, 'timebin':timebin}
    
    logger.info("... running CASA "+'split('+', '.join("{!s}={!r}".format(k, split_params[k]) for k in split_params.keys())+')')
    
    casaStuff.split(**split_params)
    
    # Re-weight the data if desired.
    
    if do_statwt:
        logger.info("Using statwt to re-weight the data.")
        statwt_params = {'vis':outfile, 'datacolumn':'DATA'}
        logger.info("... running CASA "+'statwt('+', '.join("{!s}={!r}".format(k, statwt_params[k]) for k in statwt_params.keys())+')')
        casaStuff.statwt(vis = outfile, 
                         datacolumn = 'DATA')
    
    if os.path.isdir(outfile+'.touch'):
        os.rmdir(outfile+'.touch') # mark the end of our processing

    return()

def concat_ms(
    infile_list=None,  
    outfile=None,  
    freqtol='',
    dirtol='',
    copypointing=True,
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
    concat_params = {'vis':infile_list, 'concatvis':outfile, 'copypointing':copypointing}
    if freqtol is not None and freqtol != '': concat_params['freqtol'] = freqtol
    if dirtol is not None and dirtol != '': concat_params['dirtol'] = dirtol
    logger.info("... running CASA "+'concat('+', '.join("{!s}={!r}".format(k, concat_params[k]) for k in concat_params.keys())+')')
    
    if not os.path.isdir(outfile+'.touch'):
        os.mkdir(outfile+'.touch') # mark the beginning of our processing
    
    casaStuff.concat(**concat_params)
    
    if os.path.isdir(outfile+'.touch'):
        os.rmdir(outfile+'.touch') # mark the end of our processing

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

    # uvcontsub, this outputs infile+'.contsub'
    
    uvcontsub_params = {'vis':infile, 'fitspw':spw_flagging_string, 'excludechans':True, 'combine':combine, 
                        'fitorder':fitorder, 'solint':solint, 'want_cont':False}
    logger.info("... running CASA "+'uvcontsub('+', '.join("{!s}={!r}".format(k, uvcontsub_params[k]) for k in uvcontsub_params.keys())+')')

    if not os.path.isdir(infile+'.contsub'+'.touch'):
        os.mkdir(infile+'.contsub'+'.touch') # mark the beginning of our processing

    casaStuff.uvcontsub(**uvcontsub_params)

    if os.path.isdir(infile+'.contsub'+'.touch'):
        os.rmdir(infile+'.contsub'+'.touch') # mark the end of our processing

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
    require_data = False,
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
    logger.debug("... line: %s, line freq: %.6f - %.6f, rest-freq: %.6f"%(line, line_low_ghz, line_high_ghz, restfreq_ghz))

    # If channel width restrictions are in place, calculate the
    # implied channel width requirement in GHz.

    if max_chanwidth_kms is not None:
        line_freq_ghz = (line_low_ghz+line_high_ghz)*0.5

        # Using RADIO convention for velocities:

        max_chanwidth_ghz = restfreq_ghz*max_chanwidth_kms/sol_kms

        # using the high-z convention used to be this (delete eventually)
        #max_chanwidth_ghz = line_freq_ghz*max_chanwidth_kms/sol_kms

        logger.debug("... max_chanwidth_kms: %.3f, max_chanwidth_ghz: %.6f"%(max_chanwidth_kms, max_chanwidth_ghz))
    else:
        max_chanwidth_ghz = None

    # Work out which spectral windows contain the line by looping over
    # SPWs one at a time.

    spw_list = []
    
    logger.debug("... vm = au.ValueMapping(infile) ...")
    vm = au.ValueMapping(infile)
    logger.debug("... vm = au.ValueMapping(infile) done")

    for this_spw in vm.spwInfo.keys():
        
        spw_high_ghz = np.max(vm.spwInfo[this_spw]['edgeChannels'])/1e9
        spw_low_ghz = np.min(vm.spwInfo[this_spw]['edgeChannels'])/1e9
        logger.debug("... spw: %s, freq: %.6f - %.6f GHz"%(this_spw, spw_low_ghz, spw_high_ghz))

        if spw_high_ghz < line_low_ghz:
            continue

        if spw_low_ghz > line_high_ghz:
            continue
        
        if max_chanwidth_ghz is not None:
            spw_chanwidth_ghz = abs(vm.spwInfo[this_spw]['chanWidth'])/1e9
            if spw_chanwidth_ghz > max_chanwidth_ghz:
                continue

        if require_data:
            if len(vm.scansForSpw[spw]) == 0:
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

def find_spws_for_science(
    infile = None, 
    require_data = False,
    exit_on_error = True, 
    as_list = False,
    ):
    """
    List all spectral windows that we judge likely to be used for
    science. Mostly wraps analysisUtils rather than reinventing the
    wheel.
    """

    # Check inputs
    
    if infile is None:
        logging.error("Please specify infile.")
        raise Exception("Please specify infile.")

    # Verify file existence

    if not os.path.isdir(infile):
        logger.error('Error! The input uv data measurement set "'+infile+'"does not exist!')
        raise Exception('Error! The input uv data measurement set "'+infile+'"does not exist!')

    # Call the analysisUtil version.

    spw_string = au.getScienceSpws(infile, intent = 'OBSERVE_TARGET*')
    spw_list = []
    for this_spw_string in spw_string.split(','):
        spw_list.append(int(this_spw_string))
    
    # Shouldn't get here, I think, because of the analysisUtils logic

    if require_data:
        vm = au.ValueMapping(infile)

        for spw in spw_list:
            if len(vm.scansForSpw[spw]) == 0:
                spw_list.remove(spw)
    # Return
    
    if len(spw_list) == 0:
        logger.warning('No science spectral windows found.')
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
        half_chan = abs(freq_axis[1]-freq_axis[0])*0.5
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
    line_low_ghz, line_high_ghz = lines.get_ghz_range_for_line(
        line=line_name, vsys_kms=vsys_kms, vwidth_kms=vwidth_kms,
        vlow_kms=vlow_kms, vhigh_kms=vhigh_kms)
    
    line_freq_ghz = (line_high_ghz+line_low_ghz)/2.0

    coarsest_channel = None
    for this_infile in infile_list:
        # Find spws for line
        spw_list_string = find_spws_for_line(this_infile, line, vsys_kms = vsys_kms, vwidth_kms = vwidth_kms)
        
        chan_widths_hz = au.getChanWidths(this_infile, spw_list_string)

        # Convert to km/s and return
        for this_chan_width_hz in chan_widths_hz:

            # Using RADIO convention for velocities:
            chan_width_kms = abs(this_chan_width_hz / (line_freq_ghz*1e9)*sol_kms)
            if coarsest_channel is None:
                coarsest_channel = chan_width_kms
            else:
                if chan_width_kms > coarsest_channel:
                    coarsest_channel = chan_width_kms

    return(coarsest_channel)

#######################################################
# Extract a single-line, common grid measurement set. #
#######################################################

def batch_extract_line(
    infile_list = [],
    outfile = None,
    target_chan_kms = None,
    restfreq_ghz = None,
    line = None,
    vsys_kms=None, vwidth_kms=None, vlow_kms=None, vhigh_kms=None,
    method = 'regrid_then_rebin',
    exact = False,
    freqtol = '',
    clear_pointing = True,
    overwrite = False,
    ):
    """
    Run a batch line extraction.
    """

    # Check that we have an output file defined.

    if outfile is None:
        logging.error("Please specify an output file.")
        raise Exception("Please specify an output file.")

    # Check existence of output data and abort if found and overwrite is off

    if os.path.isdir(outfile) and not os.path.isdir(outfile+'.touch'):            
        if not overwrite:
            logger.warning('... found existing output data "'+outfile+'", will not overwrite it.')
            return()

    # Else, clear all previous files and temporary files
    
    for suffix in ['', '.flagversions', '.touch', '.temp*']:
        if not (suffix.find('*') >= 0):
            if os.path.isdir(outfile+suffix):
                logger.debug('... shutil.rmtree(%r)'%(outfile+suffix))
                shutil.rmtree(outfile+suffix)
        else:
            for temp_outfile in glob.glob(outfile+suffix):
                if os.path.isdir(temp_outfile):
                    logger.debug('... shutil.rmtree(%r)'%(temp_outfile))
                    shutil.rmtree(temp_outfile)
    
    # Feed directly to generate an extraction scheme. This does a lot
    # of the error checking.

    schemes = suggest_extraction_scheme(
        infile_list = infile_list,
        target_chan_kms = target_chan_kms, method = method, exact = exact,
        restfreq_ghz = restfreq_ghz, line = line,
        vsys_kms=vsys_kms, vwidth_kms=vwidth_kms, vlow_kms=vlow_kms, vhigh_kms=vhigh_kms,        
        )

    # Execute the extraction scheme
    split_file_list = []
    for this_infile in schemes.keys():
        
        for this_spw in schemes[this_infile].keys():
            this_scheme = schemes[this_infile][this_spw]

            # Specify output file and check for existence
            this_outfile = this_infile+'.temp_spw'+str(this_spw).strip()

            this_scheme['outfile'] = this_outfile
            this_scheme['overwrite'] = overwrite
            split_file_list.append(this_outfile)

            # Execute line extraction
            
            del this_scheme['chan_width_kms']
            del this_scheme['chan_width_ghz']
            extract_line(**this_scheme)

            # Deal with pointing table - testing shows it to be a
            # duplicate for each SPW here, so we remove all rows for
            # all SPWs except the first one.

            if clear_pointing:
                # This didn't work:
                # os.system('rm -rf '+this_outfile+'/POINTING')

                # This zaps the whole table:
                au.clearPointingTable(this_outfile)

    # Concatenate and combine the output data sets

    if clear_pointing:
        copy_pointing = False
    else:
        copy_pointing = True
    
    concat_ms(
        infile_list = split_file_list,
        outfile = outfile,
        overwrite = overwrite,
        freqtol = freqtol, 
        copypointing = copy_pointing)

    # Clean up, deleting intermediate files

    for this_file in split_file_list:
        shutil.rmtree(this_file)
    
    return()

def choose_common_res(
    vals=[],
    epsilon=1e-4):
    """
    Choose a common resolution given a list and an inflation
    parameter epsilon. Returns max*(1+epsilon).).
    """
    if len(vals) == 0:
        return(None)
    ra = np.array(np.abs(vals))
    common_res = np.max(ra)*(1.+epsilon)
    return(common_res)

def suggest_extraction_scheme(
    infile_list = [],
    target_chan_kms = None,
    restfreq_ghz = None,
    line = None,
    vsys_kms=None, vwidth_kms=None, vlow_kms=None, vhigh_kms=None,
    method = 'regrid_then_rebin',
    exact = False,
    ):
    """
    Recommend extraction parameters given an input list of files, a
    desired target channel width, and a preferred algorithm. Returns a
    dictionary suitable for putting into the extraction routine.
    """

    # Check inputs
    
    if infile_list is None:
        logging.error("Please specify a list of infiles.")
        raise Exception("Please specify a list of infiles.")

    # make sure the input infile_list is a list
    
    if np.isscalar(infile_list):
        infile_list = [infile_list]    

    # Require a valid method choice

    valid_methods = ['regrid_then_rebin','rebin_then_regrid','just_regrid','just_rebin']
    if method.lower().strip() not in valid_methods:
        logger.error("Not a valid line extraction medod - "+str(method))
        raise Exception("Please specify a valid line extraction method.")

    # Get the line name and rest-frame frequency in the line_list
    # module for the input line

    if restfreq_ghz is None:
        if line is None:
            logging.error("Specify a line name or provide a rest frequency in GHz.")
            raise Exception("No rest frequency specified.")
        
        restfreq_ghz = (lines.get_line_name_and_frequency(line, exit_on_error=True))[1]

    # Work out the frequencies at the line edes.

    line_low_ghz, line_high_ghz = lines.get_ghz_range_for_line(
        restfreq_ghz=restfreq_ghz,  vsys_kms=vsys_kms, vwidth_kms=vwidth_kms,
        vlow_kms=vlow_kms, vhigh_kms=vhigh_kms)
    
    line_freq_ghz = 0.5*(line_low_ghz+line_high_ghz)

    # ----------------------------------------------------------------
    # Loop over infiles and spectral windows and record information
    # ----------------------------------------------------------------

    scheme = {}
    chan_width_list = []
    binfactor_list = []

    for this_infile in infile_list:

        vm = au.ValueMapping(this_infile)
        
        spw_list = vm.spwInfo.keys()
        
        scheme[this_infile] = {}

        for this_spw in spw_list:

            if len(vm.scansForSpw[this_spw]) == 0: 
                continue
            
            chan_width_ghz = np.abs(vm.spwInfo[this_spw]['chanWidth'])/1e9

            # Using RADIO convention:
            chan_width_kms = chan_width_ghz/restfreq_ghz * sol_kms
            # using the old relative convention was this (can delete eventually)
            #chan_width_kms = chan_width_ghz/line_freq_ghz * sol_kms 

            if chan_width_kms > target_chan_kms:
                logger.warning("Channel too big for SPW "+str(this_spw))
                continue
            
            chan_width_list.append(chan_width_kms)
            this_binfactor = int(np.floor(target_chan_kms/chan_width_kms))
            binfactor_list.append(this_binfactor)

            # Record basic file information
            scheme[this_infile][this_spw] = {}
            scheme[this_infile][this_spw]['infile'] = this_infile
            scheme[this_infile][this_spw]['spw'] = str(this_spw)

            # Record the line information
            scheme[this_infile][this_spw]['restfreq_ghz'] = restfreq_ghz
            scheme[this_infile][this_spw]['line'] = line
            scheme[this_infile][this_spw]['vlow_kms'] = vlow_kms
            scheme[this_infile][this_spw]['vhigh_kms'] = vhigh_kms
            scheme[this_infile][this_spw]['vsys_kms'] = vsys_kms
            scheme[this_infile][this_spw]['vwidth_kms'] = vwidth_kms

            # Record method information            
            scheme[this_infile][this_spw]['method'] = method
            scheme[this_infile][this_spw]['binfactor'] = this_binfactor
            scheme[this_infile][this_spw]['target_chan_kms'] = None

            # Record channel width information
            scheme[this_infile][this_spw]['chan_width_kms'] = chan_width_kms
            scheme[this_infile][this_spw]['chan_width_ghz'] = chan_width_ghz

    # ----------------------------------------------------------------
    # Figure out the strategy
    # ----------------------------------------------------------------
    
    # ... for rebinning, just do the naive division of floor(target / current)
    if method == 'just_rebin':
        return(scheme)

    # ... for regridding, just regrid to the desired width
    if method == 'just_regrid':
        for this_infile in scheme.keys():
            for this_spw in scheme[this_infile].keys():
                if exact:
                    scheme[this_infile][this_spw]['target_chan_kms'] = target_chan_kms
                else:
                    # Could revise this ... not positive of the correct choice
                    scheme[this_infile][this_spw]['target_chan_kms'] = target_chan_kms
    
    # ... for rebin-then-regrid, first rebin by the naive amount. Then
    # regrid either to the final value (if exact) or to a common
    # resolution determined by the actual channels and rebinnning.

    if method == 'rebin_then_regrid':
        for this_infile in scheme.keys():
            for this_spw in scheme[this_infile].keys():
                if exact:
                    scheme[this_infile][this_spw]['target_chan_kms'] = target_chan_kms
                else:
                    common_res = choose_common_res(
                        vals=np.array(chan_width_list)*np.array(binfactor_list),
                        epsilon=3e-4)
                    scheme[this_infile][this_spw]['target_chan_kms'] = common_res
                    #print("chan_width_list: ", chan_width_list)
                    #print("binfactor_list: ", binfactor_list)
                    #print("common res/target channel: ", common_res)

    # ... for regrid-then-rebin        

    if method == 'regrid_then_rebin':
        for this_infile in scheme.keys():
            for this_spw in scheme[this_infile].keys():
                if exact:
                    scheme[this_infile][this_spw]['target_chan_kms'] = \
                        target_chan_kms / scheme[this_infile][this_spw]['binfactor']
                else:
                    common_res = choose_common_res(
                        vals=np.array(chan_width_list)*np.array(binfactor_list),
                        epsilon=3e-4)
                    this_target_chan_kms = common_res / scheme[this_infile][this_spw]['binfactor']
                    scheme[this_infile][this_spw]['target_chan_kms'] = \
                        this_target_chan_kms
                    #print("chan_width_list: ", chan_width_list)
                    #print("binfactor_list: ", binfactor_list)
                    #print("common res: ", common_res)
                    #print("target channel: ", this_target_chan_kms)

    # Return

    return(scheme)
    
def extract_line(
    infile = None, 
    outfile = None, 
    spw = None,
    restfreq_ghz = None,
    line = 'co21', 
    vlow_kms = None,
    vhigh_kms = None,
    vsys_kms = None,
    vwidth_kms = None, 
    method = 'regrid_then_rebin',
    target_chan_kms = None,
    nchan = None,
    binfactor = None,
    overwrite = False, 
    ):
    """
    Line extraction routine. Takes infile, outfile, line of interest,
    and algorithm, along with algorithm tuning parameters.
    """

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

    if os.path.isdir(outfile) and not os.path.isdir(outfile+'.touch'):            
        if not overwrite:
            logger.warning('... found existing output data "'+outfile+'", will not overwrite it.')
            return()

    # Else, clear all previous files and temporary files

    # TBD suffixes/syntax need updating
    #<DL> updated this TBD
    for suffix in ['', '.flagversions', '.touch', '.temp*']:
        if not (suffix.find('*') >= 0):
            if os.path.isdir(outfile+suffix):
                logger.debug('... shutil.rmtree(%r)'%(outfile+suffix))
                shutil.rmtree(outfile+suffix)
        else:
            for temp_outfile in glob.glob(outfile+suffix):
                if os.path.isdir(temp_outfile):
                    logger.debug('... shutil.rmtree(%r)'%(temp_outfile))
                    shutil.rmtree(temp_outfile)
    
    # Create touch file to mark that we are processing this data
    if not os.path.isdir(outfile+'.touch'):
        os.mkdir(outfile+'.touch')
    
    # Get the line name and rest-frame frequency in the line_list
    # module for the input line.

    if restfreq_ghz is None:
        if line is None:
            logging.error("Specify a line name or provide a rest frequency in GHz.")
            raise Exception("No rest frequency specified.")

        restfreq_ghz = (lines.get_line_name_and_frequency(line, exit_on_error=True))[1]

    # Handle velocity windows, etc.
            
    if vsys_kms is None and vwidth_kms is None:
        if vlow_kms is not None and vhigh_kms is not None:
            vsys_kms = 0.5*(vhigh_kms+vlow_kms)
            vwidth_kms = (vhigh_kms - vlow_kms)
        else:
            if method == 'just_regrid' or method == 'regrid_then_rebin' or \
                    method == 'rebin_then_regrid':
                logging.error("I need a velocity width for a regridding step.")
                raise Exception("Need velocity width for regridding.")

    # ... should only reach this next block in the "just_rebinning" case

    if (vsys_kms is None) and (vwidth_kms is None) and \
            (vlow_kms is None) and (vhigh_kms is None):
        logging.error("Missing velocities. Setting to zero as a placeholder.")
        vsys_kms = 0.0
        vwidth_kms = 0.0

    # ... if now SPW selection string is provided then note whether we
    # have multiple windows.

    if spw is None:
        spw_list = find_spws_for_line(
            infile = line, restfreq_ghz = restfreq_ghz,
            vsys_kms=vsys_kms, vwidth_kms=vwidth_kms, vlow_kms=vlow_kms, vhigh_kms=vhigh_kms,
            require_data = True, exit_on_error = True, as_list = True,
            )
        if spw_list is None:
            logging.error("No SPWs for selected line and velocity range.")
            return()
        spw = spw_list.join(',')
    else:
        spw_list = spw.split(',')

    multiple_spws = len(spw_list) > 1

    # ............................................
    # Initialize the calls
    # ............................................
        
    if method == 'just_regrid' or method == 'regrid_then_rebin' or \
            method == 'rebin_then_regrid':

        if target_chan_kms is None:
            logger.warning('Need a target channel width to enable regridding.')
            return()

        vstart_kms = vsys_kms - vwidth_kms/2.0

        regrid_params, regrid_msg =  build_mstransform_call(
            infile=infile, outfile=outfile, restfreq_ghz=restfreq_ghz, spw=spw,
            vstart_kms=vstart_kms, vwidth_kms=vwidth_kms,
            target_chan_kms=target_chan_kms, nchan=nchan, 
            method='regrid')
                               
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

    n_calls = len(params_list)
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
         
        if os.path.isdir(this_params['outputvis']):
            shutil.rmtree(this_params['outputvis'])
        #if os.path.isdir(this_params['outputvis']+'.touch'):
        #    shutil.rmtree(this_params['outputvis']+'.touch')
        #    #if overwrite:
        #    #    shutil.rmtree(this_params['outputvis'])
        #    #else:
        #    #    logger.error("Intermediate file in place and overwrite=False. Returning.")
        #    #    return()
        #    #<DL> we already have overwrite check above, here we will just delete any existing old intermediate files

        logger.info("... "+this_msg)

        # in the case where we are in subsequent split, we expect a
        # single SPW and to use the data column. 

        if kk > 0:
            this_params['spw']=''
            this_params['datacolumn']='DATA'
        
        logger.info("... running CASA "+'mstransform('+', '.join("{!s}={!r}".format(t, this_params[t]) for t in this_params.keys())+')')

        if not os.path.isdir(this_params['outputvis']+'.touch'):
            os.mkdir(this_params['outputvis']+'.touch') # mark the beginning of our processing
        
        casaStuff.mstransform(**this_params)
        
        if os.path.isdir(this_params['outputvis']+'.touch') and this_params['outputvis'] != outfile:
            os.rmdir(this_params['outputvis']+'.touch') # mark the end of our processing

    # ............................................
    # Clean up leftover files
    # ............................................

    # TBD revisit

    if os.path.isdir(outfile):
        logger.info("... deleting temporary files")
        for kk in range(len(params_list)):
            for suffix in ['.temp%d'%(kk), '.temp%d.flagversions'%(kk), '.temp%d.touch'%(kk)]:
                if os.path.isdir(outfile+suffix):
                    shutil.rmtree(outfile+suffix)
    
    # Remove touch file to mark that we are have done the processing of this data
    if os.path.isdir(outfile+'.touch'):
        os.rmdir(outfile+'.touch')
    
    return()

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
            logger.info("... Defaulting to all SPW selections.")
            spw = ''

    # Determine the column to use

    if datacolumn is None:
        casaStuff.tb.open(infile, nomodify = True)
        colnames = casaStuff.tb.colnames()
        if 'CORRECTED_DATA' in colnames:
            logger.info("... Data has a CORRECTED column. Will use that.")
            datacolumn = 'CORRECTED'
        else:
            logger.info("... Data lacks a CORRECTED column. Will use DATA column.")
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
        start_vel_string =  ("{:12.8f}".format(vstart_kms)+'km/s').strip()
        
        # Check that we have a velocity width

        if vwidth_kms is None:
            if nchan is None:
                logger.error("Please specify a velocity width in km/s or number of channels.")
                raise Exception("No velocity width specified.")
            else:
                vwidth_kms = nchan*target_chan_kms

        # Figure out the current channel spacing
        
        line_low_ghz, line_high_ghz = lines.get_ghz_range_for_line(
            restfreq_ghz=restfreq_ghz, 
            vlow_kms=vstart_kms, vhigh_kms=vstart_kms+vwidth_kms)
        line_freq_ghz = (line_low_ghz+line_high_ghz)*0.5
        
        max_chan_ghz = np.max(np.abs(au.getChanWidths(infile, spw)))/1e9

        # Using RADIO velocity convention:
        current_chan_kms = max_chan_ghz/restfreq_ghz*sol_kms

        # Old approach using the relative/high-z convention was:
        #current_chan_kms = max_chan_ghz/line_freq_ghz*sol_kms
        
        skip_width = False
        if target_chan_kms is None:
            logger.warning('Target channel not set. Using current channel.')
            target_chan_kms = current_chan_kms
            skip_width = True
        elif current_chan_kms > target_chan_kms:
            logger.warning('Target channel less than current channel:')
            logger.warning('Asked for '+str(target_chan_kms)+' current '+str(current_chan_kms))
            target_chan_kms = current_chan_kms
            skip_width = True

        chanwidth_string =  ("{:12.8f}".format(target_chan_kms)+'km/s').strip()

        # Figure the number of channels if not supplied

        if nchan is None:
            nchan = int(np.max(np.ceil(vwidth_kms / target_chan_kms)))
        
        params.update(
            {'combinespws': False, 'regridms': True, 'chanaverage': False,
             'mode': 'velocity', 'interpolation': 'cubic', 
             'outframe': 'lsrk', 'veltype': 'radio', 'restfreq': restfreq_string, 
             'start': start_vel_string, 'nchan': nchan, 'width': chanwidth_string })
                      
        if skip_width:
            del params['width']

        message = '... regrid channel width '+chanwidth_string+' and nchan '+str(nchan)

    # ............................................
    # Rebin
    # ............................................
        
    if method == 'rebin':

        params.update({'combinespws': False, 'regridms': False, 'chanaverage' : True, 
                       'chanbin': binfactor })

        if binfactor == 1:
            params['chanaverage'] = False

        message = '... rebin by a factor of '+str(binfactor)

    # ............................................
    # Combine SPWs
    # ............................................
    
    if method == 'combine':
        
        params.update({'combinespws': True, 'regridms': False, 'chanaverage' : False, 
                       'keepflags': False })
        
        message = '... combine attempting to merge spectral windows.'

    return(params, message)

def reweight_data(
    infile = None, 
    edge_kms = None,
    edge_chans = None,
    overwrite = False, 
    datacolumn = None,
    ):
    """
    Use statwt to empirically re-weight data (mostly for
    lines). Accepts an "edge" definition in either channels or km/s.
    """
    
    # Check input

    if infile is None:
        logging.error("Please specify an input file.")
        raise Exception("Please specify an input file.")

    if not os.path.isdir(infile):
        logger.error('The input measurement set "'+infile+'"does not exist.')
        raise Exception('The input measurement set "'+infile+'"does not exist.')

    # Determine column to use

    if datacolumn is None:
        casaStuff.tb.open(infile, nomodify = True)
        colnames = casaStuff.tb.colnames()
        if 'CORRECTED_DATA' in colnames:
            logger.info("... Data has a CORRECTED column. Will use that.")
            datacolumn = 'CORRECTED'
        else:
            logger.info("... Data lacks a CORRECTED column. Will use DATA column.")
            datacolumn = 'DATA'
        casaStuff.tb.close()

    # Figure out the channel selection string

    exclude_str = ''
    if (edge_chans is not None) or (edge_kms is not None):

        vm = au.ValueMapping(infile)

        first = True
        for this_spw in vm.spwInfo.keys():

            if edge_kms is not None:
                spw_high_ghz = np.max(vm.spwInfo[this_spw]['edgeChannels'])/1e9
                spw_low_ghz = np.min(vm.spwInfo[this_spw]['edgeChannels'])/1e9
                spw_chanwidth_ghz = abs(vm.spwInfo[this_spw]['chanWidth'])/1e9

                mean_freq_ghz = 0.5*(spw_high_ghz+spw_low_ghz)

                # Here we COULD convert to RADIO convention:

                # mean_chanwidth_kms = spw_chanwidth_ghz/restfreq_ghz*sol_kms

                # BUT the rest frequency is not defined and this is an
                # approximate case, so would propose to leave this.

                # Note use of high redshift convention.

                logger.warning("... be aware that (only) statwt uses the high-z velocity convention.")
                logger.warning("... the rest of the pipleine uses radio convention.")
                mean_chanwidth_kms = spw_chanwidth_ghz/mean_freq_ghz*sol_kms
                
                edge_chans = int(np.ceil(edge_kms / mean_chanwidth_kms))

            nchan = vm.spwInfo[this_spw]['numChannels']

            if edge_chans*2 > nchan:
                logger.warning("... Too many edge channels for given spw: "+str(this_spw))
                logger.warning("... By default we will not exclude ANY channels.")
                continue

            low = int(np.ceil(edge_chans-1))
            if low < 0:
                low = 0

            high = int(np.floor(nchan-edge_chans-2))
            if high > nchan-1:
                high = nchan-1

            if high == low:
                logger.warning("Too many edge channels for given spw: "+str(this_spw))
                logger.warning("By default we will not exclude ANY channels.")
                continue

            if first:
                exclude_str += str(this_spw)+':'+str(low)+'~'+str(high)
                first=False
            else:                
                exclude_str += ','+str(this_spw)+':'+str(low)+'~'+str(high)
    
    if exclude_str != '':
        logger.info("... running statwt with exclusion: "+exclude_str)

    # Build the statwt call

    if 'fitspw' in inspect.getargspec(casaStuff.statwt)[0]:
        # CASA version somewhat >= 5.5.0
        if exclude_str == '':
            excludechans = False
        else:
            excludechans = True
        statwt_params = {'vis': infile, 'timebin': '0.001s', 'slidetimebin': False, 'chanbin': 'spw', 
                         'statalg': 'classic', 'datacolumn': datacolumn, 
                         'fitspw': exclude_str, 'excludechans': excludechans}
    else:
        # CASA version <= 5.4.1
        statwt_params = {'vis': infile, 'timebin': '0.001s', 'slidetimebin': False, 'chanbin': 'spw', 
                         'statalg': 'classic', 'datacolumn': datacolumn, 
                         'excludechans': exclude_str}
    
    # Run the call
    if not os.path.isdir(infile+'.touch'):
        os.mkdir(infile+'.touch') # no overwrite checks here
    
    logger.info("... running CASA "+'statwt('+', '.join("{!s}={!r}".format(t, statwt_params[t]) for t in statwt_params.keys())+')')
    
    casaStuff.statwt(**statwt_params)
    
    logger.info("... statwt done")
    
    if os.path.isdir(infile+'.touch'):
        os.rmdir(infile+'.touch') # no overwrite checks here

    return()

########################################
# Extract a continuum measurement set. #
########################################

def batch_extract_continuum(
    infile_list = [],
    outfile = None,
    ranges_to_exclude = [],
    lines_to_flag = [],
    vsys_kms=None, vwidth_kms=None, 
    vlow_kms=None, vhigh_kms=None,    
    do_statwt = False, 
    do_collapse = True, 
    clear_pointing = True,
    overwrite = False,
    ):
    """
    Run a batch continuum extraction.
    """

    # Check that we have an output file defined.

    if outfile is None:
        logging.error("Please specify an output file.")
        raise Exception("Please specify an output file.")

    # Check existence of output data and abort if found and overwrite is off

    if os.path.isdir(outfile) and not os.path.isdir(outfile+'.touch'):            
        if not overwrite:
            logger.warning('Found existing output data "'+outfile+'", will not overwrite it.')
            return()

    # Else, clear all previous files and temporary files

    for suffix in ['', '.flagversions', '.touch', '.temp*']:
        if not (suffix.find('*') >= 0):
            if os.path.isdir(outfile+suffix):
                shutil.rmtree(outfile+suffix)
        else:
            for temp_outfile in glob.glob(outfile+suffix):
                if os.path.isdir(temp_outfile):
                    shutil.rmtree(temp_outfile)
                    logger.debug('... shutil.rmtree(%r)'%(temp_outfile))
    
    # Create touch file to mark that we are processing this data
    if not os.path.isdir(outfile+'.touch'):
        os.mkdir(outfile+'.touch')
    
    # 
    
    split_file_list = []
    for this_infile in infile_list:

        # Specify output file and check for existence
        this_outfile = this_infile+'.temp_cont'

        split_file_list.append(this_outfile)
        
        extract_continuum(
            infile = this_infile, 
            outfile = this_outfile, 
            ranges_to_exclude = ranges_to_exclude,
            lines_to_flag = lines_to_flag,
            vsys_kms=vsys_kms, vwidth_kms=vwidth_kms, 
            vlow_kms=vlow_kms, vhigh_kms=vhigh_kms,    
            do_statwt = do_statwt,
            do_collapse = do_collapse, 
            overwrite = overwrite, 
            )

        # Deal with pointing table - testing shows it to be a
        # duplicate for each SPW here, so we remove all rows for
        # all SPWs except the first one.

        if clear_pointing:
            au.clearPointingTable(this_outfile)

    # Concatenate and combine the output data sets

    if clear_pointing:
        copy_pointing = False
    else:
        copy_pointing = True
    
    # ... concat

    concat_ms(
        infile_list = split_file_list,
        outfile = outfile,
        overwrite = overwrite,
        copypointing = copy_pointing)
    
    # ... statwt
    
    #if do_statwt:
    #    statwt_params = {'vis':outfile, 'timebin':'0.001s', 'slidetimebin':False, 'chanbin':'spw', 'statalg':'classic'}
    #    logger.info("... running CASA "+'statwt('+', '.join("{!s}={!r}".format(t, statwt_params[t]) for t in statwt_params.keys())+')')
    #    casaStuff.statwt(**statwt_params)

    # Clean up, deleting intermediate files

    for this_file in split_file_list:
        shutil.rmtree(this_file)
    
    # Remove touch file to mark that we are have done the processing of this data
    if os.path.isdir(outfile+'.touch'):
        os.rmdir(outfile+'.touch')
    
    return()

########################################
# Extract a continuum measurement set. #
########################################
                      
def extract_continuum(
    infile = None, 
    outfile = None, 
    ranges_to_exclude = [],    
    lines_to_flag = [],
    vsys_kms=None, vwidth_kms=None, 
    vlow_kms=None, vhigh_kms=None,    
    do_statwt = False, 
    do_collapse = True, 
    overwrite = False, 
    ):
    """
    Extract continuum uv data from a measurement set. Optionally
    reweights and collapses the data to a single channel per spw.
        
    Args:
    
    infile (str): The input measurement set data with suffix ".ms".
    outfile (str): The output measurement set data with suffix ".ms".
    do_statwt (bool): 
    do_collapse (bool): Always True to produce the single-channel continuum data.
    
    Inputs:
    
    Outputs:
    
    """
    
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

    # Check existing output data
    
    if os.path.isdir(outfile) and not os.path.isdir(outfile+'.touch'):
        if not overwrite:
            logger.warning('Found existing output data "'+outfile+'", will not overwrite it.')
            return()

    # Evaluate line flagging logic. If frequency ranges are supplied,
    # use those. If not but we do have line and velocity data, use
    # those to generate frequency ranges to flag.

    if len(ranges_to_exclude) == 0 and len(lines_to_flag) > 0: 
        vsys_method = (vsys_kms is not None) and (vwidth_kms is not None)
        vlow_method = (vlow_kms is not None) and (vhigh_kms is not None)

        if vsys_method or vlow_method:
            ranges_to_exclude = lines.get_ghz_range_for_list(
                line_list=lines_to_exclude, 
                vsys_kms=vsys, vwidth_kms=vwidth, vlow_kms=vlow_kms, vhigh_kms=vhigh_kms)

    # if overwrite, then delete existing output data.

    for suffix in ['', '.flagversions', '.touch',
                   '.temp', '.temp.flagversions', 
                   '.temp_copy', '.temp_copy.flagversions', 
                   ]:
        if os.path.isdir(outfile+suffix):
            shutil.rmtree(outfile+suffix)

    # find_spw_channels_for_provided frequency ranges

    spw_flagging_string = spw_string_for_freq_ranges(
        infile = infile, 
        freq_ranges_ghz = ranges_to_exclude,
        fail_on_empty = True,
        )

    # Make a copy of the data
    
    if not os.path.isdir(outfile+'.touch'):
        os.mkdir(outfile+'.touch')
    shutil.copytree(infile, outfile)
    if os.path.isdir(outfile+'.touch'):
        os.rmdir(outfile+'.touch')

    # Flag the relevant line-affected channels.
    if spw_flagging_string != '':
        if not os.path.isdir(outfile+'.touch'):
            os.mkdir(outfile+'.touch')
        casaStuff.flagdata(vis=outfile,
                           spw=spw_flagging_string,
                           )
        if os.path.isdir(outfile+'.touch'):
            os.rmdir(outfile+'.touch')

    # If requested, re-weight the data using statwt. This can be VERY
    # slow.

    # Here - this command needs to be examined and refined in CASA
    # 5.6.1 to see if it can be sped up. Right now things are
    # devastatingly slow.

    if do_statwt:

        casaStuff.tb.open(outfile, nomodify = True)
        colnames = casaStuff.tb.colnames()
        if 'CORRECTED_DATA' in colnames:
            logger.info("... Data has a CORRECTED column. Will use that.")
            datacolumn = 'CORRECTED'
        else:
            logger.info("... Data lacks a CORRECTED column. Will use DATA column.")
            datacolumn = 'DATA'
        casaStuff.tb.close()

        logger.info("... deriving empirical weights using STATWT.")
        if not os.path.isdir(outfile+'.touch'):
            os.mkdir(outfile+'.touch')
        casaStuff.statwt(vis=outfile,
                         timebin='0.001s',
                         slidetimebin=False,
                         chanbin='spw',
                         statalg='classic',
                         datacolumn=datacolumn,
                         )
        if os.path.isdir(outfile+'.touch'):
            os.rmdir(outfile+'.touch')

    # collapse

    if do_collapse:
        logger.info("... Collapsing each continuum SPW to a single channel.")
        
        if os.path.isdir(outfile):
            shutil.move(outfile, outfile+'.temp_copy')
        if os.path.isdir(outfile+'.flagversions'):
            shutil.move(outfile+'.flagversions', outfile+'.temp_copy'+'.flagversions')
        
        if not os.path.isdir(outfile+'.touch'):
            os.mkdir(outfile+'.touch')

        casaStuff.tb.open(outfile+'.temp_copy', nomodify = True)
        colnames = casaStuff.tb.colnames()
        if 'CORRECTED_DATA' in colnames:
            logger.info("... Data has a CORRECTED column. Will use that.")
            datacolumn = 'CORRECTED'
        else:
            logger.info("... Data lacks a CORRECTED column. Will use DATA column.")
            datacolumn = 'DATA'
        casaStuff.tb.close()

        casaStuff.split(vis=outfile+'.temp_copy',
                        outputvis=outfile,
                        width=100000,
                        datacolumn=datacolumn,
                        keepflags=False)
                        #<TODO><20200210># num_chan or width

        if os.path.isdir(outfile+'.touch'):
            os.rmdir(outfile+'.touch')

        # clean up
        if os.path.isdir(outfile+'.temp_copy'):
            shutil.rmtree(outfile+'.temp_copy')
        if os.path.isdir(outfile+'.temp_copy.flagversions'):
            shutil.rmtree(outfile+'.temp_copy.flagversions')

    return()

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








