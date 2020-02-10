"""
Standalone routines to analyze and manipulate visibilities.
"""

#region Imports and definitions

import os, sys, re, shutil, inspect, copy
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



def split_science_targets(
    in_file, 
    out_file, 
    do_split = True, 
    do_statwt = False, 
    use_symlink = True, 
    overwrite = False, 
    quiet = False, 
    ):
    """Split science targets from the input ALMA measurement set to a new measurement set.
    
    Args:
        in_file (str): The input measurement set data with suffix ".ms".
        
        out_file (str): The output measurement set data with suffix ".ms".
        
        do_split (bool): Set to False to only make a copy of the original measurement set without 
        splitting science targets data. The default is True, splitting science targets data. 
        
        use_symlink (bool): Set to False to make a copy of the original measurement set so as to 
        absolutely avoid any touching of the original data. The default is True, to trust 
        CASA to not modify our original data.
        
        overwrite (bool): Set to True to overwrite existing output data. The default is False, not 
        overwriting anything. 
    
    Inputs:
        in_file: ALMA measurement set data folder.
    
    Outputs:
        out_file: ALMA measurement set data folder.
    
    """
    
    # 
    # This funcion was part of the copy_data() function in older version phangsPipeline.py.
    # 
    
    # 
    # check input ms data dir
    if os.path.isdir(in_file):
        in_file = in_file
    else:
        logger.error('Error! The input uv data measurement set "'+in_file+'"does not exist!')
        raise Exception('Error! The input uv data measurement set "'+in_file+'"does not exist!')
    
    # 
    # check output suffix
    if re.match(r'^(.*)\.ms$', out_file, re.IGNORECASE):
        out_name = re.sub(r'^(.*)\.ms$', r'\1', out_file, re.IGNORECASE)
        out_file = out_name + '.ms'
    else:
        out_name = out_file
        out_file = out_name + '.ms'
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
    in_file, 
    vsys = None, 
    gal = None, 
    key_handler = None, 
    quiet = False, 
    output_spw = False, 
    ):    
    """List spectral windows in the input measurement set that contain specific spectral lines in the line_list module. 
    
    Args:
        in_file (str): The input measurement set data with suffix ".ms".
        vsys (float): Galaxy systematic velocity in units of km/s. 
        gal (str): Galaxy name, optional if vsys is given.
        key_handler (object): Our keyHandler object, optional if vsys is given, otherwise gal and key_handler should both be set. 
        output_spw (bool): Set to True to output not only found line names but also corresponding spectral window (spw) number. 
    
    Returns:
        lines_in_ms (list): A list of found line names.
        spws_in_ms (list): A list of the spectral windows (spws) corresponding to the found lines in lines_in_ms. It is a list of lists. 
    
    """
    
    # 
    # This funcion is modified from the list_lines_in_ms() function in older version phangsPipeline.py.
    # 
    
    # 
    # check input ms data dir
    if os.path.isdir(in_file):
        in_file = in_file
    else:
        logger.error('Error! The input uv data measurement set "'+in_file+'"does not exist!')
        raise Exception('Error! The input uv data measurement set "'+in_file+'"does not exist!')
    
    # 
    # if user has input gal and key_handler, find vsys
    gal_vsys = vsys
    if gal is not None and key_handler is not None:
        if gal in key_handler._target_dict:
            gal_vsys = key_handler._target_dict[gal]['vsys']
            #gal_vwidth = key_handler._target_dict[gal]['vwidth']
            if vsys is None:
                vsys = gal_vsys
            elif not np.isclose(vsys, gal_vsys):
                # if user has input a vsys, use it instead of the one in the key_handler, but report warning if the values are different
                logger.warning('Warning! User has input a vsys of '+str(vsys)+' km/s which is different from the vsys of '+gal_vsys+' km/s for the galaxy "'+gal+'" in the key_handler as defined in "target_definitions.txt"!')
        else:
            logger.error('Error! Could not find the input galaxy "'+gal+'" in the key_handler as defined in "target_definitions.txt"!')
            raise Exception('Error! Could not find the input galaxy "'+gal+'" in the key_handler as defined in "target_definitions.txt"!')
    # 
    # check vsys
    if vsys is None:
        logger.error('Error! Please input a `vsys` for the galaxy systematic velocity in units of km/s, or input a `gal` and `key_handler`.')
        raise Exception('Error! Please input vsys, or gal and key_handler.')
    # 
    # find lines in the ms data
    lines_in_ms = None
    spws_in_ms = None
    logger.info('Finding spectral windows for '+str(len(line_list.line_list))+' lines in our line_list module')
    for line in line_list.line_list.keys():
        # 
        # get line rest-frame frequencies
        restfreq_ghz = line_list.line_list[line]
        # 
        # work out the frequency of the line and the line wings
        target_freq_ghz = restfreq_ghz*(1.-vsys/sol_kms)
        # 
        # run analysisUtils.getScienceSpwsForFrequency() to get the corresponding spectral window (spw)
        this_spw_list = au.getScienceSpwsForFrequency(in_file, target_freq_ghz*1e9)
        # 
        if len(this_spw_list) == 0:
            continue
        # 
        if lines_in_ms is None:
            lines_in_ms = []
        if spws_in_ms is None:
            spws_in_ms = []
        lines_in_ms.append(line)
        spws_in_ms.append(this_spw_list)
    # 
    logger.info('Found '+str(len(lines_in_ms))+' lines from our line_list module that are within the spectral windows of "'+in_file+'"')
    # 
    if output_spw:
        return lines_in_ms, spws_in_ms
    else:
        return lines_in_ms




def chanwidth_for_line(
    in_file, 
    line, 
    vsys = None, 
    vwidth = None, 
    gal = None, 
    key_handler = None, 
    quiet = False, 
    ): 
    """Calculates the coarsest channel width among all spectral windows in the input measurement set that contain the input line.
    
    Args:
        in_file (str): The input measurement set data with suffix ".ms".
        line (str): Line name. 
        gal (str): Galaxy name, optional if vsys is given.
        key_handler (object): Our keyHandler object, optional if vsys is given, otherwise gal and key_handler should both be set. 
        output_spw (bool): Set to True to output not only found line names but also corresponding spectral window (spw) number. 
    
    Returns:
        chan_width_kms (float): 
    
    """
    
    # 
    # This funcion is modified from the chanwidth_for_line() function in older version phangsPipeline.py.
    # 
    
    # 
    # check input ms data dir
    if os.path.isdir(in_file):
        in_file = in_file
    else:
        logger.error('Error! The input uv data measurement set "'+in_file+'"does not exist!')
        raise Exception('Error! The input uv data measurement set "'+in_file+'"does not exist!')
    
    # 
    # if user has input gal and key_handler, find vsys
    gal_vsys = vsys
    if gal is not None and key_handler is not None:
        if gal in key_handler._target_dict:
            gal_vsys = key_handler._target_dict[gal]['vsys']
            gal_vwidth = key_handler._target_dict[gal]['vwidth']
            if vsys is None:
                vsys = gal_vsys
            elif not np.isclose(vsys, gal_vsys):
                # if user has input a vsys, use it instead of the one in the key_handler, but report warning if the values are different
                logger.warning('Warning! User has input a vsys of '+str(vsys)+' km/s which is different from the vsys of '+str(gal_vsys)+' km/s for the galaxy "'+gal+'" in the key_handler as defined in "target_definitions.txt"!')
            if vwidth is None:
                vwidth = gal_vwidth
            elif not np.isclose(vwidth, gal_vwidth):
                # if user has input a vwidth, use it instead of the one in the key_handler, but report warning if the values are different
                logger.warning('Warning! User has input a vwidth of '+str(vwidth)+' km/s which is different from the vwidth of '+str(gal_vwidth)+' km/s for the galaxy "'+gal+'" in the key_handler as defined in "target_definitions.txt"!')
        else:
            logger.error('Error! Could not find the input galaxy "'+gal+'" in the key_handler as defined in "target_definitions.txt"!')
            raise Exception('Error! Could not find the input galaxy "'+gal+'" in the key_handler as defined in "target_definitions.txt"!')
    # 
    # check vsys
    if vsys is None:
        logger.error('Error! Please input a `vsys` for the galaxy systematic velocity in units of km/s, or input a `gal` and `key_handler`.')
        raise Exception('Error! Please input vsys, or gal and key_handler.')
    # 
    # check vwidth
    if vwidth is None:
        vwidth = 500.0
        logger.info('Setting line width `vwidth` to 500 km/s as the default value. Please consider setting a better value in "target_definitions.txt" and input gal and key_handler to this function.')
    # 
    # try to match the input line in the line_list module
    matched_line = None
    if matched_line is None:
        if line in line_list.line_list:
            matched_line = line
    if matched_line is None:
        if line.lower() in line_list.line_list:
            matched_line = line.lower()
    if matched_line is None:
        logger.error('Error! Could not find the input line "'+line+'" in our line_list module. Candiate line names are: '+str(line_list.line_list.keys()))
        raise Exception('Error! Could not find the input line "'+line+'" in our line_list module. Candiate line names are: '+str(line_list.line_list.keys()))
    # 
    # get line center rest-frame frequencies
    restfreq_ghz = line_list.line_list[line]
    # 
    # Work out which spectral windows contain the line contain
    target_freq_ghz = restfreq_ghz*(1.-vsys/sol_kms)
    target_freq_high = restfreq_ghz*(1.-(vsys-0.5*vwidth)/sol_kms)
    target_freq_low = restfreq_ghz*(1.-(vsys+0.5*vwidth)/sol_kms)
    # 
    spw_list = []
    # 
    # loop over line left wing, center and right wing frequencies to match the spectral window (spw)
    for target_freq in [target_freq_high, target_freq_ghz, target_freq_low]:
        # run analysisUtils.getScienceSpwsForFrequency() to get the corresponding spectral window (spw)
        this_spw_list = au.getScienceSpwsForFrequency(in_file, target_freq*1e9)    
        spw_list.extend(this_spw_list)
    # 
    if len(spw_list) == 0:
        logger.error('No spectral windows contain the input line "'+line+'" with vwidth '+str(vwidth)+' km/s at a vsys of '+str(vsys)+' km/s. The ms data is "'+in_file+'".')
        return None
    # 
    # sort and remove duplicates
    spw_list = sorted(list(set(spw_list)))
    # 
    # make spw_list_string
    spw_list_string = ','.join(np.array(spw_list).astype(str))
    # 
    # Figure out how much averaging is needed to reach the target resolution
    chan_width_hz = au.getChanWidths(in_file, spw_list_string)
    # 
    # Convert to km/s and return
    chan_width_kms = abs(chan_width_hz / (restfreq_ghz*1e9)*sol_kms)
    # 
    # Return
    return chan_width_kms



def extract_line(
    in_file, 
    out_file, 
    line = 'co21', 
    vsys = None, 
    vwidth = None, 
    gal = None, 
    key_handler = None, 
    chan_fine = 0.5, 
    rebin_factor = 5, 
    do_statwt = False, 
    edge_for_statwt = -1, 
    quiet = False, 
    overwrite = False, 
    ):
    """Extract a spectral line uv data from a measurement set with optimized regridding. 
    
    Extract a spectral line from a measurement set and regrid onto a
    new velocity grid with the desired spacing. This doesn't
    necessarily need the PHANGS keys in place and may be a general
    purpose utility. There are some minor subtleties here related to
    regridding and rebinning.
    
    Args:
        in_file (str): The input measurement set data with suffix ".ms".
        out_file (str): The output measurement set data with suffix ".ms".
        line (str): Line name. 
        gal (str): Galaxy name, optional if vsys and vwidth are given. 
        key_handler (object): Our keyHandler object, optional if vsys and vwidth are given, otherwise gal and key_handler should both be set. 
        chan_fine (float): 
        rebin_factor (float): 
        do_statwt (bool): 
        edge_for_statwt (int): 
        quiet (bool): True for quiet and False for verbose. Default is False. 
    
    Inputs:
        in_file: ALMA measurement set data folder.
    
    Outputs:
        out_file: ALMA measurement set data folder.
    
    """
    
    # 
    # This funcion is modified from the extract_line() function in older version phangsPipeline.py.
    # 
    
    # 
    # check input ms data dir
    if os.path.isdir(in_file):
        in_file = in_file
    else:
        logger.error('Error! The input uv data measurement set "'+in_file+'"does not exist!')
        raise Exception('Error! The input uv data measurement set "'+in_file+'"does not exist!')
    
    # 
    # check output suffix
    if re.match(r'^(.*)\.ms$', out_file, re.IGNORECASE):
        out_name = re.sub(r'^(.*)\.ms$', r'\1', out_file, re.IGNORECASE)
        out_file = out_name + '.ms'
    else:
        out_name = out_file
        out_file = out_name + '.ms'
    # 
    # check existing output data in the imaging directory
    if os.path.isdir(out_file):
        if not overwrite:
            logger.warning('Found existing output data "'+out_file+'", will not overwrite it.')
            return
        else:
            shutil.rmtree(out_file)
            if os.path.isdir(out_file+'.flagversions'):
                shutil.rmtree(out_file+'.flagversions')
            if os.path.isdir(out_file+'.temp'):
                shutil.rmtree(out_file+'.temp')
            if os.path.isdir(out_file+'.temp.flagversions'):
                shutil.rmtree(out_file+'.temp.flagversions')
            if os.path.isdir(out_file+'.temp2'):
                shutil.rmtree(out_file+'.temp2')
            if os.path.isdir(out_file+'.temp2.flagversions'):
                shutil.rmtree(out_file+'.temp2.flagversions')
    # 
    # if user has input gal and key_handler, find vsys and vwidth
    gal_vsys = vsys
    if gal is not None and key_handler is not None:
        if gal in key_handler._target_dict:
            gal_vsys = key_handler._target_dict[gal]['vsys']
            gal_vwidth = key_handler._target_dict[gal]['vwidth']
            if vsys is None:
                vsys = gal_vsys
            elif not np.isclose(vsys, gal_vsys):
                # if user has input a vsys, use it instead of the one in the key_handler, but report warning if the values are different
                logger.warning('Warning! User has input a vsys of '+str(vsys)+' km/s which is different from the vsys of '+str(gal_vsys)+' km/s for the galaxy "'+gal+'" in the key_handler as defined in "target_definitions.txt"!')
            if vwidth is None:
                vwidth = gal_vwidth
            elif not np.isclose(vwidth, gal_vwidth):
                # if user has input a vwidth, use it instead of the one in the key_handler, but report warning if the values are different
                logger.warning('Warning! User has input a vwidth of '+str(vwidth)+' km/s which is different from the vwidth of '+str(gal_vwidth)+' km/s for the galaxy "'+gal+'" in the key_handler as defined in "target_definitions.txt"!')
        else:
            logger.error('Error! Could not find the input galaxy "'+gal+'" in the key_handler as defined in "target_definitions.txt"!')
            raise Exception('Error! Could not find the input galaxy "'+gal+'" in the key_handler as defined in "target_definitions.txt"!')
    # 
    # check vsys
    if vsys is None:
        logger.error('Error! Please input a `vsys` for the galaxy systematic velocity in units of km/s, or input a `gal` and `key_handler`.')
        raise Exception('Error! Please input vsys, or gal and key_handler.')
    # 
    # check vwidth
    if vwidth is None:
        vwidth = 500.0
        logger.info('Setting line width `vwidth` to 500 km/s as the default value. Please consider setting a better value in "target_definitions.txt" and input gal and key_handler to this function.')
    # 
    # try to match the input line in the line_list module
    matched_line = None
    if matched_line is None:
        if line in line_list.line_list:
            matched_line = line
    if matched_line is None:
        if line.lower() in line_list.line_list:
            matched_line = line.lower()
    if matched_line is None:
        logger.error('Error! Could not find the input line "'+line+'" in our line_list module. Candiate line names are: '+str(line_list.line_list.keys()))
        raise Exception('Error! Could not find the input line "'+line+'" in our line_list module. Candiate line names are: '+str(line_list.line_list.keys()))
    # 
    # get line center rest-frame frequencies
    restfreq_ghz = line_list.line_list[line]
    # 
    # Work out which spectral windows contain the line contain
    target_freq_ghz = restfreq_ghz*(1.-vsys/sol_kms)
    target_freq_high = restfreq_ghz*(1.-(vsys-0.5*vwidth)/sol_kms)
    target_freq_low = restfreq_ghz*(1.-(vsys+0.5*vwidth)/sol_kms)
    # 
    # ... Below are the same as in the above function chanwidth_for_line() ...
    # 
    spw_list = []
    # 
    # loop over line left wing, center and right wing frequencies to match the spectral window (spw)
    for target_freq in [target_freq_high, target_freq_ghz, target_freq_low]:
        # run analysisUtils.getScienceSpwsForFrequency() to get the corresponding spectral window (spw)
        this_spw_list = au.getScienceSpwsForFrequency(in_file, target_freq*1e9)    
        spw_list.extend(this_spw_list)
    # 
    if len(spw_list) == 0:
        logger.error('No spectral windows contain the input line "'+line+'" with vwidth '+str(vwidth)+' km/s at a vsys of '+str(vsys)+' km/s. The ms data is "'+in_file+'".')
        return None
    # 
    # sort and remove duplicates
    spw_list = sorted(list(set(spw_list)))
    # 
    # make spw_list_string
    spw_list_string = ','.join(np.array(spw_list).astype(str))
    # 
    # ... Above are the same as in the above function chanwidth_for_line() ...
    # 
    # print starting message
    if not quiet:
        logger.info("--------------------------------------")
        logger.info("EXTRACT_LINE begins:")
    # 
    if not quiet:
        logger.info("... line: "+line)
        logger.info("... spectral windows to consider: "+spw_list_string)
    # 
    # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
    # STEP 1. Shift the zero point AND change the channel width (slightly).
    # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
    # 
    start_vel_kms = (vsys - vwidth/2.0)
    chan_width_hz = au.getChanWidths(in_file, spw_list_string)
    current_chan_width_kms = abs(chan_width_hz / (restfreq_ghz*1e9)*sol_kms)        
    if chan_fine == -1:
        nchan_for_recenter = int(np.max(np.ceil(vwidth / current_chan_width_kms)))
    else:
        nchan_for_recenter = int(np.max(np.ceil(vwidth / chan_fine)))
    # 
    # Cast to text with specified precision.
    restfreq_string = "{:12.8f}".format(restfreq_ghz)+'GHz'
    start_vel_string =  "{:12.8f}".format(start_vel_kms)+'km/s'
    chanwidth_string =  "{:12.8f}".format(chan_fine)+'km/s'
    # 
    if not quiet:
        logger.info("... shifting the fine grid (before any regridding)")
        logger.info("... rest frequency: "+restfreq_string)
        logger.info("... new starting velocity: "+start_vel_string)
        logger.info("... original velocity width: "+str(current_chan_width_kms))
        logger.info("... target velocity width: "+str(chan_fine))
        logger.info("... number of channels at this stage: "+str(nchan_for_recenter))
    # 
    #os.system('rm -rf '+out_file+'.temp')
    #os.system('rm -rf '+out_file+'.temp.flagversions')
    # 
    if chan_fine == -1:
        casaStuff.mstransform(vis=in_file,
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
        casaStuff.mstransform(vis=in_file,
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

    if not quiet:
        logger.info("... channel averaging")
        logger.info("... rebinning factor: "+str(rebin_factor))

    if rebin_factor > 1:
        #os.system('rm -rf '+out_file+'.temp2')
        #os.system('rm -rf '+out_file+'.temp2.flagversions')
        casaStuff.mstransform(vis=current_file,
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
    
    if not quiet:
        logger.info("... combining spectral windows")

    #os.system('rm -rf '+out_file)
    #os.system('rm -rf '+out_file+'.flagversions')
    casaStuff.mstransform(vis=current_file,
                          outputvis=out_file,
                          spw='',
                          datacolumn='DATA',
                          regridms=False,
                          chanaverage=False,
                          combinespws=True
                          )    
    
    if not quiet:
        logger.info("... deleting old files")
    
    # Clean up
    #os.system('rm -rf '+out_file+'.temp')
    #os.system('rm -rf '+out_file+'.temp.flagversions')
    #os.system('rm -rf '+out_file+'.temp2')
    #os.system('rm -rf '+out_file+'.temp2.flagversions')
    if os.path.isdir(out_file+'.temp'):
        shutil.rmtree(out_file+'.temp')
    if os.path.isdir(out_file+'.temp.flagversions'):
        shutil.rmtree(out_file+'.temp.flagversions')
    if os.path.isdir(out_file+'.temp2'):
        shutil.rmtree(out_file+'.temp2')
    if os.path.isdir(out_file+'.temp2.flagversions'):
        shutil.rmtree(out_file+'.temp2.flagversions')
    
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
        
        logger.info("... running statwt with exclusion: "+exclude_str)
        
        # This needs to revert to oldstatwt, it seems not to work in the new form
        
        if 'fitspw' in inspect.getargspec(casaStuff.statwt)[0]:
            # CASA version somewhat >= 5.5.0
            test = casaStuff.statwt(vis=out_file,
                                    timebin='0.001s',
                                    slidetimebin=False,
                                    chanbin='spw',
                                    statalg='classic',
                                    datacolumn='data',
                                    fitspw=exclude_str,
                                    excludechans=True,
                                    )
        else:
            # CASA version <= 5.4.1
            test = casaStuff.statwt(vis=out_file,
                                    timebin='0.001s',
                                    slidetimebin=False,
                                    chanbin='spw',
                                    statalg='classic',
                                    datacolumn='data',
                                    excludechans=exclude_str,
                                    )
    
    logger.info("--------------------------------------")
    
    return



def concat_ms(
    in_file_list,  
    out_file,  
    do_chan0 = False, 
    quiet = False, 
    overwrite = False, 
    ):
    """Concatenate a list of measurement sets into one measurement set, and collapse channel0 if required. 
    
    Args:
        in_file_list (list or str): The input list of measurement sets. 
        out_file (str): The output measurement set data with suffix ".ms".
        do_chan0 (bool): True to collapse channel0.
        quiet (bool): True for quiet and False for verbose. Default is False. 
    
    Inputs:
        in_file: ALMA measurement set data folder.
    
    Outputs:
        out_file: ALMA measurement set data folder.
    
    """
    
    # 
    # This funcion is paritially replacing the call of concat() in the concat_line_for_one_gal() and concat_cont_for_one_gal() functions in older version phangsPipeline.py.
    # 
    
    # 
    # make sure the input in_file_list is a list
    if np.isscalar(in_file_list):
        in_file_list = [in_file_list]
    
    # 
    # check input ms data dir
    for in_file in in_file_list:
        if not os.path.isdir(in_file):
            logger.error('Error! The input uv data measurement set "'+in_file+'"does not exist!')
            raise Exception('Error! The input uv data measurement set "'+in_file+'"does not exist!')
    
    # 
    # check output suffix
    if re.match(r'^(.*)\.ms$', out_file, re.IGNORECASE):
        out_name = re.sub(r'^(.*)\.ms$', r'\1', out_file, re.IGNORECASE)
        out_file = out_name + '.ms'
    else:
        out_name = out_file
        out_file = out_name + '.ms'
    # 
    # check existing output data in the imaging directory
    if os.path.isdir(out_file):
        if not overwrite:
            logger.warning('Found existing output data "'+out_file+'", will not overwrite it.')
            return
    # 
    # remove existing data 
    if os.path.isdir(out_file):
        shutil.rmtree(out_file)
    if os.path.isdir(out_file+'.flagversions'):
        shutil.rmtree(out_file+'.flagversions')
    if os.path.isdir(out_file+'_chan0'):
        shutil.rmtree(out_file+'_chan0')
    if os.path.isdir(out_file+'_chan0.flagversions'):
        shutil.rmtree(out_file+'_chan0.flagversions')
    # 
    # 
    # Concatenate all of the relevant files
    casaStuff.concat(vis = in_file_list, 
                     concatvis = out_file)
                     #<TODO># what about freqtol?
    # 
    # Collapse to form a "channel 0" measurement set
    if do_chan0:
        # 
        chan0_vis = out_name + '_chan0.ms'
        if os.path.isdir(chan0_vis):
            logger.info('Found existing output data "'+chan0_vis+'". Overwriting it.')
            shutil.rmtree(chan0_vis)
        if os.path.isdir(chan0_vis+'.flagversions'):
            shutil.rmtree(chan0_vis+'.flagversions')
        # 
        # find out the channel number
        casaStuff.tb.open(out_file+os.sep+'SPECTRAL_WINDOW', nomodify = True)
        chan_width_list = casaStuff.tb.getcell('CHAN_WIDTH', 0)
        num_chan = len(chan_width_list)
        logger.info('Collapsing '+str(num_chan)+' channels in "'+os.path.basename(out_file)+'" to make "'+os.path.basename(chan0_vis)+'"')
        casaStuff.tb.close()
        # 
        if num_chan < 1:
            logger.error('Error! Could not get channel number from "'+out_file+os.sep+'SPECTRAL_WINDOW'+'"!')
            raise Exception('Error! Could not get channel number from "'+out_file+os.sep+'SPECTRAL_WINDOW'+'"!')
        # 
        casaStuff.split(vis = out_file, 
                        datacolumn = 'DATA', 
                        spw = '', 
                        outputvis = chan0_vis, 
                        width = num_chan)
    

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
    in_file, 
    out_file, 
    lines_to_flag = None, 
    vsys = None, 
    vwidth = None, 
    gal = None, 
    key_handler = None, 
    do_statwt = False, 
    do_collapse = True, 
    quiet = False, 
    overwrite = False, 
    ):
    """Extract continuum uv data from a measurement set and collapse into a single channel. 
    
    Extract a continuum measurement set, flagging any specified lines,
    reweighting using statwt, and then collapsing to a single "channel
    0" measurement.
    
    Args:
        in_file (str): The input measurement set data with suffix ".ms".
        out_file (str): The output measurement set data with suffix ".ms".
        lines_to_flag (list): A list of line names to flag. Lines names must be in our line_list module. If it is None, then we use all 12co, 13co and c18o lines.
        gal (str): Galaxy name, optional if vsys and vwidth are given. 
        key_handler (object): Our keyHandler object, optional if vsys and vwidth are given, otherwise gal and key_handler should both be set. 
        do_statwt (bool): 
        do_collapse (bool): Always True.
        quiet (bool): True for quiet and False for verbose. Default is False. 
    
    Inputs:
        in_file: ALMA measurement set data folder.
    
    Outputs:
        out_file: ALMA measurement set data folder.
    
    """
    
    # 
    # This funcion is modified from the extract_continuum() function in older version phangsPipeline.py.
    # 
    
    # 
    # check input ms data dir
    if os.path.isdir(in_file):
        in_file = in_file
    else:
        logger.error('Error! The input uv data measurement set "'+in_file+'"does not exist!')
        raise Exception('Error! The input uv data measurement set "'+in_file+'"does not exist!')
    
    # 
    # check output suffix
    if re.match(r'^(.*)\.ms$', out_file, re.IGNORECASE):
        out_name = re.sub(r'^(.*)\.ms$', r'\1', out_file, re.IGNORECASE)
        out_file = out_name + '.ms'
    else:
        out_name = out_file
        out_file = out_name + '.ms'
    # 
    # check existing output data in the imaging directory
    if os.path.isdir(out_file):
        if not overwrite:
            logger.warning('Found existing output data "'+out_file+'", will not overwrite it.')
            return
    if os.path.isdir(out_file):
        shutil.rmtree(out_file)
    if os.path.isdir(out_file+'.flagversions'):
        shutil.rmtree(out_file+'.flagversions')
    if os.path.isdir(out_file+'.temp'):
        shutil.rmtree(out_file+'.temp')
    if os.path.isdir(out_file+'.temp.flagversions'):
        shutil.rmtree(out_file+'.temp.flagversions')
    if os.path.isdir(out_file+'.temp2'):
        shutil.rmtree(out_file+'.temp2')
    if os.path.isdir(out_file+'.temp2.flagversions'):
        shutil.rmtree(out_file+'.temp2.flagversions')
    # 
    # if user has input gal and key_handler, find vsys and vwidth
    gal_vsys = vsys
    if gal is not None and key_handler is not None:
        if gal in key_handler._target_dict:
            gal_vsys = key_handler._target_dict[gal]['vsys']
            gal_vwidth = key_handler._target_dict[gal]['vwidth']
            if vsys is None:
                vsys = gal_vsys
            elif not np.isclose(vsys, gal_vsys):
                # if user has input a vsys, use it instead of the one in the key_handler, but report warning if the values are different
                logger.warning('Warning! User has input a vsys of '+str(vsys)+' km/s which is different from the vsys of '+str(gal_vsys)+' km/s for the galaxy "'+gal+'" in the key_handler as defined in "target_definitions.txt"!')
            if vwidth is None:
                vwidth = gal_vwidth
            elif not np.isclose(vwidth, gal_vwidth):
                # if user has input a vwidth, use it instead of the one in the key_handler, but report warning if the values are different
                logger.warning('Warning! User has input a vwidth of '+str(vwidth)+' km/s which is different from the vwidth of '+str(gal_vwidth)+' km/s for the galaxy "'+gal+'" in the key_handler as defined in "target_definitions.txt"!')
        else:
            logger.error('Error! Could not find the input galaxy "'+gal+'" in the key_handler as defined in "target_definitions.txt"!')
            raise Exception('Error! Could not find the input galaxy "'+gal+'" in the key_handler as defined in "target_definitions.txt"!')
    # 
    # check vsys
    if vsys is None:
        logger.error('Error! Please input a `vsys` for the galaxy systematic velocity in units of km/s, or input a `gal` and `key_handler`.')
        raise Exception('Error! Please input vsys, or gal and key_handler.')
    # 
    # check vwidth
    if vwidth is None:
        vwidth = 500.0
        logger.info('Setting line width `vwidth` to 500 km/s as the default value. Please consider setting a better value in "target_definitions.txt" and input gal and key_handler to this function.')
    # 
    # set the list of lines to flag
    if lines_to_flag is None:
        lines_to_flag = line_list.line_families['co'] + line_list.line_families['13co'] + line_list.line_families['c18o']
    else:
        lines_to_flag_copied = copy.copy(lines_to_flag)
        lines_to_flag = []
        for line_to_flag_copied in lines_to_flag_copied:
            if line_to_flag_copied.lower() == 'co':
                lines_to_flag.extend(line_list.line_families['co'])
            elif line_to_flag_copied.lower() == '12co':
                lines_to_flag.extend(line_list.line_families['co'])
            elif line_to_flag_copied.lower() == '13co':
                lines_to_flag.extend(line_list.line_families['13co'])
            elif line_to_flag_copied.lower() == 'c18o':
                lines_to_flag.extend(line_list.line_families['c18o'])
            elif line_to_flag_copied.lower() in line_list.line_families:
                lines_to_flag.extend(line_list.line_families[line_to_flag_copied.lower()])
            elif line_to_flag_copied.lower() in line_list.line_list:
                lines_to_flag.append(line_to_flag_copied.lower())
            else:
                logger.error('Error! The input line "'+line_to_flag_copied+'" in the lines_to_flag array could not be found in our line_list module!')
                raise Exception('Error! The input line "'+line_to_flag_copied+'" in the lines_to_flag array could not be found in our line_list module!')
    logger.debug('Lines to flag for gal '+gal+' continuum: '+str(lines_to_flag))
    # 
    # Make a continuum copy of the data
    command = 'cp -r -H '+in_file+' '+out_file
    logger.info(command)
    var = os.system(command)
    logger.info(var)
    # 
    # Figure out the line channels and flag them
    vm = au.ValueMapping(out_file)
    # 
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
        if len(spw_list) == 0:
            continue
        
        logger.info("Found overlap for "+line)
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
    
    logger.info("... proposed flagging "+spw_flagging_string)
    
    if spw_flagging_string != '':
        casaStuff.flagdata(vis=out_file,
                           spw=spw_flagging_string,
                           )
        
    # Here - this comman needs to be examined and refined in CASA
    # 5.6.1 to see if it can be sped up. Right now things are
    # devastatingly slow.
    if do_statwt:
        logger.info("... deriving empirical weights using STATWT.")
        casaStuff.statwt(vis=out_file,
                         timebin='0.001s',
                         slidetimebin=False,
                         chanbin='spw',
                         statalg='classic',
                         datacolumn='data',
                         )
    
    if do_collapse:
        logger.info("... Collapsing the continuum to a single channel.")
        
        os.system('rm -rf '+out_file+'.temp_copy')
        os.system('rm -rf '+out_file+'.temp_copy.flagversions')
        
        command = 'mv '+out_file+' '+out_file+'.temp_copy'
        logger.info(command)
        var = os.system(command)
        logger.info(var)

        command = 'mv '+out_file+'.flagversions '+out_file+'.temp_copy.flagversions'
        logger.info(command)
        var = os.system(command)
        logger.info(var)
        
        #<TODO><20200210># num_chan
        
        casaStuff.split(vis=out_file+'.temp_copy',
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
