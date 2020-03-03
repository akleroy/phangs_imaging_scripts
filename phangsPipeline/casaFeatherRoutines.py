"""
Standalone routines that relate to single dish-interferometer
combination using CASA's feather.
"""

#region Imports and definitions

import os
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

# Other pipeline stuff
import casaCubeRoutines as ccr

# Pipeline versionining
from pipelineVersion import version as pipeVer

#endregion

#region Feathering and single dish routines

def prep_sd_for_feather(
    sdfile_in=None,
    sdfile_out=None,
    interf_file=None,
    do_import=True,
    do_dropdeg=True,
    do_align=True,
    do_checkunits=True,
    overwrite=False):
    """
    Prepare single dish data for feathering. 

    sdfile_in : the input single dish file. Can be a FITS file (with
    do_import) or a CASA image.

    sdfile_out : the output of the program. If all flags are called,
    this will be a CASA image on the same astrometric grid as the
    interferometric file with units of Jy/beam and no degenerate axes.

    interf_file : the interferometric file being used in the feather
    call. Used as a template here to astrometrically align the data.

    do_import (default True) : If True and the infile is a FITS file,
    then import the data from FITS.

    do_dropdeg (default True) : if True then pare degenerate axes from
    the signle dish file. In general this is a good idea for
    postprocessing, where mixing and matching degenerate axes causes
    many CASA routines to fail.

    do_align (default True) : If True then align the single dish data to the
    astrometric grid of the the interferometer file.

    do_checkunits (default True) : If True then check the units to make
    sure that they are in Jy/beam, which is required by feather.

    overwrite (default False) : Delete existing files. You probably
    want to set this to True but it's a user decision.
    """

    # Check inputs
    
    if (os.path.isdir(interf_file) == False):
        logger.error("Interferometric file not found: "+interf_file)
        return(None)

    if (os.path.isdir(sdfile_in) == False) and (os.path.isfile(sdfile_in) == False):
        logger.error("Single dish file not found: "+sdfile_in)
        return(None)

    if sdfile_out is None:
        logger.error("Output single dish file name not supplied via sdfile_out=")
        return(None)

    # Initialize file handling

    current_infile = sdfile_in    
    current_outfile = sdfile_out
    tempfile_name = sdfile_out+'.temp'
    
    # Import from FITS if needed. Keep blanks as not-a-numbers.
    
    if do_import:
        if ((current_infile[-4:] == 'FITS') or \
                (current_infile[-4:] == 'fits')) and \
                os.path.isfile(current_infile):
            logger.info("Importing from FITS.")
                    
            casaStuff.importfits(
                fitsimage=current_infile, 
                imagename=current_outfile,
                zeroblanks=False,
                overwrite=overwrite)
            current_infile = current_outfile

    # Drop degenerate axes using a call to imsubimage
            
    if do_dropdeg:
        if current_infile == current_outfile:
            if os.path.isdir(tempfile_name) or os.path.isfile(tempfile_name):
                if overwrite:
                    os.system('rm -rf '+tempfile_name)
                else:
                    logger.error("Temp file needed but exists and overwrite=False - "+tempfile_name)
                    return(None)
            os.system('cp -r '+current_infile+' '+tempfile_name)
            current_infile = tempfile_name            
            os.system('rm -rf '+current_outfile)

        if os.path.isdir(current_outfile) or os.path.isfile(current_outfile):
            if overwrite:
                os.system('rm -rf '+current_outfile)
            else:
                logger.error("Output file needed exists and overwrite=False - "+current_outfile)
                return(None)

        casaStuff.imsubimage(
            imagename=current_infile, 
            outfile=current_outfile,
            dropdeg=True)

        current_infile = current_outfile

    # Align the single dish data to the astrometric grid of the interferometric data

    if do_align:
        if current_infile == current_outfile:
            if os.path.isdir(tempfile_name) or os.path.isfile(tempfile_name):
                if overwrite:
                    os.system('rm -rf '+tempfile_name)
                else:
                    logger.error("Temp file needed but exists and overwrite=False - "+tempfile_name)
                    return(None)
            os.system('cp -r '+current_infile+' '+tempfile_name)
            current_infile = tempfile_name            
            os.system('rm -rf '+current_outfile)

        casaStuff.imregrid(
            imagename=current_infile,
            template=interf_file,
            output=current_outfile,
            asvelocity=True,
            axes=[-1],
            interpolation='cubic',
            overwrite=overwrite)

        current_infile = current_outfile

    # Check the units on the singledish file and convert from K to Jy/beam if needed.

    if do_checkunits:
        hdr = casaStuff.imhead(current_outfile, mode='list')
        unit = hdr['bunit']
        if unit == 'K':
            logger.info("Unit is Kelvin. Converting.")
            ccr.convert_ktojy(
                infile=current_outfile, 
                overwrite=overwrite, 
                inplace=True)

    # Remove leftover temporary files.

    if (os.path.isdir(tempfile_name) or os.path.isfile(tempfile_name)):
        if overwrite:
            os.system('rm -rf '+tempfile_name)

    return(None)

def feather_two_cubes(
    interf_file=None,
    sd_file=None,
    out_file=None,
    do_blank=False,
    do_apodize=False,
    apod_file=None,
    apod_cutoff=-1.0,
    overwrite=False,
    ):
    """
    Feather together interferometric and total power data using CASA's
    default approach. Optionally, first apply some steps to homogenize
    the two data sets. Assumes that the data have been prepared e.g.,
    using the prep_sd_for_feather routine in this module.

    interf_file : the interferometric cube to feather.

    sd_file : the single dish cube to feather.

    out_file : the output file name

    do_blank (default False) : if True then blank masked and not-a-number
    values in the cubes. The idea is to make sure missing data are
    zero for the FFT (probably not an issue) and to make sure some
    regions don't appear in one map but not the other. This is
    probably not necessary, but CASA's treatment in feather is a bit
    unclear and the case of, e.g., a much more extended single dish
    map compared to the interferometer map comes up.

    do_apodize (default False) : if True, then apodize BOTH data sets
    using the provided apodization map.

    apod_file : if do_apodize is True, this file is the map used to scale
    both the interferometer and single dish data.

    apod_cutoff : the cutoff in the apodization file below which data
    are blanked.

    overwrite (default False) : Delete existing files. You probably
    want to set this to True but it's a user decision.
    """

    # Check inputs

    if (os.path.isdir(sd_file) == False):
        logger.error("Single dish file not found: "+sd_file)
        return(False)
        
    if (os.path.isdir(interf_file) == False):
        logger.error("Interferometric file not found: "+interf_file)
        return(False)

    if do_apodize:
        if (os.path.isdir(apod_file) == False):
            logger.error("Apodization requested, but file not found: "+apod_file)
            return(False)

    # Initialize file handling

    os.system('rm -rf '+sd_file+'.temp')
    os.system('rm -rf '+interf_file+'.temp')
    os.system('rm -rf '+out_file+'.temp')
    
    os.system('rm -rf '+sd_file+'.temp.temp')
    os.system('rm -rf '+interf_file+'.temp.temp')
    os.system('rm -rf '+out_file+'.temp.temp')

    # If requested, manipulate blanked (NaN and mask) values to make
    # sure they are zeros. This should probably not be necessary, but
    # some aspects of CASA's procedures are unclear. The main thing
    # here is that the mask used is the COMBINED single dish and
    # interferometer map, so that only regions in common should
    # survive.

    if do_blank:        
        
        current_interf_file = interf_file+'.temp'
        current_sd_file = sd_file+'.temp'

        os.system('cp -rf '+interf_file+' '+current_interf_file)
        os.system('cp -rf '+sd_file+' '+current_sd_file)

        myia = au.createCasaTool(casaStuff.iatool)

        myia.open(interf_file)
        interf_mask = myia.getchunk(getmask=True)
        myia.close()

        myia.open(sd_file)
        sd_mask = myia.getchunk(getmask=True)
        myia.close()
        
        # CASA calls unmasked values True and masked values False. The
        # region with values in both cubes is the product.

        combined_mask = sd_mask*interf_mask

        # This isn't a great solution. Just zero out the masked
        # values. It will do what we want in the FFT but the CASA mask
        # bookkeeping is being left in the dust. The workaround is
        # complicated, though, because you can't directly manipulate
        # pixel masks for some reason.

        if np.sum(combined_mask == False) > 0:
            myia.open(current_interf_file)
            interf_data = myia.getchunk()
            interf_data[combined_mask == False] = 0.0
            myia.putchunk(interf_data)
            myia.close()

            myia.open(current_sd_file)
            sd_data = myia.getchunk()
            sd_data[combined_mask == False] = 0.0
            myia.putchunk(sd_data)
            myia.close()

    else:
        current_interf_file = interf_file
        current_sd_file = sd_file

    # If apodization is requested, multiply both data sets by the same
    # taper and create a new, temporary output data set.

    if do_apodize:

        casaStuff.impbcor(imagename=current_sd_file,
                     pbimage=apod_file, 
                     outfile=current_sd_file+'.temp', 
                     mode='multiply')
        current_sd_file = current_sd_file+'.temp'
        
        casaStuff.impbcor(imagename=current_interf_file,
                     pbimage=apod_file, 
                     outfile=current_interf_file+'.temp', 
                     mode='multiply')
        current_interf_file = current_interf_file+'.temp'

    # Call feather, followed by an imsubimage to deal with degenerate
    # axis stuff.

    if overwrite:        
        os.system('rm -rf '+out_file)
    os.system('rm -rf '+out_file+'.temp')
    casaStuff.feather(imagename=out_file+'.temp',
                 highres=current_interf_file,
                 lowres=current_sd_file,
                 sdfactor=1.0,
                 lowpassfiltersd=False)
    casaStuff.imsubimage(imagename=out_file+'.temp', 
                    outfile=out_file,
                    dropdeg=True)
    os.system('rm -rf '+out_file+'.temp')

    # If we apodized, now divide out the common kernel.

    if do_apodize:
        os.system('rm -rf '+out_file+'.temp')
        os.system('mv '+out_file+' '+out_file+'.temp')
        casaStuff.impbcor(imagename=out_file+'.temp',
                     pbimage=apod_file, 
                     outfile=out_file, 
                     mode='divide', 
                     cutoff=apod_cutoff)

    # Remove temporary files

    os.system('rm -rf '+sd_file+'.temp')
    os.system('rm -rf '+interf_file+'.temp')
    os.system('rm -rf '+out_file+'.temp')

    os.system('rm -rf '+sd_file+'.temp.temp')
    os.system('rm -rf '+interf_file+'.temp.temp')
    os.system('rm -rf '+out_file+'.temp.temp')

    return(True)

#endregion
