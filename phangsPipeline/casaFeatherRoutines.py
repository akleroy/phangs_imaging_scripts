"""
Standalone routines that take input and output files and manipulate
cubes. These are called as part of the PHANGS post-processing pipeline
but also may be of general utility.
"""

#region Imports and definitions

import os
import numpy as np
import pyfits # CASA has pyfits, not astropy
import glob

# Analysis utilities
import analysisUtils as au

# CASA stuff
import casaStuff as casa

# Pipeline versionining
from pipelineVersion import version as pipeVer

#endregion

#region Feathering and single dish routines

def prep_sd_for_feather(
    sdfile_in=None,
    sdfile_out=None,
    interf_file=None,
    doimport=True,
    checkunits=True,
    doalign=True,
    overwrite=False,
    quiet=False):
    """
    Prepare single dish data for feathering. Import the data from
    FITS, check the units to make sure that they are in Jy/beam, and
    align the single dish data to the interferometric grid.
    """

    this_stub = 'PREP_SD_FOR_FEATHER: '

    # Check inputs
    
    if (os.path.isdir(interf_file) == False):
        if not quiet:
            print(this_stub+"Interferometric file not found: "+interf_file)
        return(None)

    if (os.path.isdir(sdfile_in) == False) and (os.path.isfile(sdfile_in) == False):
        if not quiet:
            print(this_stub+"Single dish file not found: "+sdfile_in)
        return(None)

    if sdfile_out is None:
        if not quiet:
            print(this_stub+"Output single dish file name not supplied via sdfile_out=")
        return(None)

    tempfile_name = sdfile_out+'.temp'
    if (os.path.isdir(tempfile_name)):
        if not overwrite:
            if not quiet:
                print(this_stub+"Temporary file is present but overwrite set to False - "+tempfile_name)
            return(None)
        else:
            os.system('rm -rf '+tempfile_name)

    current_file = sdfile_in    

    # Import from FITS if needed.
    
    if doimport:
        if (sdfile_in[-4:] == 'FITS') and os.path.isfile(sdfile_in):
            if not quiet:
                print(this_stub+"Importing from FITS.")
            if (os.path.isdir(tempfile_name)):
                if overwrite:
                    os.system('rm -rf '+tempfile_name)
                else:
                    if not quiet:
                        print(this_stub+"Temporary file is present but overwrite set to False - "+tempfile_name)
                    return(None)
                    
            importfits(fitsimage=sdfile_in, imagename=tempfile_name,
                       zeroblanks=True, overwrite=overwrite)
            current_file = tempfile_name

    # Check units on the singledish file.

    if checkunits:
        hdr = casa.imhead(current_file, mode='list')
        unit = hdr['bunit']
        if unit == 'K':
            if not quiet:
                print(this_stub+"Unit is Kelvin. Converting.")
            convert_ktojy(
                infile=current_file, 
                overwrite=overwrite, 
                inplace=True)

    # Align the single dish data to the interferometric data

    if doalign:
        casa.imregrid(
            imagename=current_file,
            template=interf_in,
            output=sdfile_out,
            asvelocity=True,
            axes=[-1],
            interpolation='cubic',
            overwrite=overwrite)

    if (os.path.isdir(tempfile_name)):
        if overwrite:
            os.system('rm -rf '+tempfile_name)

    return(None)

def feather_two_cubes(   
    interf_file=None,
    sd_file=None,
    out_file=None,
    apodize=False,
    apod_file=None,
    apod_cutoff=-1.0,
    blank=False,
    overwrite=False,
    quiet=False):
    """
    Feather the interferometric and total power data. Optionally,
    first apply some steps to homogenize the two data sets.
    """

    this_stub = 'FEATHER_TWO_CUBES: '

    if (os.path.isdir(sd_file) == False):
        if not quiet:
            print(this_stub+"Single dish file not found: "+sd_file)
        return(False)
        
    if (os.path.isdir(interf_file) == False):
        if not quiet:
            print(this_stub+"Interferometric file not found: "+interf_file)
        return(False)

    if apodize:
        if (os.path.isdir(apod_file) == False):
            if not quiet:
                print(this_stub+"Apodization requested, but file not found: "+apod_file)
            return(False)

    os.system('rm -rf '+sd_file+'.temp')
    os.system('rm -rf '+interf_file+'.temp')
    os.system('rm -rf '+out_file+'.temp')

    if blank:        
        current_interf_file = interf_file+'.temp'
        current_sd_file = sd_file+'.temp'

        os.system('cp -rf '+interf_file+' '+current_interf_file)
        os.system('cp -rf '+sd_file+' '+current_sd_file)

        myia = au.createCasaTool(iatool)

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
            myia.open(interf_file)
            interf_data = myia.getchunk()
            interf_data[combined_mask == False] = 0.0
            myia.putchunk(interf_data)
            myia.close()

            myia.open(sd_file)
            sd_data = myia.getchunk()
            sd_data[combined_mask == False] = 0.0
            myia.putchunk(sd_data)
            myia.close()

    else:
        current_interf_file = interf_file
        current_sd_file = sd_file

    # If apodization is requested, multiply both data sets by the same
    # taper and create a new, temporary output data set.

    if apodize:
        casa.impbcor(imagename=current_sd_file,
                     pbimage=apod_file, 
                     outfile=current_sd_file+'.temp', 
                     mode='multiply', 
                     cutoff=cutoff)
        current_sd_file = current_sd_file+'.temp'
        
        casa.impbcor(imagename=current_interf_file,
                     pbimage=apod_file, 
                     outfile=current_interf_file+'.temp', 
                     mode='multiply', 
                     cutoff=cutoff)
        current_interf_file = current_interf_file+'.temp'

    # Call feather, followed by an imsubimage to deal with degenerate
    # axis stuff.

    if overwrite:        
        os.system('rm -rf '+out_file)
    os.system('rm -rf '+out_file+'.temp')
    casa.feather(imagename=out_file+'.temp',
                 highres=current_interf_file,
                 lowres=current_sd_file,
                 sdfactor=1.0,
                 lowpassfiltersd=False)
    casa.imsubimage(imagename=out_file+'.temp', 
                    outfile=out_file,
                    dropdeg=True)
    os.system('rm -rf '+out_file+'.temp')

    # If we apodized, now divide out the common kernel.

    if apodize:
        os.system('rm -rf '+out_file+'.temp')
        os.system('mv '+out_file+' '+out_file+'.temp')
        casa.impbcor(imagename=out_file+'.temp',
                     pbimage=apod_file, 
                     outfile=out_file, 
                     mode='divide', 
                     cutoff=cutoff)

    # Remove temporary files

    os.system('rm -rf '+sd_file+'.temp')
    os.system('rm -rf '+interf_file+'.temp')
    os.system('rm -rf '+out_file+'.temp')

    os.system('rm -rf '+sd_file+'.temp.temp')
    os.system('rm -rf '+interf_file+'.temp.temp')
    os.system('rm -rf '+out_file+'.temp.temp')

    return(True)

#endregion
