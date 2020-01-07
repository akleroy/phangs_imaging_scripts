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

import logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

# Analysis utilities
import analysisUtils as au

# CASA stuff
import casaStuff as casa

# Other pipeline stuff
import casaMaskingRoutines as cma

# Pipeline versionining
from pipelineVersion import version as pipeVer

#endregion

#region Routines to match resolution

def common_res_for_mosaic(
    infile_list = None, 
    outfile_list = None,
    target_res=None,
    doconvolve=True,
    pixel_padding=2.0,
    overwrite=False
    ):
    """
    Convolve multi-part cubes to a common res for mosaicking. It will
    calculate the common resolution based on the beam size of all of
    the input files, unless it is fed a fixed target resolution. It
    returns the target resolution. If doconvolve is True, it also
    convolves all of the input files to output files with that
    resolution. For this it needs a list of output files matched to
    the input file list, either as another list or a dictionary.

    Has a tuning keyword 'pixel_padding' to indicate how many pixels
    worth of padding is added in quadrature to the maximum beam size
    to get the common beam.
    """
    
    # Check that the input files exist

    if infile_list is None:
        logger.error("Missing required infile_list.")
        return(None)   
    
    for this_file in infile_list:
        if os.path.isdir(this_file) == False:
            logger.error("File not found "+this_file)
            return(None)
    
    # Figure out target resolution if it is not supplied by the user

    if target_res is None:
        logger.debug("Calculating target resolution ... ")

        bmaj_list = []
        pix_list = []

        for infile in infile_list:
            logger.info("Checking "+infile)

            hdr = casa.imhead(infile)

            if (hdr['axisunits'][0] != 'rad'):
                logger.error("ERROR: Based on CASA experience. I expected units of radians.")
                logger.error("I did not find this. Returning. Adjust code or investigate file "+infile)
                return(None)
            this_pixel = abs(hdr['incr'][0]/np.pi*180.0*3600.)

            if (hdr['restoringbeam']['major']['unit'] != 'arcsec'):
                logger.error("ERROR: Based on CASA experience. I expected units of arcseconds for the beam.")
                logger.error("I did not find this. Returning. Adjust code or investigate file "+infile)
                return(None)
            this_bmaj = hdr['restoringbeam']['major']['value']

            bmaj_list.append(this_bmaj)
            pix_list.append(this_pixel)
        
        max_bmaj = np.max(bmaj_list)
        max_pix = np.max(pix_list)
        target_bmaj = np.sqrt((max_bmaj)**2+(pixel_padding*max_pix)**2)
    else:
        target_bmaj = force_beam

    if not doconvolve:
        return(target_bmaj)

    # If doconvolve is True then make sure that we have output files
    # and that they match the input files.

    if outfile_list is None:
        logger.error("Missing outfile_list required for convolution.")
        return(target_bmaj)

    if (type(outfile_list) != type([])) and (type(outfile_list) != type({})):
        logger.error("outfile_list must be dictionary or list.")
        return(target_bmaj)

    if type(outfile_list) == type([]):
        if len(infile_list) != len(outfile_list):
            logger.error("Mismatch in input and output list lengths.")
            return(target_bmaj)
        outfile_dict = {}
        for ii in range(len(infile_list)):
            outfile_dict[infile_list[ii]] = outfile_list[ii]

    if type(outfile_list) == type({}):
        outfile_dict = outfile_list

    missing_keys = 0
    for infile in infile_list:
        if infile not in outfile_dict.keys():
            logger.error("Missing output file for infile: "+infile)
            missing_keys += 1
    if missing_keys > 0:
        logger.error("Missing "+str(missing_keys)+" output file names.")
        return(target_bmaj)

    # With a target resolution and matched lists we can proceed.

    for this_infile in infile_list:
        this_outfile = outfile_dict[this_infile]
        logger.debug("Convolving "+this_infile+' to '+this_outfile)
        
        casa.imsmooth(imagename=this_infile,
                      outfile=this_outfile,
                      targetres=True,
                      major=str(target_bmaj)+'arcsec',
                      minor=str(target_bmaj)+'arcsec',
                      pa='0.0deg',
                      overwrite=overwrite
                      )

    return(target_bmaj)

#endregion

#region Routines to match astrometry

def calculate_mosaic_extent(
    infile_list = None, 
    force_ra_ctr = None, 
    force_dec_ctr = None,
    ):
    """
    Given a list of input files, calculate the center and extent of
    the mosaic needed to cover them all. Optionally, force the center
    of the mosaic to some value. In that case, the extent is
    calculated to be the rectangle centered on the supplied (RA,
    Dec). 

    If the RA and Dec. center are supplied they are assumed to be in
    decimal degrees.
    """

    if infile_list is None:
        logger.error("Missing required infile_list.")
        return(None)

    for this_infile in infile_list:
        if not os.path.isdir(this_infile):
            logger.error("File not found "+this_infile+"Returning.")
            return(None)

    # The list of corner RA and Dec positions.
    ra_list = []
    dec_list = []

    # TBD - right now we assume matched frequency/velocity axis
    freq_list = []

    for this_infile in infile_list:
        this_hdr = casa.imhead(this_infile)

        if this_hdr['axisnames'][0] != 'Right Ascension':
            logger.error("Expected axis 0 to be Right Ascension. Returning.")
            return(None)
        if this_hdr['axisunits'][0] != 'rad':
            logger.error("Expected axis units to be radians. Returning.")
            return(None)
        if this_hdr['axisnames'][1] != 'Declination':
            logger.error("Expected axis 1 to be Declination. Returning.")
            return(None)
        if this_hdr['axisunits'][1] != 'rad':
            logger.error("Expected axis units to be radians. Returning.")
            return(None)

        this_shape = this_hdr['shape']
        xlo = 0
        xhi = this_shape[0]-1
        ylo = 0
        yhi = this_shape[1]-1

        pixbox = str(xlo)+','+str(ylo)+','+str(xlo)+','+str(ylo)
        blc = imval(this_infile, chans='0', box=pixbox)

        pixbox = str(xlo)+','+str(yhi)+','+str(xlo)+','+str(yhi)
        tlc = imval(this_infile, chans='0', box=pixbox)
        
        pixbox = str(xhi)+','+str(yhi)+','+str(xhi)+','+str(yhi)
        trc = imval(this_infile, chans='0', box=pixbox)

        pixbox = str(xhi)+','+str(ylo)+','+str(xhi)+','+str(ylo)
        brc = imval(this_infile, chans='0', box=pixbox)
        
        ra_list.append(blc['coords'][0][0])
        ra_list.append(tlc['coords'][0][0])
        ra_list.append(trc['coords'][0][0])
        ra_list.append(brc['coords'][0][0])

        dec_list.append(blc['coords'][0][1])
        dec_list.append(tlc['coords'][0][1])
        dec_list.append(trc['coords'][0][1])
        dec_list.append(brc['coords'][0][1])
        
    min_ra = np.min(ra_list)
    max_ra = np.max(ra_list)
    min_dec = np.min(dec_list)
    max_dec = np.max(dec_list)

    # TBD - right now we assume matched frequency/velocity axis
    min_freq = None
    max_freq = None

    if force_ra_ctr == None:
        ra_ctr = (max_ra+min_ra)*0.5
    else:
        ra_ctr = force_ra_ctr*np.pi/180.

    if force_dec_ctr == None:
        dec_ctr = (max_dec+min_dec)*0.5
    else:
        dec_ctr = force_dec_ctr*np.pi/180.

    delta_ra = 2.0*np.max([np.abs(max_ra-ra_ctr),np.abs(min_ra-ra_ctr)])
    delta_ra *= np.cos(dec_ctr)
    delta_dec = 2.0*np.max([np.abs(max_dec-dec_ctr),np.abs(min_dec-dec_ctr)])
    
    output = {
        'ra_ctr':[ra_ctr*180./np.pi,'degrees'],
        'dec_ctr':[dec_ctr*180./np.pi,'degrees'],
        'delta_ra':[delta_ra*180./np.pi*3600.,'arcsec'],
        'delta_dec':[delta_dec*180./np.pi*3600.,'arcsec'],
        }

    return(output)

def build_common_header(
    infile_list = None, 
    ra_ctr = None, 
    dec_ctr = None,
    delta_ra = None, 
    delta_dec = None,
    template_file = None,
    allowbigimage = False,
    toobig=1e4
    ):
    """
    Build a target header to be used as a template by imregrid. RA_CTR
    and DEC_CTR are assumed to be in decimal degrees. DELTA_RA and
    DELTA_DEC are assumed to be in arcseconds.
    """
    
    if infile_list is None:
        logger.error("Missing required infile_list.")
        return(None)

    for this_infile in infile_list:
        if not os.path.isdir(this_infile):
            logger.error("File not found "+this_infile+" . Returning.")
            return(None)

    if template_file is not None:
        if os.path.isdir(template_file) == False:
            logger.error("The specified template file does not exist.")
            return(None)

    # Base the target header on a template. If no template is supplied
    # than the template is the first file in the list.

    if template_file is None:
        template_file = infile_list[0]

    target_hdr = casa.imregrid(template_file, template='get')
    
    # Get the pixel scale. This makes some assumptions. We could put a
    # lot of general logic here, but we are usually working in a
    # pretty specific case.

    if (target_hdr['csys']['direction0']['units'][0] != 'rad') or \
            (target_hdr['csys']['direction0']['units'][1] != 'rad'):
        logger.error("ERROR: Based on CASA experience. I expected pixel units of radians.")
        logger.error("I did not find this. Returning. Adjust code or investigate file "+infile_list[0])
        return(None)

    # Add our target center pixel values to the header after
    # converting to radians.

    ra_ctr_in_rad = ra_ctr * np.pi / 180.
    dec_ctr_in_rad = dec_ctr * np.pi / 180.

    target_hdr['csys']['direction0']['crval'][0] = ra_ctr_in_rad
    target_hdr['csys']['direction0']['crval'][1] = dec_ctr_in_rad

    # Calculate the size of the image in pixels and set the central
    # pixel coordinate for the RA and Dec axis.
    
    ra_pix_in_as = np.abs(target_hdr['csys']['direction0']['cdelt'][0]*180./np.pi*3600.)
    ra_axis_size = np.ceil(delta_ra / ra_pix_in_as)
    new_ra_ctr_pix = ra_axis_size/2.0

    dec_pix_in_as = np.abs(target_hdr['csys']['direction0']['cdelt'][1]*180./np.pi*3600.)
    dec_axis_size = np.ceil(delta_dec / dec_pix_in_as)
    new_dec_ctr_pix = dec_axis_size/2.0
    
    # Check that the axis size isn't too big. This is likely to be a
    # bug. If allowbigimage is True then bypass this, otherwise exit.

    if not allowbigimage:
        if ra_axis_size > toobig or dec_axis_size > toobig:
            logger.error("WARNING! This is a very big image you plan to create, "+str(ra_axis_size)+ \
                             " x "+str(dec_axis_size))
            logger.error(" To make an image this big set allowbigimage=True. Returning.")
            return(None)

    # Enter the new values into the header and return.

    target_hdr['csys']['direction0']['crpix'][0] = new_ra_ctr_pix
    target_hdr['csys']['direction0']['crpix'][1] = new_dec_ctr_pix
    
    target_hdr['shap'][0] = int(ra_axis_size)
    target_hdr['shap'][1] = int(dec_axis_size)
    
    return(target_hdr)

def align_for_mosaic(
    infile_list = None,
    outfile_list = None,
    target_hdr=None,
    overwrite=False
    ):
    """
    Align a list of files to a target coordinate system.
    """

    if infile_list is None or \
            outfile_list is None or \
            target_hdr is None:
        logger.error("Missing required input.")
        return(None)

    for ii in range(len(infile_list)):
        this_infile = infile_list[ii]
        this_outfile = outfile_list[ii]        

        if os.path.isdir(this_infile) == False:
            logger.error("File "+this_infile+" not found. Continuing.")
            continue

        casa.imregrid(imagename=this_infile,
                      template=target_hdr,
                      output=this_outfile,
                      asvelocity=True,
                      axes=[-1],
                      interpolation='cubic',
                      overwrite=overwrite)

    return(None)

def common_astrometry_for_mosaic(
    infile_list = None,
    outfile_list = None,
    weightfiles_in = None,
    weightfiles_out = None,
    ra_ctr = None, 
    dec_ctr = None,
    delta_ra = None, 
    delta_dec = None,
    allowbigimage = False,
    toobig=1e4,
    overwrite=False,
    ):
    """
    Build a common astrometry for a mosaic and align all image and
    weight files to that astrometry. This wraps the other routines
    here and so can be called for an end-to-end alignment.
    """

    # Error checking - mostly the subprograms do this.

    if infile_list is None:
        logger.error("Infile list missing.")
        return(None)

    if outfile_list is None:
        logger.error("Outfile list missing.")
        return(None)

    # Determine common header using only the input images

    logger.info('Generating extent of target header.')
    extent = calculate_mosaic_extent(
        infile_list = infile_list, 
        force_ra_ctr = ra_ctr, 
        force_dec_ctr = dec_ctr)
    if ra_ctr is None:
        ra_ctr = extent['ra_ctr'][0]
    if dec_ctr is None:
        dec_ctr = extent['dec_ctr'][0]
    if delta_ra is None:
        delta_ra = extent['delta_ra'][0]
    if delta_dec is None:
        delta_dec = extent['delta_dec'][0]

    logger.info('Generating target header.')
    target_hdr = build_common_header(
        infile_list = infile_list, 
        ra_ctr = ra_ctr, 
        dec_ctr = dec_ctr,
        delta_ra = delta_ra, 
        delta_dec = delta_dec,
        allowbigimage = allowbigimage,
        toobig=toobig,
        )

    # Align the input files to the new astrometry

    logger.info('Aligning image files.')
    align_for_mosaic(
        infile_list = infile_list,
        outfile_list = outfile_list,
        target_hdr=target_hdr,
        overwrite=overwrite
        )

    # Align the weight files to the new astrometry

    if weightfiles_in is not None:
        logger.info('Aligning weighting files.')
        align_for_mosaic(
            infile_list = weightfiles_in,
            outfile_list = weightfiles_out,
            target_hdr=target_hdr,
            overwrite=overwrite
            )

    return(None)

#endregion

#region Routines to deal with weighting

def generate_weight_file(
    image_file = None,
    input_file = None,
    input_value = None,    
    input_type = 'pb',
    outfile = None,
    scale_by_noise = False,
    noise_value = None,
    scale_by_factor = None,
    overwrite=False,
    ):
    """
    Generate a weight image for use in a linear mosaic.
    """

    # Check input

    if image_file is None and input_file is None:
        logger.error("I need either an input or an image template file.")
        return(None)

    if input_file is None and input_value is None:
        logger.error("I need either an input value or an input file.")
        return(None)

    if input_file is not None and input_value is not None:
        logger.error("I need ONE OF an input value or an input file. Got both.")
        return(None)

    if outfile is None:
        logger.error("Specify output file.")
        return(None)

    if input_file is not None:
        valid_types = ['pb', 'noise', 'weight']
        if input_type not in valid_types:
            logger.error("Valid input types are :"+str(valid_types))
            return(None)

    if input_file is None and input_value is None:
        logger.error("Need either an input value or an input file.")
        return(None)

    if input_file is not None:
        if not os.path.isdir(input_file):
            logger.error("Missing input file directory - "+input_file)
            return(None)

    if image_file is not None:
        if not os.path.isdir(image_file):
            logger.error("Missing image file directory - "+image_file)
            return(None)

    # If scaling by noise is requested and no estimate is provided,
    # generate an estimate

    if scale_by_noise:

        if noise_value is None and image_file is None:
            logger.error("I can only scale by the noise if I get an image file to caluclate the noise. Returning.")
            return(None)

        if noise_value is None:
            
            pass

    # Define the template for the astrometry

    if input_file is None:
        template = image_file
    else:
        template = input_file

    # Check the output file

    if os.path.isdir(outfile) or os.path.isfile(outfile):
        if not overwrite:
            logger.error("File exists and overwrite set to false - "+outfile)
            return(None)
        os.system('rm -rf '+outfile)

    # Copy the template and read the data into memory

    os.system("cp -r "+template+" "+outfile)

    myia = au.createCasaTool(casa.iatool)
    myia.open(outfile)
    data = myia.getchunk()

    # Case 1 : We just have an input value.

    if input_file is None and input_value is not None:        

        if input_type is 'noise':
            weight_value = 1./input_value**2
        if input_type is 'pb':
            weight_value = input_value**2
        if input_type is 'weight':
            weight_value = input_value

        weight_image = data*0.0 + weight_value
    
    # Case 2 : We have an input image

    if input_file is not None:

        os.system("cp -r "+template+" "+outfile)

        if input_type is 'noise':
            weight_image = 1./data**2
        if input_type is 'pb':
            weight_image = data**2
        if input_type is 'weight':
            weight_image = data

    # If request, scale the data by a factor

    if scale_by_factor is not None:
        
        weight_image = weight_image * scale_by_factor

    # If request, scale the data by the inverse square of the noise estimate.

    if scale_by_noise:

        noise_scale_factor = 1./noise_value**2

    myia.putchunk(data)
    myia.close()

    return(None)
    

#endregion

#region Routines to carry out the mosaicking

def mosaic_aligned_data(
    infile_list = None, 
    weightfile_list = None,
    outfile = None, 
    overwrite=False
    ):
    """
    Combine a list of aligned data with primary-beam (i.e., inverse
    noise) weights using simple linear mosaicking.
    """

    if infile_list is None or \
            weightfile_list is None or \
            outfile is None:
        logger.error("Missing required input.")
        return(None)

    sum_file = outfile+'.sum'
    weight_file = outfile+'.weight'

    if (os.path.isdir(outfile) or \
            os.path.isdir(sum_file) or \
            os.path.isdir(weight_file)) and \
            (overwrite == False):
        logger.error("Output file present and overwrite off. Returning.")
        return(None)

    if overwrite:
        os.system('rm -rf '+outfile+'.temp')
        os.system('rm -rf '+outfile)
        os.system('rm -rf '+sum_file)
        os.system('rm -rf '+weight_file)
        os.system('rm -rf '+outfile+'.mask')

    imlist = infile_list[:]
    imlist.extend(weightfile_list)
    n_image = len(infile_list)
    lel_exp_sum = ''
    lel_exp_weight = ''
    first = True
    for ii in range(n_image):
        this_im = 'IM'+str(ii)
        this_wt = 'IM'+str(ii+n_image)
        this_lel_sum = '('+this_im+'*'+this_wt+'*'+this_wt+')'
        this_lel_weight = '('+this_wt+'*'+this_wt+')'
        if first:
            lel_exp_sum += this_lel_sum
            lel_exp_weight += this_lel_weight
            first=False
        else:
            lel_exp_sum += '+'+this_lel_sum
            lel_exp_weight += '+'+this_lel_weight

    casa.immath(imagename = imlist, mode='evalexpr',
                expr=lel_exp_sum, outfile=sum_file,
                imagemd = imlist[0])
    
    myia = au.createCasaTool(casa.iatool)
    myia.open(sum_file)
    myia.set(pixelmask=1)
    myia.close()

    casa.immath(imagename = imlist, mode='evalexpr',
                expr=lel_exp_weight, outfile=weight_file)
    myia.open(weight_file)
    myia.set(pixelmask=1)
    myia.close()

    casa.immath(imagename = [sum_file, weight_file], mode='evalexpr',
                expr='iif(IM1 > 0.0, IM0/IM1, 0.0)', outfile=outfile+'.temp',
                imagemd = sum_file)

    casa.immath(imagename = weight_file, mode='evalexpr',
                expr='iif(IM0 > 0.0, 1.0, 0.0)', outfile=outfile+'.mask')
    
    casa.imsubimage(imagename=outfile+'.temp', outfile=outfile,
                    mask='"'+outfile+'.mask"', dropdeg=True)

    return(None)

#endregion
