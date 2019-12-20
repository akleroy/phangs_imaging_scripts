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

#region FITS file manipulation, scaling, extraction, and unit handling

def primary_beam_correct(
    infile=None, 
    pbfile=None, 
    outfile=None, 
    cutoff=0.25, 
    overwrite=False):
    """
    Construct a primary-beam corrected image.
    """

    if infile is None or pbfile is None or outfile is None:
        print("Missing required input.")
        return

    if os.path.isdir(infile) == False:
        print("Input file missing - "+infile)
        return

    if os.path.isdir(pbfile) == False:
        print("Primary beam file missing - "+pbfile)
        return

    if overwrite:
        os.system('rm -rf '+outfile)

    casa.impbcor(imagename=infile, pbimage=pbfile, outfile=outfile, cutoff=cutoff)

    return()
    
def convolve_to_round_beam(
    infile=None, 
    outfile=None, 
    force_beam=None, 
    overwrite=False):
    """
    Convolve supplied image to have a round beam. Optionally, force
    that beam to some size, else it figures out the beam.
    """
    
    if infile is None or outfile is None:
        print("Missing required input.")
        return

    if os.path.isdir(infile) == False:
        print("Input file missing - "+infile)
        return    

    if force_beam is None:
        hdr = casa.imhead(infile)

        if (hdr['axisunits'][0] != 'rad'):
            print("ERROR: Based on CASA experience. I expected units of radians.")
            print("I did not find this. Returning. Adjust code or investigate file "+infile)
            return
        pixel_as = abs(hdr['incr'][0]/np.pi*180.0*3600.)

        if (hdr['restoringbeam']['major']['unit'] != 'arcsec'):
            print("ERROR: Based on CASA experience. I expected units of arcseconds for the beam.")
            print("I did not find this. Returning. Adjust code or investigate file "+infile)
            return    
        bmaj = hdr['restoringbeam']['major']['value']    
        target_bmaj = np.sqrt((bmaj)**2+(2.0*pixel_as)**2)
    else:
        target_bmaj = force_beam

    casa.imsmooth(imagename=infile,
                  outfile=outfile,
                  targetres=True,
                  major=str(target_bmaj)+'arcsec',
                  minor=str(target_bmaj)+'arcsec',
                  pa='0.0deg',
                  overwrite=overwrite
                  )

    return(target_bmaj)

def calc_jytok(
    hdr=None,
    infile=None):
    """
    Calculate the Jy/beam -> Kelvin conversion. Accepts a header
    already read using imhead or a file name that it will parse.
    """

    c = 2.99792458e10
    h = 6.6260755e-27
    kb = 1.380658e-16

    if hdr is None:
        if infile is None:
            print("No header and no infile. Returning.")
            return(None)
        hdr = casa.imhead(target_file, mode='list')

    if hdr['cunit3'] != 'Hz':
        print("I expected frequency as the third axis but did not find it.")
        print("Returning.")
        return(None)
    
    crpix3 = hdr['crpix3']
    cdelt3 = hdr['cdelt3']
    crval3 = hdr['crval3']
    naxis3 = hdr['shape'][2]
    faxis_hz = (np.arange(naxis3)+1.-crpix3)*cdelt3+crval3
    freq_hz = np.mean(faxis_hz)
    
    bmaj_unit = hdr['beammajor']['unit']
    if bmaj_unit != 'arcsec':
        print("Beam unit is not arcsec, which I expected. Returning.")
        print("Unit instead is "+bmaj_unit)
        return    
    bmaj_as = hdr['beammajor']['value']
    bmin_as = hdr['beamminor']['value']
    bmaj_sr = bmaj_as/3600.*np.pi/180.
    bmin_sr = bmin_as/3600.*np.pi/180.
    beam_in_sr = np.pi*(bmaj_sr/2.0*bmin_sr/2.0)/np.log(2)
    
    jytok = c**2 / beam_in_sr / 1e23 / (2*kb*freq_hz**2)

    return(jytok)

def convert_jytok(
    infile=None, 
    outfile=None, 
    overwrite=False, 
    inplace=False):
    """
    Convert a cube from Jy/beam to K.
    """

    if infile is None or (outfile is None and inplace==False):
        print("Missing required input.")
        return
    
    if os.path.isdir(infile) == False:
        print("Input file not found: "+infile)
        return
    
    if inplace == False:
        if overwrite:
            os.system('rm -rf '+outfile)
        
        if os.path.isdir(outfile):
            print("Output file already present: "+outfile)
            return

        os.system('cp -r '+infile+' '+outfile)
        target_file = outfile
    else:
        target_file = infile

    hdr = casa.imhead(target_file, mode='list')
    unit = hdr['bunit']
    if unit != 'Jy/beam':
        print("Unit is not Jy/beam. Returning.")
        return
    
    jytok = calc_jytok(hdr=hdr)

    myia = au.createCasaTool(iatool)
    myia.open(target_file)
    vals = myia.getchunk()
    vals *= jytok
    myia.putchunk(vals)
    myia.setbrightnessunit("K")
    myia.close()

    casa.imhead(target_file, mode='put', hdkey='JYTOK', hdvalue=jytok)

    return()

def convert_ktojy(
    infile=None, 
    outfile=None, 
    overwrite=False, 
    inplace=False):
    """
    Convert a cube from K to Jy/beam.
    """

    if infile is None or (outfile is None and inplace==False):
        print("Missing required input.")
        return
    
    if os.path.isdir(infile) == False:
        print("Input file not found: "+infile)
        return
    
    if inplace == False:
        if overwrite:
            os.system('rm -rf '+outfile)
        
        if os.path.isdir(outfile):
            print("Output file already present: "+outfile)
            return

        os.system('cp -r '+infile+' '+outfile)
        target_file = outfile
    else:
        target_file = infile

    hdr = casa.imhead(target_file, mode='list')
    unit = hdr['bunit']
    if unit != 'K':
        print("Unit is not K. Returning.")
        return
    
    jytok = calc_jytok(hdr=hdr)

    myia = au.createCasaTool(iatool)
    myia.open(target_file)
    vals = myia.getchunk()
    vals *= 1.0/jytok
    myia.putchunk(vals)
    myia.setbrightnessunit("Jy/beam")
    myia.close()

    casa.imhead(target_file, mode='put', hdkey='JYTOK', hdvalue=jytok)

    return()

def trim_cube(    
    infile=None, 
    outfile=None, 
    overwrite=False, 
    inplace=False, 
    min_pixperbeam=3):
    """
    Trim empty space from around the edge of a cube. Also rebin the
    cube to smaller size, while ensuring a minimum number of pixels
    across the beam. Used to reduce the volume of cubes.
    """
    
    if infile is None or outfile is None:
        print("TRIM_CUBE: Missing required input.")
        return
    
    if os.path.isdir(infile) == False:
        print("TRIM_CUBE: Input file not found: "+infile)
        return

    # First, rebin if needed
    hdr = casa.imhead(infile)
    if (hdr['axisunits'][0] != 'rad'):
        print("ERROR: Based on CASA experience. I expected units of radians.")
        print("I did not find this. Returning. Adjust code or investigate file "+infile)
        return

    pixel_as = abs(hdr['incr'][0]/np.pi*180.0*3600.)

    if (hdr['restoringbeam']['major']['unit'] != 'arcsec'):
        print("ERROR: Based on CASA experience. I expected units of arcseconds for the beam.")
        print("I did not find this. Returning. Adjust code or investigate file "+infile)
        return    
    bmaj = hdr['restoringbeam']['major']['value']    
    
    pix_per_beam = bmaj*1.0 / pixel_as*1.0
    
    if pix_per_beam > 6:
        casa.imrebin(
            imagename=infile,
            outfile=outfile+'.temp',
            factor=[2,2,1],
            crop=True,
            dropdeg=True,
            overwrite=overwrite,
            )
    else:
        os.system('cp -r '+infile+' '+outfile+'.temp')

    # Figure out the extent of the image inside the cube
    myia = au.createCasaTool(iatool)
    myia.open(outfile+'.temp')
    mask = myia.getchunk(getmask=True)    
    myia.close()

    this_shape = mask.shape

    mask_spec_x = np.sum(np.sum(mask*1.0,axis=2),axis=1) > 0
    pad = 0
    xmin = np.max([0,np.min(np.where(mask_spec_x))-pad])
    xmax = np.min([np.max(np.where(mask_spec_x))+pad,mask.shape[0]-1])

    mask_spec_y = np.sum(np.sum(mask*1.0,axis=2),axis=0) > 0
    ymin = np.max([0,np.min(np.where(mask_spec_y))-pad])
    ymax = np.min([np.max(np.where(mask_spec_y))+pad,mask.shape[1]-1])

    mask_spec_z = np.sum(np.sum(mask*1.0,axis=0),axis=0) > 0
    zmin = np.max([0,np.min(np.where(mask_spec_z))-pad])
    zmax = np.min([np.max(np.where(mask_spec_z))+pad,mask.shape[2]-1])
    
    box_string = ''+str(xmin)+','+str(ymin)+','+str(xmax)+','+str(ymax)
    chan_string = ''+str(zmin)+'~'+str(zmax)

    print("... ... ... box selection: "+box_string)
    print("... ... ... channel selection: "+chan_string)

    if overwrite:
        os.system('rm -rf '+outfile)
        casa.imsubimage(
            imagename=outfile+'.temp',
            outfile=outfile,
            box=box_string,
            chans=chan_string,
            )
    
    os.system('rm -rf '+outfile+'.temp')

    return()
    
def export_and_cleanup(
    infile=None,
    outfile=None,
    overwrite=False,    
    remove_cards=[],
    add_cards=[],
    add_history=[],
    zap_history=True,
    roundbeam_tol=0.01):
    """
    Export from a CASA image file to a FITS file, in the process
    cleaning up header keywords that are usually confusing or
    useless. Optionally add new keywords to the header and check
    whether the beam is close enough to being round that it makes
    sense to overwrite it.
    """
    
    if infile is None or outfile is None:
        print("EXPORT_AND_CLEANUP: Missing required input.")
        return

    casa.exportfits(imagename=infile, 
                    fitsimage=outfile,
                    velocity=True, overwrite=True, dropstokes=True, 
                    dropdeg=True, bitpix=-32)
    
    # Clean up headers

    hdu = pyfits.open(outfile)

    hdr = hdu[0].header
    data = hdu[0].data
    
    # Cards to remove by default

    for card in ['BLANK','DATE-OBS','OBSERVER','O_BLANK','O_BSCALE',
                 'O_BZERO','OBSRA','OBSDEC','OBSGEO-X','OBSGEO-Y','OBSGEO-Z',
                 'DISTANCE']:
        if card in hdr.keys():
            hdr.remove(card)
            
    # User cards to remove

    for card in remove_cards:
        if card in hdr.keys():
            hdr.remove(card)
            
    # Delete history
    
    if zap_history:
        while 'HISTORY' in hdr.keys():
            hdr.remove('HISTORY')

    # Add history

    for history_line in add_history:
        add_history(history_line)
            
    # Get the data min and max right

    datamax = np.nanmax(data)
    datamin = np.nanmin(data)
    hdr['DATAMAX'] = datamax
    hdr['DATAMIN'] = datamin

    # Round the beam recorded in the header if it lies within the
    # specified tolerance.

    if roundbeam_tol > 0.0:
        bmaj = hdr['BMAJ']
        bmin = hdr['BMIN']
        if bmaj != bmin:
            frac_dev = np.abs(bmaj-bmin)/bmaj
            if frac_dev <= roundbeam_tol:
                print("Rounding beam.")
                hdr['BMAJ'] = bmaj
                hdr['BMIN'] = bmaj
                hdr['BPA'] = 0.0
            else:
                print("Beam too asymmetric to round.")
                print("... fractional deviation: "+str(frac_dev))
    
    # Overwrite

    hdu.writeto(outfile, clobber=True)
        
    return()

#endregion

#region Linear mosaicking routines

def common_res_for_mosaic(
    infile_list = None, 
    outfile_list = None,
    target_res=None,
    doconvolve=True,
    pixel_padding=2.0,
    overwrite=False):
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
        print("Missing required infile_list.")
        return(None)   
    
    for this_file in infile_list:
        if os.path.isdir(this_file) == False:
            print("File not found "+this_file)
            return(None)
    
    # Figure out target resolution if it is not supplied by the user

    if target_res is None:
        print("Calculating target resolution ... ")

        bmaj_list = []
        pix_list = []

        for infile in infile_list:
            print("Checking "+infile)

            hdr = casa.imhead(infile)

            if (hdr['axisunits'][0] != 'rad'):
                print("ERROR: Based on CASA experience. I expected units of radians.")
                print("I did not find this. Returning. Adjust code or investigate file "+infile)
                return(None)
            this_pixel = abs(hdr['incr'][0]/np.pi*180.0*3600.)

            if (hdr['restoringbeam']['major']['unit'] != 'arcsec'):
                print("ERROR: Based on CASA experience. I expected units of arcseconds for the beam.")
                print("I did not find this. Returning. Adjust code or investigate file "+infile)
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
        print("Missing outfile_list required for convolution.")
        return(None)

    if (type(outfile_list) != type([])) and (type(outfile_list) != type({})):
        print("outfile_list must be dictionary or list.")
        return(None)

    if type(outfile_list) == type([]):
        if len(infile_list) != len(outfile_list):
            print("Mismatch in input and output list lengths.")
            return
        outfile_dict = {}
        for ii in range(len(infile_list)):
            outfile_dict[infile_list[ii]] = outfile_list[ii]

    if type(outfile_list) == type({}):
        outfile_dict = outfile_list

    missing_keys = 0
    for infile in infile_list:
        if infile not in outfile_dict.keys():
            print("Missing output file for infile: "+infile)
            missing_keys += 1
    if missing_keys > 0:
        print("Missing "+str(missing_keys)+" output file names.")
        return(target_bmaj)

    # With a target resolution and matched lists we can proceed.

    for this_infile in infile_list:
        this_outfile = outfile_dict[this_infile]
        print("Convolving "+this_infile+' to '+this_outfile)
        
        casa.imsmooth(imagename=this_infile,
                      outfile=this_outfile,
                      targetres=True,
                      major=str(target_bmaj)+'arcsec',
                      minor=str(target_bmaj)+'arcsec',
                      pa='0.0deg',
                      overwrite=overwrite
                      )

    return(target_bmaj)

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
        print("Missing required infile_list.")
        return(None)

    for this_infile in infile_list:
        if not os.path.isdir(this_infile):
            print("File not found "+this_infile)
            print("Returning.")
            return(None)

    # The list of corner RA and Dec positions.
    ra_list = []
    dec_list = []

    # TBD
    freq_list = []

    for this_infile in infile_list:
        this_hdr = casa.imhead(this_infile)

        if this_hdr['axisnames'][0] != 'Right Ascension':
            print("Expected axis 0 to be Right Ascension. Returning.")
            return(None)
        if this_hdr['axisunits'][0] != 'rad':
            print("Expected axis units to be radians. Returning.")
            return(None)
        if this_hdr['axisnames'][1] != 'Declination':
            print("Expected axis 1 to be Declination. Returning.")
            return(None)
        if this_hdr['axisunits'][1] != 'rad':
            print("Expected axis units to be radians. Returning.")
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

    # TBD
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
    toobig=1e4):
    """
    Build a target header to be used as a template by imregrid. RA_CTR
    and DEC_CTR are assumed to be in decimal degrees. DELTA_RA and
    DELTA_DEC are assumed to be in arcseconds.
    """
    
    if infile_list is None:
        print("Missing required infile_list.")
        return(None)

    for this_infile in infile_list:
        if not os.path.isdir(this_infile):
            print("File not found "+this_infile)
            print("Returning.")
            return(None)

    if template_file is not None:
        if os.path.isdir(template_file) == False:
            print("The specified template file does not exist.")
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
        print("ERROR: Based on CASA experience. I expected pixel units of radians.")
        print("I did not find this. Returning. Adjust code or investigate file "+infile_list[0])
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
            print("WARNING! This is a very big image you plan to create, ", ra_axis_size, " x ", dec_axis_size)
            print("To make an image this big set allowbigimage=True. Returning.")
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
    overwrite=False):
    """
    Align a list of files to a target coordinate system.
    """

    if infile_list is None or \
            outfile_list is None or \
            target_hdr is None:
        print("Missing required input.")
        return(None)

    for ii in range(len(infile_list)):
        this_infile = infile_list[ii]
        this_outfile = outfile_list[ii]        

        if os.path.isdir(this_infile) == False:
            print("File "+this_infile+" not found. Continuing.")
            continue

        casa.imregrid(imagename=this_infile,
                      template=target_hdr,
                      output=this_outfile,
                      asvelocity=True,
                      axes=[-1],
                      interpolation='cubic',
                      overwrite=overwrite)

    return(None)

def mosaic_aligned_data(
    infile_list = None, 
    weightfile_list = None,
    outfile = None, 
    overwrite=False):
    """
    Combine a list of aligned data with primary-beam (i.e., inverse
    noise) weights using simple linear mosaicking.
    """

    if infile_list is None or \
            weightfile_list is None or \
            outfile is None:
        print("Missing required input.")
        return(None)

    sum_file = outfile+'.sum'
    weight_file = outfile+'.weight'

    if (os.path.isdir(outfile) or \
            os.path.isdir(sum_file) or \
            os.path.isdir(weight_file)) and \
            (overwrite == False):
        print("Output file present and overwrite off.")
        print("Returning.")
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
    
    myia = au.createCasaTool(iatool)
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

#region Feathering and single dish routines

def prep_sd_for_feather(
    interf_file=None,
    sdfile_in=None,
    sdfile_out=None,
    doimport=True,
    checkunits=True,
    doalign=True,
    overwrite=False):
    """
    Prepare single dish data for feathering. Import the data from
    FITS, check the units to make sure that they are in Jy/beam, and
    align the single dish data to the interferometric grid.
    """

    # Check inputs
    
    if (os.path.isdir(interf_file) == False):
        print("Interferometric file not found: "+interf_file)
        return(None)

    if (os.path.isdir(sdfile_in) == False) and (os.path.isfile(sdfile_in) == False):
        print("Single dish file not found: "+sdfile_in)
        return(None)

    if sdfile_out is None:
        print("Output single dish file name not supplied via sdfile_out=")
        return(None)

    if (os.path.isdir(sdfile_out+'.temp')):
        if not overwrite:
            print("Temporary file is present but overwrite set to False. Returning.")
            return(None)
        else:
            os.system('rm -rf '+sdfile_out+'.temp')

    current_file = sdfile_in    

    # Import from FITS if needed.
    
    if doimport:
        if (sdfile_in[-4:] == 'FITS') and os.path.isfile(sdfile_in):
            print("Importing from FITS.")
            if overwrite:
                os.system('rm -rf '+sdfile_out+'.temp')
            importfits(fitsimage=sdfile_in, imagename=sdfile_out+'.temp',
                       zeroblanks=True, overwrite=overwrite)
            current_file = sdfile_out+'.temp'

    # Check units on the singledish file.

    if checkunits:
        hdr = casa.imhead(current_file, mode='list')
        unit = hdr['bunit']
        if unit == 'K':
            print("Unit is Kelvin. Converting.")
            convert_ktojy(infile=current_file, 
                          overwrite=overwrite, inplace=True)

    # Align the single dish data to the interferometric data

    if doalign:
        casa.imregrid(imagename=current_file,
                      template=interf_in,
                      output=sdfile_out,
                      asvelocity=True,
                      axes=[-1],
                      interpolation='cubic',
                      overwrite=overwrite)

    if overwrite:
        os.system('rm -rf '+sdfile_out+'.temp')

    return(None)

def feather_two_cubes(   
    interf_file=None,
    sd_file=None,
    out_file=None,
    apodize=False,
    apod_file=None,
    apod_cutoff=-1.0,
    blank=False,
    overwrite=False):
    """
    Feather the interferometric and total power data. Optionally,
    first apply some steps to homogenize the two data sets.
    """

    if (os.path.isdir(sd_file) == False):
        print("Single dish file not found: "+sd_file)
        return
        
    if (os.path.isdir(interf_file) == False):
        print("Interferometric file not found: "+interf_file)
        return

    if apodize:
        if (os.path.isdir(apod_file) == False):
            print("Apodization requested, but file not found: "+apod_file)
            return

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

#endregion
