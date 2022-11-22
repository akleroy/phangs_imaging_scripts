"""
Standalone routines that take input and output files and manipulate
cubes. These are called as part of the PHANGS post-processing pipeline
but also may be of general utility.
"""

#region Imports and definitions

import os
import glob
import logging

import numpy as np
import scipy.ndimage as nd
try:
    import pyfits  # CASA has pyfits, not astropy
except ImportError:
    import astropy.io.fits as pyfits

# Analysis utilities
import analysisUtils as au

# Pipeline versioning
from .pipelineVersion import version as pipeVer

# CASA stuff
from . import casaStuff

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

#endregion

#region Copying, scaling, etc.

def copy_dropdeg(
    infile=None,
    outfile=None,
    overwrite=False
    ):
    """
    Copy using imsubimage to drop degenerate axes. Optionally handle
    overwriting and importing from FITS.
    """

    if os.path.isdir(outfile):
        if not overwrite:
            logger.error("Output exists and overwrite set to false - "+outfile)
            return(False)
        os.system('rm -rf '+outfile)

    used_temp_outfile = False
    if (infile[-4:] == 'FITS') and os.path.isfile(infile):
        logger.info("Importing from FITS.")
        temp_outfile = outfile+'.temp'
        if os.path.isdir(temp_outfile):
            if not overwrite:
                logger.error("Temp file exists and overwrite set to false - "+temp_outfile)
                return(False)
            os.system('rm -rf '+temp_outfile)

        importfits(fitsimage=infile,
                   imagename=temp_outfile,
                   zeroblanks=False,
                   overwrite=overwrite)
        used_temp_outfile = True

    if used_temp_outfile:
        casaStuff.imsubimage(
            imagename=temp_outfile,
            outfile=outfile,
            dropdeg=True)
        os.system('rm -rf '+temp_outfile)
    else:
        casaStuff.imsubimage(
            imagename=infile,
            outfile=outfile,
            dropdeg=True)

    return(True)


def get_mask(infile, huge_cube_workaround=True):
    """
    Get a mask from a CASA image file. Includes a switch for large cubes, where getchunk can segfault.
    """

    if huge_cube_workaround:
        os.system('rm -rf ' + infile + '.temp_deg_ordered')
        casaStuff.imtrans(imagename=infile + '.temp_deg', outfile=infile + '.temp_deg_ordered',
                          order='0132')

        casaStuff.makemask(mode='copy', inpimage=infile + '.temp_deg_ordered',
                           inpmask=infile + '.temp_deg_ordered:mask0',
                           output=infile + '.temp_mask', overwrite=True)

        casaStuff.exportfits(imagename=infile + '.temp_mask',
                             fitsimage=infile + '.temp.fits',
                             stokeslast=False, overwrite=True)
        hdu = pyfits.open(infile + '.temp.fits')[0]
        mask = hdu.data.T[:, :, 0, :]

        os.system('rm -rf ' + infile + '.temp_deg_ordered')
        os.system('rm -rf ' + infile + '.temp_mask')
        os.system('rm -rf ' + infile + '.temp.fits')

    else:
        myia = au.createCasaTool(casaStuff.iatool)
        myia.open(infile+'.temp_deg')
        mask = myia.getchunk(getmask=True)
        myia.close()

    os.system('rm -rf ' + infile + '.temp_deg')

    return mask


def copy_mask(infile, outfile, huge_cube_workaround=True):
    """
    Copy a mask from infile to outfile. Includes a switch for large cubes, where getchunk/putchunk can segfault
    """

    if huge_cube_workaround:
        os.system('rm -rf ' + outfile + '/mask0')
        os.system('cp -r ' + infile + '/mask0' + ' ' + outfile + '/mask0')
    else:
        myia = au.createCasaTool(casaStuff.iatool)
        myia.open(infile)
        mask = myia.getchunk(getmask=True)
        myia.close()

        myia = au.createCasaTool(casaStuff.iatool)
        myia.open(outfile)
        mask = myia.putregion(pixelmask=mask)
        myia.close()

    return True


def multiply_cube_by_value(infile, value, brightness_unit, huge_cube_workaround=True):
    """
    Multiply a cube by some value, and update the brightness unit accordingly. Includes a switch for large cubes, where
    getchunk/putchunk may fail.
    """

    if huge_cube_workaround:
        casaStuff.exportfits(imagename=infile,
                             fitsimage=infile + '.fits',
                             overwrite=True)
        hdu = pyfits.open(infile + '.fits')[0]
        hdu.data *= value

        hdu.writeto(infile + '.fits', clobber=True)
        casaStuff.importfits(fitsimage=infile + '.fits',
                             imagename=infile,
                             overwrite=True)
        os.system('rm -rf ' + infile + '.fits')
    else:
        myia = au.createCasaTool(casaStuff.iatool)
        myia.open(infile)
        vals = myia.getchunk()
        vals *= value
        myia.putchunk(vals)
        myia.close()

    myia = au.createCasaTool(casaStuff.iatool)
    myia.open(infile)
    myia.setbrightnessunit(brightness_unit)
    myia.close()

    return True


def export_and_cleanup(
    infile=None,
    outfile=None,
    overwrite=False,
    remove_cards=[],
    add_cards={},
    add_history=[],
    zap_history=True,
    round_beam=True,
    roundbeam_tol=0.01
    ):
    """
    Export from a CASA image file to a FITS file, in the process
    cleaning up header keywords that are usually confusing or
    useless. Optionally add new keywords to the header and check
    whether the beam is close enough to being round that it makes
    sense to overwrite it.
    """

    if infile is None or outfile is None:
        logger.error("Missing required input.")
        return(False)

    if os.path.isdir(infile) == False:
        logger.error("Input file does not exist - "+infile)
        return(False)

    if os.path.isfile(outfile):
        if not overwrite:
            logger.error("Output exists and overwrite set to false - "+outfile)
            return(False)

    casaStuff.exportfits(imagename=infile,
                    fitsimage=outfile,
                    velocity=True,
                    overwrite=True,
                    dropstokes=True,
                    dropdeg=True,
                    bitpix=-32)

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

    for card in add_cards:
        hdr[card] = add_cards[card]

    # Get the data min and max right

    datamax = np.nanmax(data)
    datamin = np.nanmin(data)
    hdr['DATAMAX'] = datamax
    hdr['DATAMIN'] = datamin


    # Round the beam recorded in the header if it lies within the
    # specified tolerance.

    if round_beam:
        if roundbeam_tol > 0.0:
            bmaj = hdr['BMAJ']
            bmin = hdr['BMIN']
            if bmaj != bmin:
                frac_dev = np.abs(bmaj-bmin)/bmaj
                if frac_dev <= roundbeam_tol:
                    logger.info("Rounding beam.")
                    hdr['BMAJ'] = bmaj
                    hdr['BMIN'] = bmaj
                    hdr['BPA'] = 0.0
                else:
                    logger.info("Beam too asymmetric to round.")
                    logger.info("... fractional deviation: "+str(frac_dev))

    # Never forget where you came from
    hdr['COMMENT'] = 'Produced with PHANGS-ALMA pipeline version ' + pipeVer

    # Overwrite
    hdu.writeto(outfile, clobber=True)

    return()

def trim_cube(
        infile=None,
        outfile=None,
        overwrite=False,
        inplace=False,
        min_pixperbeam=3,
        pad=1
    ):
    """
    Trim empty space from around the edge of a cube. Also rebin the
    cube to smaller size, while ensuring a minimum number of pixels
    across the beam. Used to reduce the volume of cubes.
    """

    if infile is None or outfile is None:
        logger.error("Missing required input.")
        return(False)

    if os.path.isdir(infile) == False:
        logger.error("Input file not found: "+infile)
        return(False)

    # First, rebin if needed
    hdr = casaStuff.imhead(infile)
    if (hdr['axisunits'][0] != 'rad'):
        logger.error("ERROR: Based on CASA experience. I expected units of radians. I did not find this. Returning.")
        logger.error("Adjust code or investigate file "+infile)
        return(False)

    pixel_as = abs(hdr['incr'][0]/np.pi*180.0*3600.)

    if (hdr['restoringbeam']['major']['unit'] != 'arcsec'):
        logger.error("ERROR: Based on CASA experience. I expected units of arcseconds for the beam. I did not find this. Returning.")
        logger.error("Adjust code or investigate file "+infile)
        return(False)
    bmaj = hdr['restoringbeam']['major']['value']

    pix_per_beam = bmaj*1.0 / pixel_as*1.0

    if pix_per_beam > 6:
        casaStuff.imrebin(
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

    myia = au.createCasaTool(casaStuff.iatool)
    myia.open(outfile + '.temp')
    myia.adddegaxes(outfile=outfile + '.temp_deg', stokes='I', overwrite=True)
    myia.close()

    mask = get_mask(outfile, huge_cube_workaround=True)

    this_shape = mask.shape

    mask_spec_x = np.sum(np.sum(mask*1.0,axis=2),axis=1) > 0
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

    logger.info("... box selection: "+box_string)
    logger.info("... channel selection: "+chan_string)

    if overwrite:
        os.system('rm -rf '+outfile)
        casaStuff.imsubimage(
            imagename=outfile+'.temp',
            outfile=outfile,
            box=box_string,
            chans=chan_string,
            )

    os.system('rm -rf '+outfile+'.temp')

    return(True)

def trim_rind(
        infile=None,
        outfile=None,
        overwrite=False,
        inplace=True,
        pixels=1
        ):

    if infile is None or outfile is None:
        logger.error("Missing required input.")
        return(False)

    if os.path.isdir(infile) == False:
        logger.error("Input file not found: "+infile)
        return(False)

    if inplace == False:
        if os.path.isdir(outfile) or os.path.isfile(outfile):
            if overwrite:
                os.system('rm -rf '+outfile)
            else:
                logger.error("Output file already present: "+outfile)
                return(False)

        os.system('cp -r '+infile+' '+outfile)
        target_file = outfile
    else:
        target_file = infile

    # Figure out the extent of the image inside the cube
    myia = au.createCasaTool(casaStuff.iatool)
    myia.open(target_file)
    mask = myia.getchunk(getmask=True)
    elt = nd.generate_binary_structure(2,1)
    if pixels > 1:
        elt = nd.iterate_structure(elt, pixels-1)
    mask = nd.binary_erosion(mask, elt[:,:,np.newaxis, np.newaxis])
    myia.putregion(pixelmask=mask)
    myia.close()
    return(True)

def primary_beam_correct(
    infile=None,
    pbfile=None,
    outfile=None,
    cutoff=0.25,
    overwrite=False
    ):
    """
    Construct a primary-beam corrected image.
    """

    if infile is None or pbfile is None or outfile is None:
        logger.error("Missing required input.")
        return(False)

    if os.path.isdir(infile) == False:
        logger.error("Input file missing - "+infile)
        return(False)

    if os.path.isdir(pbfile) == False:
        logger.error("Primary beam file missing - "+pbfile)
        return(False)

    if os.path.isfile(outfile) or os.path.isdir(outfile):
        if overwrite:
            os.system('rm -rf '+outfile)
        else:
            logger.error("Output exists and overwrite set to false - "+outfile)
            return(False)

    casaStuff.impbcor(imagename=infile, pbimage=pbfile, outfile=outfile, cutoff=cutoff)

    return(True)

#endregion

#region Convolution and alignment

def align_to_target(
    infile=None,
    outfile=None,
    template=None,
    interpolation='cubic',
    asvelocity=True,
    overwrite=False,
    axes=[-1]
    ):
    """
    Align one cube to another, creating a copy. Right now a thin
    wrapper to imregrid used mostly to avoid exposing CASA directly
    into the postprocessHandler. Might evolve in the future.
    """

    if infile is None or template is None or outfile is None:
        logger.error("Missing required input.")
        return(False)

    if os.path.isdir(infile) == False and os.path.isfile(infile) == False:
        logger.error("Input file missing - "+infile)
        return(False)

    if os.path.isdir(template) == False and os.path.isfile(template) == False:
        logger.error("Template file missing - "+template)
        return(False)

    if os.path.isfile(outfile) or os.path.isdir(outfile):
        if overwrite:
            os.system('rm -rf '+outfile)
        else:
            logger.error("Output exists and overwrite set to false - "+outfile)
            return(False)

    casaStuff.imregrid(
        imagename=infile,
        template=template,
        output=outfile,
        interpolation=interpolation,
        asvelocity=asvelocity,
        axes=axes,
        overwrite=True)

    return(True)

def convolve_to_round_beam(
    infile=None,
    outfile=None,
    force_beam=None,
    overwrite=False
    ):
    """
    Convolve supplied image to have a round beam. Optionally, force
    that beam to some size, else it figures out the beam.
    """

    if infile is None or outfile is None:
        logger.error("Missing required input.")
        return(None)

    if os.path.isdir(infile) == False:
        logger.error("Input file missing - "+infile)
        return(None)

    hdr = casaStuff.imhead(infile)

    if (hdr['axisunits'][0] != 'rad'):
        logger.error("Based on CASA experience. I expected units of radians. I did not find this.")
        logger.error("Adjust code or investigate file "+infile)
        return(None)
    pixel_as = abs(hdr['incr'][0]/np.pi*180.0*3600.)

    if 'perplanebeams' in hdr.keys():
        beamnames = hdr['perplanebeams']['beams'].keys()
        majorlist = [hdr['perplanebeams']['beams'][b]['*0']['major']['value']
                     for b in beamnames]
        bmaj = np.max(majorlist)
    else:
        if (hdr['restoringbeam']['major']['unit'] != 'arcsec'):
            logger.error("Based on CASA experience. I expected units of arcseconds for the beam. I did not find this.")
            logger.error("Adjust code or investigate file "+infile)
            return(None)
        bmaj = hdr['restoringbeam']['major']['value']

    if force_beam is None:
        target_bmaj = np.sqrt((bmaj)**2+(2.0*pixel_as)**2)
    else:
        min_bmaj = np.sqrt((bmaj)**2+(2.0*pixel_as)**2)
        if force_beam < min_bmaj:
            logger.warning("Requested beam is too small for convolution.")
            return(None)
        target_bmaj = force_beam

    casaStuff.imsmooth(imagename=infile,
                       outfile=outfile,
                       targetres=True,
                       major=str(target_bmaj) + 'arcsec',
                       minor=str(target_bmaj) + 'arcsec',
                       pa='0.0deg',
                       overwrite=overwrite
                       )

    # Copy over mask
    copy_mask(infile, outfile, huge_cube_workaround=True)

    return(target_bmaj)

#endregion

#region Units stuff

def calc_jytok(
    hdr=None,
    infile=None
    ):
    """
    Calculate the Jy/beam -> Kelvin conversion. Accepts a header
    already read using imhead or a file name that it will parse.
    """

    c = 2.99792458e10
    kb = 1.380658e-16

    if hdr is None:
        if infile is None:
            logger.error("No header and no infile. Returning.")
            return(None)
        hdr = casaStuff.imhead(infile, mode='list')

    for ii in range(len(hdr['shape'])):
        if hdr['cunit{}'.format(ii+1)] == 'Hz':
            break
        if ii == len(hdr['shape'])-1:
            logger.error("I expected a frequency axis but did not find it. Returning.")
            return(None)
    crpix = hdr['crpix{}'.format(ii+1)]
    cdelt = hdr['cdelt{}'.format(ii+1)]
    crval = hdr['crval{}'.format(ii+1)]
    naxis = hdr['shape'][ii]
    faxis_hz = (np.arange(naxis)+1.-crpix)*cdelt+crval
    freq_hz = np.mean(faxis_hz)

    bmaj_unit = hdr['beammajor']['unit']
    if bmaj_unit != 'arcsec':
        logger.error("Beam unit is not arcsec, which I expected. Returning. Unit instead is "+bmaj_unit)
        return(None)
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
    inplace=False
    ):
    """
    Convert a cube from Jy/beam to K.
    """

    if infile is None or (outfile is None and inplace==False):
        logger.error("Missing required input.")
        return(False)

    if os.path.isdir(infile) == False:
        logger.error("Input file not found: "+infile)
        return(False)

    if inplace == False:
        if os.path.isdir(outfile) or os.path.isfile(outfile):
            if overwrite:
                os.system('rm -rf '+outfile)
            else:
                logger.error("Output file already present: "+outfile)
                return(False)

        os.system('cp -r '+infile+' '+outfile)
        target_file = outfile
    else:
        target_file = infile

    hdr = casaStuff.imhead(target_file, mode='list')
    unit = hdr['bunit']
    if unit != 'Jy/beam':
        logger.error("Input unit is not Jy/beam for file "+target_file+" . Instead found "+unit)
        return(False)

    jytok = calc_jytok(hdr=hdr)

    multiply_cube_by_value(target_file, jytok, brightness_unit='K', huge_cube_workaround=True)

    casaStuff.imhead(target_file, mode='put', hdkey='JYTOK', hdvalue=jytok)

    return(True)

def convert_ktojy(
    infile=None,
    outfile=None,
    overwrite=False,
    inplace=False
    ):
    """
    Convert a cube from K to Jy/beam.
    """

    if infile is None or (outfile is None and inplace==False):
        logger.error("Missing required input.")
        return(False)

    if os.path.isdir(infile) == False:
        logger.error("Input file not found: "+infile)
        return(False)

    if inplace == False:
        if os.path.isdir(outfile) or os.path.isfile(outfile):
            if overwrite:
                os.system('rm -rf '+outfile)
            else:
                logger.error("Output file already present: "+outfile)
                return(False)

        os.system('cp -r '+infile+' '+outfile)
        target_file = outfile
    else:
        target_file = infile

    hdr = casaStuff.imhead(target_file, mode='list')
    unit = hdr['bunit']
    if unit != 'K':
        logger.error("Input unit is not K for file "+target_file+" . Instead found "+unit)
        return(False)

    jytok = calc_jytok(hdr=hdr)

    multiply_cube_by_value(target_file, 1/jytok, 'Jy/beam', huge_cube_workaround=True)

    casaStuff.imhead(target_file, mode='put', hdkey='JYTOK', hdvalue=jytok)

    return(True)

#endregion
