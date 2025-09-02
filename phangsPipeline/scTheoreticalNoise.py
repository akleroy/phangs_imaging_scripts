# Routines to build a theoretical noise map

import logging

import numpy as np

import scipy.ndimage as nd
import scipy.ndimage.morphology as morph
import scipy.stats as ss
from scipy.signal import savgol_coeffs

import astropy.wcs as wcs
import astropy.units as u
from astropy.convolution import convolve, Gaussian2DKernel
from astropy.io import fits
from astropy.stats import mad_std
from astropy.table import Table

from matplotlib import pyplot as plt

from spectral_cube import SpectralCube

from .utilsImages import *

from .pipelineVersion import tableversion, version

np.seterr(divide='ignore', invalid='ignore')

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

def make_theoretical_noise_from_aot(
        infile=None, array='12m', freq_ghz=230.538,
        **kwargs,
        ):
    """Make a theoretical noise map from an AOT export file in decimal
    degrees. Also accepts a list of files to generate an exposure map
    for a mosaic. Extra keyword arguments passed to
    make_theoeratical_noise.

    infile : string or list of AOT-exported files with absolute
    positions in decimal degrees.

    array : '12m' (default) or '7m' used for primary beam

    freq_ghz : used for primary beam

    """

    if isinstance(infile,str):
        # AOT export format
        field_tab = Table.read(infile, format='ascii.csv', data_start=2)

        # Could add some error checking
        field_ra = field_tab['RA'].data
        field_dec = field_tab['Dec'].data    
    elif isinstance(infile,list):
        field_ra_list = []
        field_dec_list = []
        for this_infile in infile:
            this_field_tab = Table.read(
                this_infile, format='ascii.csv', data_start=2)
            this_field_ra = this_field_tab['RA'].data
            this_field_dec = this_field_tab['Dec'].data
            field_ra_list.append(this_field_ra)
            field_dec_list.append(this_field_dec)
        field_ra = np.concatenate(field_ra_list)
        field_dec = np.concatenate(field_dec_list)        
    else:
        return(np.nan)
        
    # Primary beam - default to 12m
    
    beam_fwhm_arcsec = 19.0*(300./freq_ghz)
    if array == '7m':
        beam_fwhm_arcsec = 33.0*(300./freq_ghz)
        
    output = make_theoretical_noise(
        field_ra = field_ra, field_dec = field_dec,
        beam_fwhm_arcsec = beam_fwhm_arcsec,
        **kwargs)

    return(output)

def make_theoretical_noise(
        field_ra = None, field_dec = None,
        beam_fwhm_arcsec = None,
        exposure = None, noiseamp = None,
        template_hdr = None, pix_scale = None,
        support_factor = 1.0,
        outfile = None, overwrite = False,
        asnoise = True, normalize = False,
):
    """Make a theoretical noise map from a set of field centers.

    field_ra : RA of field centers

    field_dec : Dec of field centers

    beam_fwhm_arcsec : the FWHM of the fields in arcsec

    exposure : the exposure time of the pointing. Assume fixed across
    all pointings if not supplied.

    noiseamp : the amplitude of the noise (e.g., due to Tsys
    variations) at each pointing. Assumed fixed if across all
    pointings not supplied.

    template_hdr : a header used to construct the output map. If not
    supplied, one is constructed from the field centers.

    pix_scale : pixel scale of the output image. Only used if a new
    header is constructed. Defaults to a fraction of the primary beam
    if not supplied.

    outfile : if supplied, the output map is written to disk in this
    file.

    overwrite : if True will overwrite previous images. Default False.

    asnoise : returns the results as a predicted noise map. Otherwise,
    returns an exposure map ( propto 1/noise^0.5).

    """

    # Make sure we have a working header
    if template_hdr is None:

        max_ra = np.nanmax(field_ra)
        min_ra = np.nanmin(field_ra)
        max_dec = np.nanmax(field_dec)
        min_dec = np.nanmin(field_dec)

        if pix_scale is None:
            pix_scale = beam_fwhm_arcsec / 10. / 3600.

        pad_arcsec = beam_fwhm_arcsec*2.0
        
        template_hdr = make_simple_header_from_box(
            tlc_coord = (max_ra, max_dec), brc_coord = (min_ra, min_dec),
            pad_arcsec = pad_arcsec, pix_scale_deg = pix_scale,
        )
        
    # Get images of RA and Dec
    ra_img, dec_img = make_axes(header=template_hdr)

    # Make sure we have an initialized set of weights
    if exposure is None:
        exposure = np.ones_like(field_ra, dtype=float)

    if noiseamp is None:
        noiseamp = np.ones_like(field_ra, dtype=float)

    # This is an effective integration map
    weight = exposure/noiseamp**2
    
    # Loop over fields
    exposure_map = np.zeros_like(ra_img, dtype=float)

    # Beam size in sigma units
    beam_sigma_arcsec = beam_fwhm_arcsec / 2.355

    for this_ra, this_dec, this_weight in \
        zip(field_ra, field_dec, weight):
        
        # Assume Euclidean
        offset_arcsec = np.sqrt(np.cos(this_dec/180.*np.pi)**2*(ra_img - this_ra)**2 + \
                                (dec_img - this_dec)**2)*3600.

        # Add the effective expose to the accumulation
        pb_response = np.exp(-0.5*(offset_arcsec/beam_sigma_arcsec)**2)
        this_exposure_map = pb_response**2 * this_weight * \
            (offset_arcsec <= support_factor*beam_fwhm_arcsec)

        #plt.imshow(ra_img, origin='lower')
        #plt.show()
        
        #plt.imshow(offset_arcsec, origin='lower')
        #plt.show()

        #plt.imshow(this_exposure_map, origin='lower')
        #plt.show()
       
        exposure_map += this_exposure_map

    # If requested invert to return as noise
    if asnoise:
        output_map = exposure_map**(-0.5)
    else:
        output_map = exposure_map

    if normalize:
        output_map /= np.nanmax(output_map)
        
    # Return or write
    if outfile is not None:
        this_hdu = fits.PrimaryHDU(output_map, template_hdr)
        this_hdu.writeto(outfile, overwrite=overwrite)
        
    # Return the map
    return(output_map)
