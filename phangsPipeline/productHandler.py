"""productHandler

This module creates signal masks based on image cubes, and then applies
the masks to make moment maps. This is done for each galaxy at multiple
spatial scales.

This code needs to be run inside CASA.

There should not be any direct calls to CASA in here. Eventually, this
should be able to run without CASA enabled (though it won't be able to
call any of the CASA-specific routines). Right now, just avoid direct
calls to CASA from this class.
    
 Example:
    $ casa
    from phangsPipeline import keyHandler as kh
    from phangsPipeline import productHandler as prh
    this_kh = kh.KeyHandler(master_key = 'phangsalma_keys/master_key.txt')
    this_prh = prh.ProductHandler(key_handler = this_kh)
    this_prh.set_targets(only = ['ngc4321'])
    this_prh.set_no_interf_configs(no_interf = True)
    this_prh.set_feather_configs(only = ['7m+tp'])
    this_prh.set_line_products(only = ['co21'])
    this_prh.set_no_cont_products(no_cont = True)
    this_prh.loop_product()

 """

import os, sys, re, shutil
import glob
import numpy as np
#import pyfits
from astropy.io import fits
from astropy.wcs import WCS
from spectral_cube import SpectralCube, Projection

import logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

casa_enabled = (sys.argv[0].endswith('start_casa.py'))

if casa_enabled:
    logger.debug('casa_enabled = True')
    import casaCubeRoutines as ccr
    import casaMosaicRoutines as cmr
    import casaFeatherRoutines as cfr
    import casaImagingRoutines as imr
    import casaMaskingRoutines as msr
    reload(imr) #<TODO><DEBUG>#
    reload(msr) #<TODO><DEBUG># #<TODO># we still need imr.cleanCall for dry_run?
    reload(imr) #<TODO><DEBUG>#
else:
    logger.debug('casa_enabled = False')
    sys.path.append(os.path.dirname(os.path.abspath(__file__)))

import utils
import line_list
import handlerTemplate
import scMaskingRoutines as scmasking
import scProductRoutines as scproduct


class ProductHandler(handlerTemplate.HandlerTemplate):
    """
    Class to create signal masks based on image cubes, and then applies
    the masks to make moment maps. This is done for each galaxy at multiple
    spatial scales.
    """
    
    ############
    # __init__ #
    ############
    
    def __init__(
        self,
        key_handler = None,
        dry_run = False,
        ):
        
        # inherit template class
        handlerTemplate.HandlerTemplate.__init__(self,key_handler = key_handler, dry_run = dry_run)
    

    ################
    # loop_product #
    ################

    def loop_product(
        self,
        config_lowresmask = "7m+tp",
        res_lowresmask = "10p72",
        ):
        """
        """
        if len(self.get_targets()) == 0:            
            logger.error("Need a target list.")
            return(None)
 
        if len(self.get_all_products()) == 0:            
            logger.error("Need a products list.")
            return(None)

        for this_target, this_product, this_config in \
            self.looper(do_targets=True,do_products=True,do_configs=True):

            ### step from build_products_12m.pro; if requested, wipe previous versions of the convolution
            ### step from build_products_12m.pro; convolve to specific resolution
            # done by postprocessHandler.py

            res_list = self._kh.get_res_for_config(this_config)
            if res_list is None:
                logger.error("No target resolutions found for config"+this_config)
                return()

            indir = self._kh.get_postprocess_dir_for_target(
                target=this_target, changeto=False)
            outdir = self._kh.get_product_dir_for_target(
                target=this_target, changeto=False)

            ### step from build_products_12m.pro; generate a low resolution mask from the flat cube
            tag = 'pbcorr_trimmed_k'
            lowres_pbcorr_trimmed_k_file = self._kh.get_cube_filename(
                target = this_target, config = config_lowresmask, product = this_product,
                ext = tag+"_res"+res_lowresmask,
                casa = True,
                casaext = '.fits')

            there = glob.glob(indir+lowres_pbcorr_trimmed_k_file)
            if there:
                lowres_cube_data, lowres_cube_wcs, lowres_cube_noise, lowres_cube_mask = \
                    self.recipe_make_mask_one_beam(indir+lowres_pbcorr_trimmed_k_file)
            else:
                res_tag = self._kh.get_tag_for_res(res_list.max()) # <TODO> does this work if config_lowresmask = "12m" or "12m+7m" or "12m+7m+tp"?
                lowres_pbcorr_trimmed_k_file = self._kh.get_cube_filename(
                    target = this_target, config = config_lowresmask, product = this_product,
                    ext = tag+"_res"+res_tag,
                    casa = True,
                    casaext = '.fits')
                lowres_cube_data, lowres_cube_wcs, lowres_cube_noise, lowres_cube_mask = \
                    self.recipe_make_mask_one_beam(indir+lowres_pbcorr_trimmed_k_file)

            ### step from build_products_12m.pro; estimate the noise for each cube
            ### step from build_products_12m.pro; build masks holding bright signal at each resolution
            for this_res in res_list:
                res_tag = self._kh.get_tag_for_res(this_res)

                tag = 'pbcorr_trimmed_k'
                pbcorr_trimmed_k_file = self._kh.get_cube_filename(
                    target = this_target, config = this_config, product = this_product,
                    ext = tag+"_res"+res_tag,
                    casa = True,
                    casaext = '.fits')

                there = glob.glob(indir+pbcorr_trimmed_k_file)
                if there:
                    cube_data, cube_wcs, cube_noise, cube_mask = \
                        self.recipe_make_mask_one_beam(indir+pbcorr_trimmed_k_file)

            ### step from build_products_12m.pro; hybridize the masks
            # hybridmasked cubes
                    combined_mask = cube_mask + lowres_cube_mask
                    hybrid_mask = np.where(combined_mask>=1, 1, 0)
                    cube_data_hybridmasked = cube_data * hybrid_mask
                    cube_noise_hybridmasked = cube_noise * hybrid_mask

                    broadcube_data = SpectralCube(data=cube_data_hybridmasked, wcs=cube_wcs)
                    broadcube_noise = SpectralCube(data=cube_noise_hybridmasked, wcs=cube_wcs)

                    os.system("rm -rf " + outdir + pbcorr_trimmed_k_file.replace(".fits","_hybridmask.fits"))
                    hybrid_mask_spectralcube = SpectralCube(data=hybrid_mask, wcs=cube_wcs)
                    hybrid_mask_spectralcube.write(outdir + pbcorr_trimmed_k_file.replace(".fits","_hybridmask.fits"), format="fits")

            # signal masked cubes
                    cube_data_signalmasked = cube_data * cube_mask
                    cube_noise_signalmasked = cube_noise * cube_mask

                    strictcube_data = SpectralCube(data=cube_data_signalmasked, wcs=cube_wcs)
                    strictcube_noise = SpectralCube(data=cube_noise_signalmasked, wcs=cube_wcs)

                    os.system("rm -rf " + outdir + pbcorr_trimmed_k_file.replace(".fits","_signalmask.fits"))
                    cube_mask_spectralcube = SpectralCube(data=cube_mask.astype(int), wcs=cube_wcs)
                    cube_mask_spectralcube.write(outdir + pbcorr_trimmed_k_file.replace(".fits","_singalmask.fits"), format="fits")

            ### step from build_products_12m.pro; collapse into a simple set of moment maps
            # broad map creation
                    scproduct.write_moment0(
                        cube = broadcube_data,
                        outfile = outdir + pbcorr_trimmed_k_file.replace(".fits","_broad_mom0.fits"),
                        rms = broadcube_noise,
                        )

                    scproduct.write_moment1(
                        cube = broadcube_data,
                        outfile = outdir + pbcorr_trimmed_k_file.replace(".fits","_broad_mom1.fits"),
                        rms = broadcube_noise,
                        )

                    scproduct.write_moment2(
                        cube = broadcube_data,
                        outfile = outdir + pbcorr_trimmed_k_file.replace(".fits","_broad_mom2.fits"),
		        rms = broadcube_noise,
                        )

                    scproduct.write_ew(
                        cube = broadcube_data,
                        outfile = outdir + pbcorr_trimmed_k_file.replace(".fits","_broad_ew.fits"),
                        rms = broadcube_noise,
                        )

                    scproduct.write_tmax(
                        cubein = broadcube_data,
                        outfile = outdir + pbcorr_trimmed_k_file.replace(".fits","_broad_tmax.fits"),
                        rms = broadcube_noise,
                        )

                    """ValueError: All-NaN slice encountered
                    scproduct.write_vmax(
                        cubein = broadcube_data,
                        outfile = outdir + pbcorr_trimmed_k_file.replace(".fits","_broad_vmax.fits"),
                        rms = broadcube_noise,
                        )
                    """
                    """ValueError: All-NaN slice encountered
                    scproduct.write_vquad(
                        cubein = broadcube_data,
                        outfile = outdir + pbcorr_trimmed_k_file.replace(".fits","_broad_vquad.fits"),
                        rms = broadcube_noise,
                        )
                    """ 

            # strict map creation
                    scproduct.write_moment0(
                        cube = strictcube_data,
                        outfile = outdir + pbcorr_trimmed_k_file.replace(".fits","_strict_mom0.fits"),
                        rms = strictcube_noise,
                        )

                    scproduct.write_moment1(
                        cube = strictcube_data,
                        outfile = outdir + pbcorr_trimmed_k_file.replace(".fits","_strict_mom1.fits"),
                        rms = strictcube_noise,
                        )

                    scproduct.write_moment2(
                        cube = strictcube_data,
                        outfile = outdir + pbcorr_trimmed_k_file.replace(".fits","_strict_mom2.fits"),
                        rms = strictcube_noise,
                        )

                    scproduct.write_ew(
                        cube = strictcube_data,
                        outfile = outdir + pbcorr_trimmed_k_file.replace(".fits","_strict_ew.fits"),
                        rms = strictcube_noise,
                        )

                    scproduct.write_tmax(
                        cubein = strictcube_data,
                        outfile = outdir + pbcorr_trimmed_k_file.replace(".fits","_strict_tmax.fits"),
                        rms = strictcube_noise,
                        )

                    """ValueError: All-NaN slice encountered
                    scproduct.write_vmax(
                        cubein = strictcube_data,
                        outfile = outdir + pbcorr_trimmed_k_file.replace(".fits","_strict_vmax.fits"),
                        rms = strictcube_noise,
                        )
                    """
                    """ValueError: All-NaN slice encountered
                    scproduct.write_vquad(
                        cubein = strictcube_data,
                        outfile = outdir + pbcorr_trimmed_k_file.replace(".fits","_strict_vquad.fits"),
                        rms = strictcube_noise,
                        )
                    """ 

            ### IDL step; make maps using more sophisticated masking techniques


    #############################
    # recipe_make_mask_one_beam #
    #############################

    def recipe_make_mask_one_beam(
        self,
        fitsimage = None,
        ):
        """build masks holding bright signal at each resolution
        """
        hdulist = fits.open(fitsimage)
        cube_data = hdulist[0].data
        cube_wcs = WCS(hdulist[0].header) #, naxis=3)
        cube_noise = scmasking.noise_cube(cube_data)
        cube_mask = scmasking.simple_mask(cube_data, cube_noise)

        return cube_data, cube_wcs, cube_noise, cube_mask

