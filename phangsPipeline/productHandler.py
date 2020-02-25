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
    this_prh.set_interf_configs(only = ['7m+tp'])
    this_prh.set_line_products(only = ['co21'])
    this_prh.set_no_cont_products(no_cont = True)
    this_prh.loop_product()

 """

import os, sys, re, shutil
import glob
import numpy as np
import pyfits
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

            res_list = self._kh.get_res_for_config(this_config)
            if res_list is None:
                logger.error("No target resolutions found for config"+this_config)
                return()

            #lowest_res = self._kh.get_tag_for_res(res_list.max())
            #print(lowest_res)

            for this_res in res_list:
                res_tag = self._kh.get_tag_for_res(this_res)

                tag = 'pbcorr_trimmed_k'
                pbcorr_trimmed_k_file = self._kh.get_cube_filename(
                    target = this_target, config = this_config, product = this_product,
                    ext = tag+"_res"+res_tag,
                    casa = True,
                    casaext = '.fits')

                indir = self._kh.get_postprocess_dir_for_target(
                    target=this_target, changeto=False)
                outdir = self._kh.get_product_dir_for_target(
                    target=this_target, changeto=False)

                there = glob.glob(indir+pbcorr_trimmed_k_file)
                if there:
                    cube_data, cube_noise, cube_mask = \
                        self.recipe_make_mask_one_beam(indir+pbcorr_trimmed_k_file)

                    cube_data_hybridmasked = cube_data * cube_mask #* lowres_cube_mask
                    cube_noise_hybridmasked = cube_noise * cube_mask #* lowres_cube_mask

                    spectralcube_data = SpectralCube(data=cube_data_hybridmasked) #, wcs=cube_wcs)
                    spectralcube_noise = SpectralCube(data=cube_noise_hybridmasked) #, wcs=cube_wcs)

                    scproduct.write_moment0(
                        cube = spectralcube_data,
                        outfile = outdir + pbcorr_trimmed_k_file.replace(".fits","_test.fits"),
                        rms = spectralcube_noise,
                        )


            """
            this_dir = self._kh.get_postprocess_dir_for_target(
                target=this_target, changeto=True)

            ### IDL step; if requested, wipe previous versions of the convolution
            #listfiles = glob.glob(this_dir + "/" + this_target + "_" + this_config + "_*pc.fits") # <TODO>; need to check sub-directory
            #for i in (len(listfiles)):
            #    os.system("rm -rf " + listfiles[i])

            ### IDL step; convolve to specific resolution
            # postprocessHandler does this

            ### IDL step; generate a low resolution mask from the flat cube
            # <TODO> choose convolved cube here (33 arcsec => clean mask), 10-15 arcsec
            tag = 'pbcorr_trimmed_k'
            pbcorr_trimmed_k_file = self._kh.get_cube_filename(
                target = this_target, config = this_config, product = this_product,
                ext = 'pbcorr_trimmed_k'+extra_ext,
                casa = True,
                casaext = '.image')

            _, _, _, lowres_cube_mask = \
                self.recipe_make_mask_one_beam(pbcorr_trimmed_k_file)

            # <TODO> loop against spatial scale hereafter
            ### IDL step; estimate the noise for each cube
            ### IDL step; build masks holding bright signal at each resolution
            cube_data, cube_wcs, cube_noise, cube_mask = \
                self.recipe_make_mask_one_beam('some_file.fits')

            ### IDL step; hybridize the masks
            cube_data_hybridmasked = cube_data * cube_mask * lowres_cube_mask
            cube_noise_hybridmasked = cube_noise * cube_mask * lowres_cube_mask

            # convert to SpectralCube format
            spectralcube_data = SpectralCube(data=cube_data_hybridmasked, wcs=cube_wcs)
            spectralcube_noise = SpectralCube(data=cube_noise_hybridmasked, wcs=cube_wcs)

            ### IDL step; collapse into a simple set of moment maps
            scproduct.write_moment0(
                cube = spectralcube_data,
                outfile = "test", # <TODO>
                rms = spectralcube_noise,
                )

            ### IDL step; make maps using more sophisticated masking techniques
            # <TODO>; which module should I use?
            """

    #############################
    # recipe_make_mask_one_beam #
    #############################

    def recipe_make_mask_one_beam(
        self,
        fitsimage = None,
        ):
        """build masks holding bright signal at each resolution
        """
        hdulist = pyfits.open(fitsimage)
        cube_data = hdulist[0].data
        cube_wcs = WCS(hdulist[0].header) #, naxis=3)
        cube_noise = scmasking.noise_cube(cube_data)
        cube_mask = scmasking.simple_mask(cube_data, cube_noise)

        return cube_data, cube_noise, cube_mask
