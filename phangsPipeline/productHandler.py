"""productHandler

This module creates signal masks based on image cubes, and then applies
the masks to make moment maps. This is done for each galaxy at multiple
spatial scales.

 Example:
    $ ipython
    from phangsPipeline import keyHandler as kh
    from phangsPipeline import productHandler as prh
    this_kh = kh.KeyHandler(master_key = 'phangsalma_keys/master_key.txt')
    this_prh = prh.ProductHandler(key_handler = this_kh)
    this_prh.set_targets(only = ['ngc0628', 'ngc2997', 'ngc4321'])
    this_prh.set_no_interf_configs(no_interf = False)
    this_prh.set_interf_configs(only = ['7m'])
    this_prh.set_feather_configs(only = ['7m+tp'])
    this_prh.set_line_products(only = ['co21'])
    this_prh.set_no_cont_products(no_cont = True)
    this_prh.loop_product()

 """

import os, sys, re, shutil
import glob
import numpy as np
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
        res_lowresmask = "10p72",
        make_directories=True,
        ):
        """
        """
        if len(self.get_targets()) == 0:            
            logger.error("Need a target list.")
            return(None)
 
        if len(self.get_all_products()) == 0:            
            logger.error("Need a products list.")
            return(None)

        if make_directories:
            self._kh.make_missing_directories(product = True)

        for this_target, this_product, this_config in \
            self.looper(do_targets=True,do_products=True,do_configs=True):

            res_list = self._kh.get_res_for_config(this_config)
            if res_list is None:
                logger.error("No target resolutions found for config"+this_config)
                return()

            indir = self._kh.get_postprocess_dir_for_target(
                target=this_target, changeto=False)
            outdir = self._kh.get_product_dir_for_target(
                target=this_target, changeto=False)

            ### step from build_products_12m.pro; if requested, wipe previous versions of the convolution
            ### step from build_products_12m.pro; convolve to specific resolution
            # done by postprocessHandler.py

            ### step from build_products_12m.pro; generate a low resolution mask from the flat cube
            tag = 'pbcorr_trimmed_k'
            lowres_pbcorr_trimmed_k_file = self._kh.get_cube_filename(
                target = this_target, config = this_config.replace("12m+",""), product = this_product,
                ext = tag+"_res"+res_lowresmask,
                casa = True,
                casaext = '.fits')

            # use 10p72 cube as a "low-resolution" cube
            # if not present, use available lowest resoltuion cube instead
            there = glob.glob(indir+lowres_pbcorr_trimmed_k_file)
            if there:
                lowres_cube_data, lowres_cube_wcs, lowres_cube_noise, lowres_cube_mask = \
                    self.recipe_simple_masking(indir+lowres_pbcorr_trimmed_k_file)
            else:
                for this_res in np.sort(res_list)[::-1]:
                    res_tag = self._kh.get_tag_for_res(this_res)
                    lowres_pbcorr_trimmed_k_file = self._kh.get_cube_filename(
                        target = this_target, config = this_config.replace("12m+",""), product = this_product,
                        ext = tag+"_res"+res_tag,
                        casa = True,
                        casaext = '.fits')
                    there = glob.glob(indir+lowres_pbcorr_trimmed_k_file)
                    if there:
                        lowres_cube_data, lowres_cube_wcs, lowres_cube_noise, lowres_cube_mask = \
                            self.recipe_simple_masking(indir+lowres_pbcorr_trimmed_k_file)
                        break

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
                        self.recipe_simple_masking(indir+pbcorr_trimmed_k_file)

                    ### step from build_products_12m.pro; hybridize the masks
                    # hybridmasked cubes
                    broadcube_data, broadcube_noise = self.recipe_hybrid_masking(
                       cube_data, cube_wcs, cube_noise, cube_mask, lowres_cube_mask,
                       outfitsfile =  outdir+pbcorr_trimmed_k_file.replace(".fits","_hybridmask.fits").replace("_"+tag,""))

                    # signal masked cubes
                    strictcube_data, strictcube_noise = self.recipe_signal_masking(
                       cube_data, cube_wcs, cube_noise, cube_mask,
                       outfitsfile =  outdir + pbcorr_trimmed_k_file.replace(".fits","_signalmask.fits").replace("_"+tag,""))

                    ### step from build_products_12m.pro; collapse into a simple set of moment maps
                    # broad map creation
                    self.recipe_products(
                        cube = broadcube_data,
                        rms = broadcube_noise,
                        commonoutfitsfile = outdir + pbcorr_trimmed_k_file.replace(".fits","_broad.fits").replace("_"+tag,""),
                        do_vmax = False,
                        do_vquad = False)

                    # strict map creation
                    self.recipe_products(
                        cube = strictcube_data,
                        rms = strictcube_noise,
                        commonoutfitsfile = outdir + pbcorr_trimmed_k_file.replace(".fits","_strict.fits").replace("_"+tag,""),
                        do_vmax = False,
                        do_vquad = False)

            ### IDL step; make maps using more sophisticated masking techniques


    #########################
    # recipe_simple_masking #
    #########################

    def recipe_simple_masking(
        self,
        fitsimage = None,
        ):
        """build mask based on the local noise
        """
        hdulist = fits.open(fitsimage)
        cube_data = hdulist[0].data
        cube_wcs = WCS(hdulist[0].header)
        cube_noise = scmasking.noise_cube(cube_data)
        cube_mask = scmasking.simple_mask(cube_data, cube_noise)

        return cube_data, cube_wcs, cube_noise, cube_mask


    #########################
    # recipe_hybrid_masking #
    #########################

    def recipe_hybrid_masking(
        self,
        cube_data,
        cube_wcs,
        cube_noise,
        cube_mask,
        lowres_cube_mask,
        outfitsfile,
        ):
        """build hybrid mask and export to FITS format
        """
        combined_mask = cube_mask + lowres_cube_mask
        hybrid_mask = np.where(combined_mask>=1, 1, 0)
        cube_data_hybridmasked = cube_data * hybrid_mask
        cube_noise_hybridmasked = cube_noise * hybrid_mask

        broadcube_data = SpectralCube(data=cube_data_hybridmasked, wcs=cube_wcs)
        broadcube_noise = SpectralCube(data=cube_noise_hybridmasked, wcs=cube_wcs)

        os.system("rm -rf " + outfitsfile)
        hybrid_mask_spectralcube = SpectralCube(data=hybrid_mask, wcs=cube_wcs)
        hybrid_mask_spectralcube.write(outfitsfile, format="fits")

        return broadcube_data, broadcube_noise


    #########################
    # recipe_signal_masking #
    #########################

    def recipe_signal_masking(
        self,
        cube_data,
        cube_wcs,
        cube_noise,
        cube_mask,
        outfitsfile,
        ):
        """build signal mask and export to FITS format
        """
        cube_data_signalmasked = cube_data * cube_mask
        cube_noise_signalmasked = cube_noise * cube_mask

        strictcube_data = SpectralCube(data=cube_data_signalmasked, wcs=cube_wcs)
        strictcube_noise = SpectralCube(data=cube_noise_signalmasked, wcs=cube_wcs)

        os.system("rm -rf " + outfitsfile)
        cube_mask_spectralcube = SpectralCube(data=cube_mask.astype(int), wcs=cube_wcs)
        cube_mask_spectralcube.write(outfitsfile, format="fits")

        return strictcube_data, strictcube_noise


    ###################
    # recipe_products #
    ###################

    def recipe_products(
        self,
        cube,
        rms,
        commonoutfitsfile,
        do_moment0 = True,
        do_moment1 = True,
        do_moment2 = True,
        do_ew = True,
        do_tmax = True,
        do_vmax = True,
        do_vquad = True,
        do_errormaps = False,
        ):
        """ collapse masked cube
        """
        if do_moment0==True:
            if do_errormaps==True:
                errorfile_moment0 = commonoutfitsfile.replace(".fits","_emom0.fits")
            else:
                errorfile_moment0 = None

            scproduct.write_moment0(
                cube = cube,
                outfile = commonoutfitsfile.replace(".fits","_mom0.fits"),
                errorfile = errorfile_moment0,
                rms = rms)

        if do_moment1==True:
            if do_errormaps==True:
                errorfile_moment1 = commonoutfitsfile.replace(".fits","_emom1.fits")
            else:
                errorfile_moment1 = None
 
            scproduct.write_moment1(
                cube = cube,
                outfile = commonoutfitsfile.replace(".fits","_mom1.fits"),
                errorfile = errorfile_moment1,
                rms = rms)

        if do_moment2==True:
            if do_errormaps==True:
                errorfile_moment2 = commonoutfitsfile.replace(".fits","_emom2.fits")
            else:
                errorfile_moment2 = None
 
            scproduct.write_moment2(
                cube = cube,
                outfile = commonoutfitsfile.replace(".fits","_mom2.fits"),
                errorfile = errorfile_moment2,
		rms = rms)

        if do_ew==True:
            if do_errormaps==True:
                errorfile_ew = commonoutfitsfile.replace(".fits","_eew.fits")
            else:
                errorfile_ew = None
 
            scproduct.write_ew(
                cube = cube,
                outfile = commonoutfitsfile.replace(".fits","_ew.fits"),
                errorfile = errorfile_ew,
                rms = rms)

        if do_tmax==True:
            if do_errormaps==True:
                errorfile_tmax = commonoutfitsfile.replace(".fits","_etmax.fits")
            else:
                errorfile_tmax = None
 
            scproduct.write_tmax(
                cubein = cube,
                outfile = commonoutfitsfile.replace(".fits","_tmax.fits"),
                errorfile = errorfile_tmax,
                rms = rms)

        if do_vmax==True: # <TODO> ValueError: All-NaN slice encountered
            if do_errormaps==True:
                errorfile_vmax = commonoutfitsfile.replace(".fits","_evmax.fits")
            else:
                errorfile_vmax = None
 
            scproduct.write_vmax(
                cubein = cube,
                outfile = commonoutfitsfile.replace(".fits","_vmax.fits"),
                errorfile = errorfile_vmax,
                rms = rms)

        if do_vquad==True: # <TODO> ValueError: All-NaN slice encountered
            if do_errormaps==True:
                errorfile_vquad = commonoutfitsfile.replace(".fits","_evquad.fits")
            else:
                errorfile_vquad = None
 
            scproduct.write_vquad(
                cubein = cube,
                outfile = commonoutfitsfile.replace(".fits","_vquad.fits"),
                errorfile = errorfile_vquad,
                rms = rms)

