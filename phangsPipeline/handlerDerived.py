"""DerivativeHandler

This module creates signal masks based on image cubes, and then applies
the masks to make moment maps. This is done for each galaxy at multiple
spatial scales.

Example:
    $ ipython
    from phangsPipeline import handlerKeys as kh
    from phangsPipeline import handlerDerived as dh
    this_kh = kh.KeyHandler(master_key = 'phangsalma_keys/master_key.txt')
    this_dh = dh.ProductHandler(key_handler = this_kh)
    this_dh.set_targets(only = ['ngc0628', 'ngc2997', 'ngc4321'])
    this_dh.set_no_interf_configs(no_interf = False)
    this_dh.set_interf_configs(only = ['7m'])
    this_dh.set_feather_configs(only = ['7m+tp'])
    this_dh.set_line_products(only = ['co21'])
    this_dh.set_no_cont_products(no_cont = True)
    this_dh.loop_product()

 """

import os, sys, re, shutil
import glob
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
import astropy.units as u
from spectral_cube import SpectralCube, Projection
from spectral_cube.masks import BooleanArrayMask

import logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

casa_enabled = (sys.argv[0].endswith('start_casa.py'))

if casa_enabled:
    logger.debug('casa_enabled = True')
else:
    logger.debug('casa_enabled = False')
    sys.path.append(os.path.dirname(os.path.abspath(__file__)))

# import utils
import utilsResolutions
import utilsFilenames
import line_list
import handlerTemplate
import scMaskingRoutines as scmasking
import scDerivativeRoutines as scderiv
from scMoments import moment_generator

class DerivedHandler(handlerTemplate.HandlerTemplate):
    """
    Class to create signal masks based on image cubes, and then apply
    the masks to make moment maps. This is done for each galaxy at
    multiple spatial/angular scales.
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
    
    ########################
    # Main processing loop #
    ########################

    def loop_make_products(
        self,
        do_noise=True,
        do_strict_masking=True,
        do_hybrid_masking=True,
        do_moments=True,
        make_directories=True,
        ):
        """
        Loops over the full set of targets, spectral products (note
        the dual definition here), and configurations to do the
        imaging. Toggle the parts of the loop using the do_XXX
        booleans. Other choices affect algorithms used.
        """
        
        # Error checking
        
        if len(self.get_targets()) == 0:            
            logger.error("Need a target list.")
            return(None)
 
        if len(self.get_all_products()) == 0:            
            logger.error("Need a products list.")
            return(None)

        # If requested, make the directories

        if make_directories:
            self._kh.make_missing_directories(product = True)

        # Loop over target, product, config combinations

        for this_target, this_product, this_config in \
            self.looper(do_targets=True,do_products=True,do_configs=True):

            # Either put these flags here or wrap a
            # "phangs_product_recipe" - but try to make the parts much
            # more modular below?
            
            if do_noise:

                # AKL - add the resolution loop HERE

                pass

            if do_strict_masking:

                # AKL - add the resolution loop HERE

                pass

            if do_hybrid_masking:

                # AKL- Resolution loop not needed, I think?

                pass

            if do_moments:

                # AKL - add the resolution loop HERE
                pass

            # generate a low resolution mask from the flat cube and estimate the local noise
            lowres_cube_mask, \
            lowres_cube_tag = \
            self.recipe_building_low_resolution_mask(target=this_target, 
                                                     product=this_product, 
                                                     config=this_config)
            # 
            # build masks for all resolutions
            self.recipe_building_all_resolution_masks(target=this_target, 
                                                      product=this_product, 
                                                      config=this_config, 
                                                      lowres_cube_mask=lowres_cube_mask, 
                                                      lowres_cube_tag=lowres_cube_tag)
            # 
            # end of loop

    ###########################################
    # Defined file names for various products #
    ###########################################
    
    def _fname_dict(
        self,
        target=None,
        config=None,
        product=None,
        extra_ext='',
        res_lowresmask='',
        ):
        """
        Function to define file names used in other functions.
        """

        # Error checking

        if target is None:
            logger.error("Need a target.")
            return()

        if product is None:
            logger.error("Need a product.")
            return()

        if config is None:
            logger.error("Need a config.")
            return()
        
        # The output is a nested dictionary structure, for each cube
        # resolution (res_tag)

        fname_dict = {} 
        
        indir = self._kh.get_postprocess_dir_for_target(target=target, changeto=False)
        indir = os.path.abspath(indir)
        logger.debug('indir: '+indir)
        
        outdir = self._kh.get_derived_dir_for_target(target=target, changeto=False)
        outdir = os.path.abspath(outdir)
        logger.debug('outdir: '+outdir)
        
        if not os.path.isdir(outdir):
            os.makedirs(outdir)
        
        res_list = self._kh.get_res_for_config(config)
        if res_list is None:
            logger.error('No target resolutions found for target '+target+' and config'+config)
            raise Exception('No target resolutions found for target '+target+' and config'+config)
        
        lowest_res_tag = '' # we will find the lowest-resolution cube for res_lowresmask, if res_lowresmask is not given.
        lowest_res = None
        
        for this_res in res_list:
            
            #res_tag = self._kh.get_tag_for_res(this_res)
            res_tag = utilsResolutions.get_tag_for_res(this_res) # this will be like either '5p00' or '80pc'

            cube_filename = utilsFilenames.get_cube_filename(target = target, 
                                                             config = config, 
                                                             product = product,
                                                             ext = 'pbcorr_trimmed_k'+'_res'+res_tag,
                                                             casa = False)
            
            if os.path.isfile(os.path.join(indir, cube_filename)):
                fname_dict[res_tag] = {}
                fname_dict[res_tag]['pbcorr_trimmed_k'] = os.path.join(indir, cube_filename)
                # 
                if lowest_res is None:
                    lowest_res = this_res
                    lowest_res_tag = res_tag
                elif utilsResolutions.get_angular_resolution_for_res(this_res) > utilsResolutions.get_angular_resolution_for_res(lowest_res):
                    lowest_res = this_res
                    lowest_res_tag = res_tag

                # for tag in ['hybridmask', 'signalmask']:
                #     derived_name = utilsFilenames.get_derived_rootname(target=target,
                #                                                        config=config,
                #                                                        product=product,
                #                                                        ext='pbcorr_trimmed_k',
                #                                                        res_tag=res_tag,
                #                                                        derived=tag)
                #     fname_dict[res_tag][tag] = os.path.join(outdir, derived_name)
                for tag in ['broad', 'strict']:
                    derived_name = utilsFilenames.get_derived_rootname(target=target,
                                                                       config=config,
                                                                       product=product,
                                                                       ext='pbcorr_trimmed_k',
                                                                       res_tag=res_tag,
                                                                       derived=tag)
                    fname_dict[res_tag][tag] = os.path.join(
                        outdir, derived_name)
            else:
                # file not found
                this_res_in_arcsec = utilsResolutions.get_angular_resolution_for_res(this_res)
                logger.warning('Cube with tag '+res_tag+' at '+str(np.round(this_res_in_arcsec,2))+' arcsec resolution was not found: "'+os.path.join(indir, cube_filename)+'"')
        
        # if user has input a res_lowresmask, and the file exist, use it, otherwise use the lowest res data

        if (res_lowresmask != '') and (res_lowresmask in fname_dict):
            fname_dict['res_lowresmask'] = res_lowresmask
        else:
            fname_dict['res_lowresmask'] = lowest_res_tag
        
        return(fname_dict)        
        
    ##################################################################
    # Tasks - discrete steps on target, product, config combinations #
    ##################################################################

    def task_estimate_noise(
        self,
        target = None, 
        config = None, 
        product = None,
        res_tag = None,         
        ):
        """
        Placeholder for a task to generate a noise cube.
        """
        pass

    def task_generate_mask(
        self,
        target = None, 
        config = None, 
        product = None, 
        res_tag = None,
        ):
        """
        Placeholder for a task to generate a mask given a file name and noise cube.
        """
        pass

    def task_write_mask(
        self, 
        mask = None, 
        wcs = None, 
        header = None, 
        comments = None, 
        outfile = None, 
        ):
        """
        Write a cube mask 3D array to disk.
        """
        if mask is not None and wcs is not None and outfile is not None:
            if os.path.isfile(outfile):
                os.remove(outfile)
            header_to_write = wcs.to_header()
            if header is not None:
                for key in ['BMAJ', 'BMIN', 'BPA']:
                    if key in header:
                        header_to_write[key] = header[key]
            if comments is not None:
                if np.isscalar(comments):
                    comments = [comments]
                header_to_write['COMMENT'] = ''
                for comment in comments:
                    header_to_write['COMMENT'] = str(comment).strip()
                header_to_write['COMMENT'] = ''
            mask_spectralcube = SpectralCube(data=mask.astype(int), wcs=wcs, header=header_to_write)
            mask_spectralcube.write(outfile, format="fits")

    def recipe_build_strict_moments(self,
                                    target=None,
                                    product=None,
                                    config=None,
                                    ):

        if target is None:
            logger.error('Please input a target.')
            raise Exception('Please input a target.')
        if product is None:
            logger.error('Please input a product.')
            raise Exception('Please input a product.')
        if config is None:
            logger.error('Please input a config.')
            raise Exception('Please input a config.')
        
        # get fname dict for this target and config
        fname_dict = self._fname_dict(target=target,
                                      config=config,
                                      product=product)

        res_list = self._kh.get_res_for_config(config)
        for this_res in res_list:
            res_tag = utilsResolutions.get_tag_for_res(this_res)
            if res_tag in fname_dict:
                cube = SpectralCube.read(
                    fname_dict[res_tag]['pbcorr_trimmed_k'])
                root_name = fname_dict[res_tag]['strict']
                moment_generator(cube,
                                 root_name=root_name, 
                                 generate_mask=True)

    def recipe_build_broad_moments(self,
                                   target=None,
                                   product=None,
                                   config=None,
                                   ):
        if target is None:
            logger.error('Please input a target.')
            raise Exception('Please input a target.')
        if product is None:
            logger.error('Please input a product.')
            raise Exception('Please input a product.')
        if config is None:
            logger.error('Please input a config.')
            raise Exception('Please input a config.')

        fname_dict = self._fname_dict(target=target,
                                      config=config,
                                      product=product)
        res_list = self._kh.get_res_for_config(config)
        for this_res in res_list:
            res_tag = utilsResolutions.get_tag_for_res(this_res)
            if res_tag in fname_dict:
                hires_root_name = fname_dict[res_tag]['strict']
                lores_cube_tag = fname_dict['res_lowresmask']
                if lores_cube_tag == res_tag:
                    continue
                lores_root_name = fname_dict[lores_cube_tag]['strict']

                cube = SpectralCube.read(
                    fname_dict[res_tag]['pbcorr_trimmed_k'])
                broad_mask = scmasking.recipe_hybridize_mask(hires_root_name
                                                             + '_signalmask.fits',
                                                             lores_root_name
                                                             + '_signalmask.fits')
                root_name = fname_dict[res_tag]['broad']
                rms = SpectralCube.read(hires_root_name + '_noise.fits')
                moment_generator(cube,
                                 rms=rms,
                                 mask=broad_mask,
                                 root_name=root_name,
                                 generate_mask=False,
                                 generate_noise=False)

    def recipe_building_low_resolution_mask(
        self,
        target = None, 
        product = None, 
        config = None, 
        res_lowresmask = '10p72',
        ):
        """
        Recipe to build the low-resolution mask.
        """
        # 
        # TODO: we will always try to use 10.72 arcsec resolution cube to 
        #       make the low-resolution mask. If not found then we will use 
        #       the lowest resolution cube to make the low-resolution mask.
        # 
        # check input
        if target is None:
            logger.error('Please input a target.')
            raise Exception('Please input a target.')
        if product is None:
            logger.error('Please input a product.')
            raise Exception('Please input a product.')
        if config is None:
            logger.error('Please input a config.')
            raise Exception('Please input a config.')
        # 
        # get fname dict for this target and config
        fname_dict = self._fname_dict(target = target, 
                                      config = config, 
                                      product = product, 
                                      res_lowresmask = res_lowresmask)
        
        lowres_cube_tag = fname_dict['res_lowresmask']
        lowres_cube_data = fits.getdata(fname_dict[lowres_cube_tag]['pbcorr_trimmed_k'])
        lowres_cube_noise = scmasking.noise_cube(lowres_cube_data)
        lowres_cube_mask = scmasking.simple_mask(lowres_cube_data, lowres_cube_noise)
        
        return lowres_cube_mask, lowres_cube_tag
        
    # AKL - consider to break this apart as suggested above. Try to
    # encapsulate the individual steps in smaller tasks so that we can
    # reuse them in different circumstances.

    def recipe_building_all_resolution_masks(
        self,
        target = None, 
        product = None, 
        config = None, 
        lowres_cube_mask = None, 
        lowres_cube_tag = None, 
        ):
        """
        Recipe to build masks for all resolution cubes.
        """
        # 
        # check input
        if target is None:
            logger.error('Please input a target.')
            raise Exception('Please input a target.')
        if product is None:
            logger.error('Please input a product.')
            raise Exception('Please input a product.')
        if config is None:
            logger.error('Please input a config.')
            raise Exception('Please input a config.')
        if lowres_cube_mask is None or lowres_cube_tag is None:
            lowres_cube_mask, \
            lowres_cube_tag = \
            self.recipe_building_low_resolution_mask(target=this_target, 
                                                     product=this_product, 
                                                     config=this_config)
        # 
        # get fname dict for this target and config
        fname_dict = self._fname_dict(target = target, 
                                      config = config, 
                                      product = product)
        # 
        # get resolution list
        res_list = self._kh.get_res_for_config(config)
        # 
        # loop all resolution cubes
        # build masks holding bright signal at each resolution
        for this_res in res_list:
            # 
            #res_tag = self._kh.get_tag_for_res(this_res)
            res_tag = utilsResolutions.get_tag_for_res(this_res) # this will be like either '5p00' or '80pc'
            # 
            if res_tag in fname_dict:
                # 
                # open cube fits data
                with fits.open(fname_dict[res_tag]['pbcorr_trimmed_k']) as hdulist:
                    cube_data = hdulist[0].data
                    cube_header = hdulist[0].header
                    cube_wcs = WCS(cube_header)
                    # 
                    # build the simple mask and estimate the local noise for this cube
                    cube_noise = scmasking.noise_cube(cube_data)
                    cube_signal_mask = scmasking.simple_mask(cube_data, cube_noise)
                    # 
                    # combine the cube signal mask and low-resolution mask to get the hybridmask
                    cube_hybrid_mask = np.logical_or(cube_signal_mask, lowres_cube_mask) #<TODO>#
                    #cube_hybrid_mask = scmasking.hybridize_mask(hires_in=cube_signal_mask, lores_in=lowres_cube_mask) # this requires SpectralCube type, but the masks are np.array <TODO>
                    # 
                    # write the hybrid (broad) mask to disk
                    self.task_write_mask(\
                        mask = cube_hybrid_mask, 
                        wcs = cube_wcs, 
                        header = cube_header, 
                        comments = ['This hybridmask is made from the OR combination of the signalmask and lowresmask.', 
                                    'This signalmask is made with the '+res_tag+'arcsec resolution cube.',
                                    'This lowresmask is made with the '+lowres_cube_tag+'arcsec resolution cube.'], 
                        outfile = fname_dict[res_tag]['hybridmask'], 
                        )
                        #<TODO># We can add some comments here into the output FITS file header.
                    # 
                    # write the signal (strict) mask to disk
                    self.task_write_mask(
                        mask = cube_signal_mask, 
                        wcs = cube_wcs, 
                        header = cube_header, 
                        comments = ['This signalmask is made with the '+res_tag+'arcsec resolution cube.'], 
                        outfile = fname_dict[res_tag]['signalmask'], 
                        )
                        #<TODO># We can add some comments here into the output FITS file header.
                    # 
                    # apply hybrid (broad) mask to the cube data
                    broadcube_data = SpectralCube(data=cube_data, wcs=cube_wcs, mask=BooleanArrayMask(cube_hybrid_mask, wcs=cube_wcs))
                    broadcube_noise = SpectralCube(data=cube_noise, wcs=cube_wcs, mask=BooleanArrayMask(cube_hybrid_mask, wcs=cube_wcs))
                    # 
                    # apply signal (strict) mask to the cube data
                    strictcube_data = SpectralCube(data=cube_data, wcs=cube_wcs, mask=BooleanArrayMask(cube_signal_mask, wcs=cube_wcs))
                    strictcube_noise = SpectralCube(data=cube_noise, wcs=cube_wcs, mask=BooleanArrayMask(cube_signal_mask, wcs=cube_wcs))
                    # 
                    # make moment maps and other products with the hybrid (broad) mask
                    self.task_write_products(
                        cube = broadcube_data,
                        rms = broadcube_noise,
                        outfile = fname_dict[res_tag]['broad'], 
                        res_tag = res_tag, 
                        )
                    # 
                    # make moment maps and other products with the signal (strict) mask
                    self.task_write_products(
                        cube = strictcube_data,
                        rms = strictcube_noise,
                        outfile = fname_dict[res_tag]['strict'], 
                        res_tag = res_tag, 
                        )
            # 
            

        





