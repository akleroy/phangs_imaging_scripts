"""DerivativeHandler

This module creates signal masks based on image cubes, and then applies
the masks to make moment maps. This is done for each galaxy at multiple
spatial scales.

Example:
    $ ipython
    from phangsPipeline import handlerKeys as kh
    from phangsPipeline import handlerDerived as dh
    this_kh = kh.KeyHandler(master_key = 'phangsalma_keys/master_key.txt')
    this_dh = dh.DerivedHandler(key_handler = this_kh)
    this_dh.set_targets(only = ['ngc0628', 'ngc2997', 'ngc4321'])
    this_dh.set_no_interf_configs(no_interf = False)
    this_dh.set_interf_configs(only = ['7m'])
    this_dh.set_feather_configs(only = ['7m+tp'])
    this_dh.set_line_products(only = ['co21'])
    this_dh.set_no_cont_products(no_cont = True)
    this_dh.loop_make_products()

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

casa_enabled = ((sys.argv[0].endswith('start_casa.py'))
                or (sys.argv[0].endswith('casa')))

if casa_enabled:
    logger.debug('casa_enabled = True')
else:
    logger.debug('casa_enabled = False')
    sys.path.append(os.path.dirname(os.path.abspath(__file__)))

# import phangs pipeline stuff
import utilsResolutions
import utilsFilenames
import utilsLines
import handlerTemplate

from scConvolution import smooth_cube
from scNoiseRoutines import recipe_phangs_noise
from scMaskingRoutines import recipe_phangs_strict_mask, recipe_phangs_broad_mask

#import scDerivativeRoutines as scderiv
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

    def loop_derive_products(
        self,
        do_all=False,
        do_convolve=False,
        do_noise=False,
        do_strictmask=False,
        do_broadmask=False,
        do_moments=False,
        make_directories=True, 
        extra_ext_in='', 
        extra_ext_out='', 
        overwrite=True, 
        ):
        """
        Loops over the full set of targets, spectral products (note
        the dual definition of "product" here), and configurations to
        do the imaging. Toggle the parts of the loop using the do_XXX
        booleans. Other choices affect algorithms used.
        """
        
        if do_all:
            do_convolve = True
            do_noise = True
            do_strictmask = True
            do_broadmask = True
            do_moments = True

        # Error checking
        
        if len(self.get_targets()) == 0:            
            logger.error("Need a target list.")
            return(None)
 
        if len(self.get_all_products()) == 0:            
            logger.error("Need a products list.")
            return(None)

        # If requested, make the directories

        if make_directories:
            self._kh.make_missing_directories(derived = True)

        # Convolve the data to all requested angular and physical resolutions.
        
        if do_convolve:

            for this_target, this_product, this_config in \
                    self.looper(do_targets=True,do_products=True,do_configs=True):

                # Always start with the native resolution
                           
                self.task_convolve(
                    target=this_target, config=this_config, product=this_product,
                    just_copy = True, overwrite=overwrite)

                # Loop over all angular and physical resolutions.
                
                res_dict = self._kh.get_ang_res_dict(
                    config=this_config,product=this_product)
                res_list = list(res_dict)
                if len(res_list) > 0:
                    res_list.sort()
                for this_res_tag in res_list:
                    this_res_value = res_dict[this_res_tag]
                    self.task_convolve(
                        target=this_target, config=this_config, product=this_product,
                        res_tag=this_res_tag,res_value=this_res_value,res_type='ang', 
                        overwrite=overwrite)

                res_dict = self._kh.get_phys_res_dict(
                    config=this_config,product=this_product)
                res_list = list(res_dict)
                if len(res_list) > 0:
                    res_list.sort()
                for this_res_tag in res_list:
                    this_res_value = res_dict[this_res_tag]
                    self.task_convolve(
                        target=this_target, config=this_config, product=this_product,
                        res_tag=this_res_tag,res_value=this_res_value,res_type='phys', 
                        overwrite=overwrite)

        # Estimate the noise for each cube.
        
        if do_noise:

            for this_target, this_product, this_config in \
                    self.looper(do_targets=True,do_products=True,do_configs=True):

                # Always start with the native resolution

                self.task_estimate_noise(
                    target=this_target, config=this_config, product=this_product,
                    overwrite=overwrite)

                # Loop over all angular and physical resolutions.
                
                res_dict = self._kh.get_ang_res_dict(
                    config=this_config,product=this_product)
                res_list = list(res_dict)
                if len(res_list) > 0:
                    res_list.sort()
                for this_res_tag in res_list:

                    self.task_estimate_noise(
                        target=this_target, config=this_config, product=this_product,
                        res_tag=this_res_tag, overwrite=overwrite)

                res_dict = self._kh.get_phys_res_dict(
                    config=this_config,product=this_product)
                res_list = list(res_dict)
                if len(res_list) > 0:
                    res_list.sort()
                for this_res_tag in res_list:

                    self.task_estimate_noise(
                        target=this_target, config=this_config, product=this_product,
                        res_tag=this_res_tag, overwrite=overwrite)

        # Make "strict" signal masks for each cube
        
        if do_strictmask:

            for this_target, this_product, this_config in \
                    self.looper(do_targets=True,do_products=True,do_configs=True):

                # Always start with the native resolution

                self.task_build_strict_mask(
                    target=this_target, config=this_config, product=this_product,
                    overwrite=overwrite, res_tag=None)

                # Loop over all angular and physical resolutions.

                for this_res in self._kh.get_ang_res_dict(
                    config=this_config,product=this_product):

                    self.task_build_strict_mask(
                        target=this_target, config=this_config, product=this_product,
                        overwrite=overwrite, res_tag=this_res)

                for this_res in self._kh.get_phys_res_dict(
                    config=this_config,product=this_product):

                    self.task_build_strict_mask(
                        target=this_target, config=this_config, product=this_product,
                        overwrite=overwrite, res_tag=this_res)

        # Make "broad" combination masks.
        
        if do_broadmask:

            for this_target, this_product, this_config in \
                    self.looper(do_targets=True,do_products=True,do_configs=True):
                        
                # Only build one broad mask that covers all resolutions

                self.task_build_broad_mask(
                    target=this_target, config=this_config, product=this_product,
                    overwrite=overwrite, res_tag=None)

        # Make "moments" - derived data products.
        
        if do_moments:

            for this_target, this_product, this_config in \
                    self.looper(do_targets=True,do_products=True,do_configs=True):
                        
                # Always start with the native resolution

                self.task_generate_moments(
                    target=this_target, product=this_product, config=this_config,
                    res_tag=None, overwrite=overwrite)

                # Loop over all angular and physical resolutions.

                for this_res in self._kh.get_ang_res_dict(
                    config=this_config,product=this_product):

                    self.task_generate_moments(
                        target=this_target, product=this_product, config=this_config,
                        res_tag=this_res, overwrite=overwrite)

                for this_res in self._kh.get_phys_res_dict(
                    config=this_config,product=this_product):

                    self.task_generate_moments(
                        target=this_target, product=this_product, config=this_config,
                        res_tag=this_res, overwrite=overwrite)

# end of loop


    ###########################################
    # Defined file names for various products #
    ###########################################

    def _fname_dict(
        self,
        target=None,
        config=None,
        product=None,
        res_tag=None,
        extra_ext_in='',
        extra_ext_out='', 
        ):
        """
        Function to define file names used in other functions.
        """

        # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
        # Error checking
        # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

        if target is None:
            raise Exception("Need a target.")
        if product is None:
            raise Exception("Need a product.")
        if config is None:
            raise Exception("Need a config.")
        
        # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
        # Initialize
        # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

        # The output is a nested dictionary structure, for each cube
        # resolution (res_tag)

        fname_dict = {}

        # Resolution string (if any, not required)
        fname_dict['res_tag'] = res_tag
        if res_tag is None:
            res_tag = ''

        # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
        # Original files
        # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

        orig_filename = utilsFilenames.get_cube_filename(
            target = target,  config = config, product = product,
            ext = 'pbcorr_trimmed_k'+extra_ext_in,
            casa = False)

        fname_dict['orig'] = orig_filename

        # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
        # Output Convolved Cubes
        # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
       
        cube_filename = utilsFilenames.get_cube_filename(
            target = target, config = config, product = product,
            ext = res_tag+extra_ext_out,
            casa = False)

        fname_dict['cube'] = cube_filename

        # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
        # Noise Cubes
        # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

        noise_filename = utilsFilenames.get_cube_filename(
            target = target, config = config, product = product,
            ext = res_tag+extra_ext_out+'_noise',
            casa = False)

        fname_dict['noise'] = noise_filename

        # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
        # Signal Mask
        # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

        # This differs by resolution

        strictmask_filename = utilsFilenames.get_cube_filename(
            target = target, config = config, product = product,
            ext = res_tag+extra_ext_out+'_strictmask',
            casa = False)

        fname_dict['strictmask'] = strictmask_filename

        # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
        # Broad / Hybrid Mask
        # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

        # Note that the broadmask is the same across all resolutions.

        broadmask_filename = utilsFilenames.get_cube_filename(
            target = target, config = config, product = product,
            ext = extra_ext_out+'_broadmask',
            casa = False)

        fname_dict['broadmask'] = broadmask_filename

        # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
        # Moments
        # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
                
        # Just note the root of the moment file name. There are a lot
        # of extensions that will be filled in by the programs.

        moment_root = utilsFilenames.get_cube_filename(
            target=target, config=config, product=product,
            ext=res_tag+extra_ext_out)        
        moment_root = moment_root.replace('.fits','')
        fname_dict['momentroot'] = moment_root
        
        # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
        # Return
        # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
        
        return(fname_dict)
            
    ##################################################################
    # Tasks - discrete steps on target, product, config combinations #
    ##################################################################

    def task_convolve(
        self,
        target = None, 
        config = None, 
        product = None, 
        res_tag = None,
        res_value = None,
        res_type = 'ang',
        just_copy = False,
        extra_ext_in = '', 
        extra_ext_out = '', 
        overwrite = False, 
        tol=0.1,
        nan_treatment='interpolate',
        ):
        """
        Convolve data to lower resolutions. Defaults to copying in some cases.
        """
        
        # Parse the input resolution
        if not just_copy:
            
            if res_value is None:
                logger.warning("Need an input resolution.")
                logger.warning("Defaulting to copy mode.")
                just_copy = True
                
            if res_type.lower() not in ['ang','phys']:
                logger.warning("Input resolution can be angular or physical, ang or phys .")
                logger.warning("Defaulting to copy mode.")
                just_copy = True

            if res_tag is None:
                logger.warning("Need a resolution tag to avoid overlapping file names.")
                logger.warning("Defaulting to copy mode.")
                just_copy = True

        # Generate file names

        indir = self._kh.get_postprocess_dir_for_target(target=target, changeto=False)
        indir = os.path.abspath(indir)+'/'

        outdir = self._kh.get_derived_dir_for_target(target=target, changeto=False)
        outdir = os.path.abspath(outdir)+'/'

        fname_dict_in = self._fname_dict(
            target=target, config=config, product=product, res_tag=None, 
            extra_ext_in=extra_ext_in)

        if just_copy:
            fname_dict_out = self._fname_dict(
                target=target, config=config, product=product, res_tag=None, 
                extra_ext_out=extra_ext_out)
        else:
            fname_dict_out = self._fname_dict(
                target=target, config=config, product=product, res_tag=res_tag, 
                extra_ext_out=extra_ext_out)

        input_file = fname_dict_in['orig']
        outfile = fname_dict_out['cube']

        # Check input file existence        
    
        if not (os.path.isfile(indir+input_file)):
            logger.warning("Missing "+indir+input_file)
            return()

        # Access keywords for mask generation
        
        convolve_kwargs = self._kh.get_derived_kwargs(
            config=config, product=product, kwarg_type='convolve_kw'
            )

        logger.info("")
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&")
        logger.info("Copying or convolving cube for:")
        logger.info(str(target)+" , "+str(product)+" , "+str(config))
        if just_copy:
            logger.info("... mode is just copying.")
        else:
            logger.info("... mode is convolving.")
            logger.info("... to resolution tag: "+res_tag)
            logger.info("... which is resolution type: "+res_type)
            logger.info("... and value: "+str(res_value))
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&")
        logger.info("")
        
        logger.info("Input file "+input_file)
        logger.info("Target file: "+outfile)
        logger.info("Keywords: "+str(convolve_kwargs))
            
        if (not self._dry_run):

            if just_copy:

                if overwrite:
                    os.system('rm -rf '+outdir+outfile)
                
                if (os.path.isfile(outdir+outfile)):
                    logger.warning("Target file already present "+outdir+outfile)
                    return()

                os.system('cp -r '+indir+input_file+' '+outdir+outfile)

            else:
                
                if 'tol' in convolve_kwargs:
                    tol = convolve_kwargs['tol']

                if 'nan_treatment' in convolve_kwargs:
                    nan_treatment = convolve_kwargs['nan_treatment']

                if res_type == 'ang':
                    input_res_value = res_value*u.arcsec
                    smooth_cube(incube=indir+input_file, outfile=outdir+outfile,
                                angular_resolution=input_res_value, 
                                tol=tol, nan_treatment=nan_treatment,
                                overwrite=overwrite)

                if res_type == 'phys':
                    this_distance = self._kh.get_distance_for_target(target)
                    if this_distance is None:
                        logger.error("No distance for target "+target)
                        return()
                    this_distance = this_distance*1e6*u.pc
                    input_res_value = res_value*u.pc
                    smooth_cube(incube=indir+input_file, outfile=outdir+outfile,
                                linear_resolution=input_res_value, distance=this_distance,
                                tol=tol, nan_treatment=nan_treatment,
                                overwrite=overwrite)

        return()

    def task_estimate_noise(
        self,
        target = None, 
        config = None, 
        product = None, 
        res_tag = None,
        extra_ext = '', 
        overwrite = False, 
        ):
        """
        Estimate the noise associated with a data cube and save it to disk.
        """

        # Generate file names

        indir = self._kh.get_derived_dir_for_target(target=target, changeto=False)
        indir = os.path.abspath(indir)+'/'

        outdir = self._kh.get_derived_dir_for_target(target=target, changeto=False)
        outdir = os.path.abspath(outdir)+'/'

        fname_dict = self._fname_dict(
            target=target, config=config, product=product, res_tag=res_tag, 
            extra_ext_in=extra_ext)

        input_file = fname_dict['cube']
        outfile = fname_dict['noise']

        # Check input file existence        
    
        if not (os.path.isfile(indir+input_file)):
            logger.warning("Missing "+indir+input_file)
            return()

        # Access keywords for noise generation
        
        noise_kwargs = self._kh.get_derived_kwargs(
            config=config, product=product, kwarg_type='noise_kw')

        # Report

        logger.info("")
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&")
        logger.info("Running a noise estimate for:")
        logger.info(str(target)+" , "+str(product)+" , "+str(config))
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&")
        logger.info("")
        
        logger.info("Input file "+input_file)
        logger.info("Target file: "+outfile)
        logger.info("Keyword arguments: "+str(noise_kwargs))
            
        # Call noise routines
    
        if (not self._dry_run):
            
            recipe_phangs_noise(
                incube=indir+input_file,
                outfile=outdir+outfile,
                noise_kwargs=noise_kwargs,
                return_spectral_cube=False,
                overwrite=overwrite)

    def task_build_strict_mask(
        self,
        target = None, 
        config = None, 
        product = None, 
        res_tag = None,
        extra_ext = '', 
        overwrite = False, 
        ):
        """
        Estimate the noise associated with a data cube and save it to disk.
        """

        # Generate file names

        indir = self._kh.get_derived_dir_for_target(target=target, changeto=False)
        indir = os.path.abspath(indir)+'/'

        outdir = self._kh.get_derived_dir_for_target(target=target, changeto=False)
        outdir = os.path.abspath(outdir)+'/'

        fname_dict = self._fname_dict(
            target=target, config=config, product=product, res_tag=res_tag, 
            extra_ext_in=extra_ext)

        input_file = fname_dict['cube']
        noise_file = fname_dict['noise']
        outfile = fname_dict['strictmask']

        # Check input file existence        
    
        if not (os.path.isfile(indir+input_file)):
            logger.warning("Missing cube: "+indir+input_file)
            return()

        if not (os.path.isfile(indir+noise_file)):
            logger.warning("Missing noise estimate: "+indir+noise_file)
            return()

        # Access keywords for mask generation
        
        strictmask_kwargs = self._kh.get_derived_kwargs(
            config=config, product=product, kwarg_type='strictmask_kw'
            )

        # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
        # Report
        # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

        logger.info("")
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&")
        logger.info("Creating a strict mask for:")
        logger.info(str(target)+" , "+str(product)+" , "+str(config))
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&")
        logger.info("")
        
        logger.info("Input file "+input_file)
        logger.info("Noise file "+noise_file)
        logger.info("Target file: "+outfile)
        logger.info("Kwargs: "+str(strictmask_kwargs))
            
        # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
        # Call the masking routines
        # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
    
        if (not self._dry_run):
            
            recipe_phangs_strict_mask(
                incube=indir+input_file,
                innoise=indir+noise_file,
                outfile=outdir+outfile,
                mask_kwargs=strictmask_kwargs,
                return_spectral_cube=False,
                overwrite=overwrite)

    def task_build_broad_mask(
        self,
        target = None, 
        config = None, 
        product = None, 
        res_tag = None,
        extra_ext = '', 
        overwrite = False, 
        ):
        """
        Estimate the noise associated with a data cube and save it to disk.
        """

        # Generate file names

        indir = self._kh.get_derived_dir_for_target(target=target, changeto=False)
        indir = os.path.abspath(indir)+'/'

        outdir = self._kh.get_derived_dir_for_target(target=target, changeto=False)
        outdir = os.path.abspath(outdir)+'/'

        fname_dict = self._fname_dict(
            target=target, config=config, product=product, res_tag=res_tag, 
            extra_ext_in=extra_ext)

        input_file = fname_dict['strictmask']
        outfile = fname_dict['broadmask']

        # Check input file existence        
    
        if not (os.path.isfile(indir+input_file)):
            logger.warning("Missing cube: "+indir+input_file)
            return()

        # Access keywords for mask generation
        
        broadmask_kwargs = self._kh.get_derived_kwargs(
            config=config, product=product, kwarg_type='broadmask_kw'
            )

        # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
        # Create the list of masks to combine
        # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

        list_of_masks = []
        
        linked_configs = self._kh.get_linked_mask_configs(
            config=config, product=product)

        if config not in linked_configs:
            linked_configs.append(config)
            
        for cross_config in linked_configs:

            fname_dict = self._fname_dict(
                target=target, config=cross_config, product=product, res_tag=None, 
                extra_ext_in=extra_ext)

            this_mask = fname_dict['strictmask']
            if this_mask not in list_of_masks:
                if os.path.isfile(indir+this_mask):
                    list_of_masks.append(indir+this_mask)
            
            # Loop over all angular and physical resolutions.

            for this_res in self._kh.get_ang_res_dict(
                config=cross_config,product=product):

                fname_dict = self._fname_dict(
                    target=target, config=cross_config, product=product, res_tag=this_res, 
                    extra_ext_in=extra_ext)
                
                this_mask = fname_dict['strictmask']
                if this_mask not in list_of_masks:
                    if os.path.isfile(indir+this_mask):
                        list_of_masks.append(indir+this_mask)

            for this_res in self._kh.get_phys_res_dict(
                config=cross_config,product=product):

                fname_dict = self._fname_dict(
                    target=target, config=cross_config, product=product, res_tag=this_res, 
                    extra_ext_in=extra_ext)
                
                this_mask = fname_dict['strictmask']
                if this_mask not in list_of_masks:
                    if os.path.isfile(indir+this_mask):
                        list_of_masks.append(indir+this_mask)

        # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
        # Report
        # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

        logger.info("")
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&")
        logger.info("Creating a broad mask for:")
        logger.info(str(target)+" , "+str(product)+" , "+str(config))
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&")
        logger.info("")
        
        logger.info("Input file "+input_file)
        logger.info("List of other masks "+str(list_of_masks))
        logger.info("Target file: "+outfile)
        logger.info("Kwargs: "+str(broadmask_kwargs))
            
        # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
        # Call the mask combining routine
        # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
    
        if (not self._dry_run):
            
            recipe_phangs_broad_mask(
                indir+input_file,
                list_of_masks=list_of_masks,
                outfile=outdir+outfile,
                #mask_kwargs=broadmask_kwargs,
                #return_spectral_cube=False,
                overwrite=overwrite)

    def task_generate_moments(
        self,
        target = None, 
        config = None, 
        product = None, 
        res_tag = None, 
        extra_ext = '', 
        overwrite = False, 
        ):
        """
        Generate moment maps.
        """
        # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
        # Look up filenames, list of moments, etc.
        # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

        # Generate file names

        indir = self._kh.get_derived_dir_for_target(target=target, changeto=False)
        indir = os.path.abspath(indir)+'/'

        outdir = self._kh.get_derived_dir_for_target(target=target, changeto=False)
        outdir = os.path.abspath(outdir)+'/'

        # Filenames

        fname_dict_nores = self._fname_dict(
            target=target, config=config, product=product, res_tag=None, 
            extra_ext_in=extra_ext)

        fname_dict = self._fname_dict(
            target=target, config=config, product=product, res_tag=res_tag, 
            extra_ext_in=extra_ext)

        # ... broad mask never has a resolution tag

        broadmask_file = fname_dict_nores['broadmask']

        # ... files with resolution tag

        input_file = fname_dict['cube']
        noise_file = fname_dict['noise']
        strictmask_file = fname_dict['strictmask']

        outroot = fname_dict['momentroot']
        
        # Check input file and mask existence        
        
        if not (os.path.isfile(indir+input_file)):
            logger.warning("Missing cube: "+indir+input_file)
            return()

        found_broadmask = (os.path.isfile(indir+broadmask_file))
        found_strictmask = (os.path.isfile(indir+strictmask_file))

        # Look up which moments to calculate

        list_of_moments = self._kh.get_moment_list(config=config, product=product)

        # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
        # Report
        # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

        logger.info("")
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&")
        logger.info("Generating moment maps for:")
        logger.info(str(target)+" , "+str(product)+" , "+str(config))
        if res_tag is not None:
            logger.info("Resolution "+str(res_tag))
        logger.info("Found a strict mask? "+str(found_strictmask))
        logger.info("Found a broad mask? "+str(found_broadmask))
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&")
        logger.info("")
        
        logger.info("... input file: "+input_file)
        logger.info("... noise file: "+noise_file)
        logger.info("... strict mask file: "+strictmask_file)
        logger.info("... broad mask file: "+broadmask_file)
        logger.info("... list of moments: "+str(list_of_moments))
        logger.info("... output root: "+outroot)
            
        # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
        # Execute
        # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

        if (not self._dry_run):
            
            for this_mom in list_of_moments:
                
                logger.info('... generating moment: '+str(this_mom))

                mom_params = self._kh.get_params_for_moment(this_mom)

                # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                # Look up mask
                # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                if mom_params['mask'] is None:
                    mask_file = None
                elif mom_params['mask'].strip().lower() == 'none':
                    mask_file = None
                elif mom_params['mask'] == 'strictmask':
                    if not found_strictmask:
                        logger.warning("Strict mask needed but not found. Skipping.")
                        continue
                    mask_file = indir+strictmask_file
                elif mom_params['mask'] == 'broadmask':
                    if not found_broadmask:
                        logger.warning("Broad mask needed but not found. Skipping.")
                        continue
                    mask_file = indir+broadmask_file
                else:
                    logger.warning("Mask choice not recognized for moment: "+str(this_mom))
                    logger.warning("Skipping.")
                    continue

                # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                # Check noise
                # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if not (os.path.isfile(indir+noise_file)):
                    logger.warning("Missing noise: "+indir+noise_file)
                    noise_in = None
                    errorfile = None
                else:
                    noise_in = indir+noise_file
                    errorfile = outdir+outroot+mom_params['ext_error']+'.fits'

                # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                # Set up output file
                # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                outfile = outdir+outroot+mom_params['ext']+'.fits'

                # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                # Call the moment generator
                # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                moment_generator(
                    indir+input_file, mask=mask_file, noise=noise_in,
                    moment=mom_params['algorithm'], momkwargs=mom_params['kwargs'],
                    outfile=outfile, errorfile=errorfile,
                    channel_correlation=None)
                
                
