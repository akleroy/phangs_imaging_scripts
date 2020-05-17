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
#from scMoments import moment_generator

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
            do_moment = True

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
                        
                pass

        # Make "moments" - derived data products.
        
        if do_moments:

            for this_target, this_product, this_config in \
                    self.looper(do_targets=True,do_products=True,do_configs=True):
                        
                # Loop over all angular and physical resolutions.

                for this_res in self._kh.get_ang_res_dict(
                    config=this_config,product=this_product):

                    pass

                for this_res in self._kh.get_phys_res_dict(
                    config=this_config,product=this_product):

                    pass

        # Loop over target, product, config combinations

#        for this_target, this_product, this_config in \
#            self.looper(do_targets=True,do_products=True,do_configs=True):
            # do signalmask moment maps for each resolution cube
#            if do_signalmask_moment_maps:
#                for this_res in self._kh.get_res_for_config(this_config):
#                    self.task_generate_moment_maps(target=this_target, product=this_product,
#                                                   config=this_config,
#                                                   res=this_res,
#                                                   extra_ext_in=extra_ext_in,
#                                                   extra_ext_out=extra_ext_out,
#                                                   overwrite=overwrite)
            
            # do hybridmask moment maps for each resolution cube, using a cube close to 10.72 arcsec resolution
#            if do_hybridmask_moment_maps:
#                lowres, lowres_tag = self._find_lowest_res(target=this_target, config=this_config, product=this_product, closeto='10.72arcsec')
#                for this_res in self._kh.get_res_for_config(this_config):
#                    self.task_hybridize_masks(target=this_target, product=this_product, config=this_config, res=this_res, lowres=lowres, extra_ext_in=extra_ext_in, extra_ext_out=extra_ext_out, overwrite=overwrite)
            
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
            ext = extra_ext_out+'_signalmask',
            casa = False)

        fname_dict['broadmask'] = broadmask_filename

        # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
        # Moments
        # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
                
        # Just note the root of the moment file name. There are a lot
        # of extensions that will be filled in by the programs.

        derived_root_strict = utilsFilenames.get_cube_filename(
            target=target, config=config, product=product,
            ext=res_tag+extra_ext_out+'_strictmask')        
        derived_root_strict = derived_root_strict.replace('.fits','')
        fname_dict['momentroot_strict'] = derived_root_strict

        derived_root_broad = utilsFilenames.get_cube_filename(
            target=target, config=config, product=product,
            ext=res_tag+extra_ext_out+'_broadmask')
        derived_root_broad = derived_root_broad.replace('.fits','')
        fname_dict['momentroot_broad'] = derived_root_broad
        
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
            
        if (not self._dry_run):

            if just_copy:

                if overwrite:
                    os.system('rm -rf '+outdir+outfile)
                
                if (os.path.isfile(outdir+outfile)):
                    logger.warning("Target file already present "+outdir+outfile)
                    return()

                os.system('cp -r '+indir+input_file+' '+outdir+outfile)

            else:
                
                if res_type == 'ang':
                    input_res_value = res_value*u.arcsec
                    smooth_cube(incube=indir+input_file, outfile=outdir+outfile,
                                angular_resolution=input_res_value, tol=tol)

                if res_type == 'phys':
                    this_distance = self._kh.get_distance_for_target(target)
                    if this_distance is None:
                        logger.error("No distance for target "+target)
                        return()
                    this_distance = this_distance*1e6*u.pc
                    input_res_value = res_value*u.pc
                    smooth_cube(incube=indir+input_file, outfile=outdir+outfile,
                                linear_resolution=input_res_value, distance=this_distance,
                                tol=tol, overwrite=overwrite)

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

        # Report

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
            
        # Call noise routines
    
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
        
        broadmask_kwargs = self._kh.get_derived_kwargs(
            config=config, product=product, kwarg_type='broadmask_kw'
            )

        # Report

        logger.info("")
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&")
        logger.info("Creating a strict mask for:")
        logger.info(str(target)+" , "+str(product)+" , "+str(config))
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&")
        logger.info("")
        
        logger.info("Input file "+input_file)
        logger.info("Noise file "+noise_file)
        logger.info("Target file: "+outfile)
        logger.info("Kwargs: "+str(broadmask_kwargs))
            
        # Call noise routines
    
        if (not self._dry_run):
            
            recipe_phangs_broad_mask(
                template=indir+input_file,
                list_of_masks=list_of_masks,
                outfile=outdir+outfile,
                mask_kwargs=broadmask_kwargs,
                return_spectral_cube=False,
                overwrite=overwrite)

    def task_generate_moment_maps(
        self,
        target = None, 
        config = None, 
        product = None, 
        res = None, 
        extra_ext_in = '', 
        extra_ext_out = '', 
        overwrite = False, 
        ):
        """
        Placeholder for a task to...
        """
        fname_dict = self._fname_dict(target=target, config=config, product=product, res=res, extra_ext_in=extra_ext_in, extra_ext_out=extra_ext_out)
        # 
        check_existence = True
        if not os.path.isfile(fname_dict['signalmask']):
            check_existence = False
        if not os.path.isfile(fname_dict['derived_strict_map']+'_noise'+'.fits'):
            check_existence = False
        else:
            for tag in ['mom0','mom1','mom2','tpeak','vpeak','ew','vquad']:
                for etag in ['','_error']:
                    if not os.path.isfile(fname_dict['derived_strict_map']+'_'+tag+etag+'.fits'):
                        check_existence = False
                        break
                if not check_existence: 
                    break
        if check_existence and not overwrite:
            logger.info('Found existing moment maps: "'+fname_dict['derived_strict_map']+'*". Will not overwrite.')
            return
        # 
        if not os.path.isfile(fname_dict['in_cube']):
            logger.warning('Input cube file not found: "'+fname_dict['in_cube']+'"')
            return
        # 
        logger.info('Generate moment maps: "'+fname_dict['derived_strict_map']+'*"')
        moment_generator(fname_dict['in_cube'], 
                         root_name = fname_dict['derived_strict_map'], 
                         generate_mask = True, 
                         mask_name = 'signalmask', 
                         generate_noise = True, 
                         )
                         # mask will have a file name: fname_dict['derived_strict_map']+'_signalmask.fits'
                         # noise will have a file name: fname_dict['derived_strict_map']+'_noise.fits'
        output_mask_file = fname_dict['derived_strict_map']+'_signalmask.fits' # according to moment_generator()
        if not os.path.isfile(output_mask_file):
            raise Exception('Error! Failed to run momemnt_generator and produce "'+output_mask_file+'"')
        if output_mask_file != fname_dict['signalmask']:
            shutil.move(output_mask_file, fname_dict['signalmask'])
            if not os.path.isfile(fname_dict['signalmask']):
                raise Exception('Error! Failed to run momemnt_generator and produce "'+fname_dict['signalmask']+'"')
    
    
    def task_hybridize_masks(
        self,
        target = None, 
        config = None, 
        product = None, 
        res = None, 
        lowres = None, 
        extra_ext_in = '', 
        extra_ext_out = '', 
        overwrite = False, 
        ):
        """
        Hybridize each res mask with lowres mask.
        """
        # 
        if lowres is None:
            logger.warning('Low-resolution cube resolution is None. Will not hybridize the mask.')
            return
        # 
        lowres_fname_dict = self._fname_dict(target=target, config=config, product=product, res=lowres, extra_ext_in=extra_ext_in, extra_ext_out=extra_ext_out)
        fname_dict = self._fname_dict(target=target, config=config, product=product, res=res, extra_ext_in=extra_ext_in, extra_ext_out=extra_ext_out)
        # 
        check_existence = True
        if not os.path.isfile(fname_dict['hybridmask']):
            check_existence = False
        #if not os.path.isfile(fname_dict['derived_broad_map']+'_noise'+'.fits'):
        #    check_existence = False
        else:
            for tag in ['mom0','mom1','mom2','tpeak','vpeak','ew','vquad']:
                for etag in ['','_error']:
                    if not os.path.isfile(fname_dict['derived_broad_map']+'_'+tag+etag+'.fits'):
                        check_existence = False
                        break
                if not check_existence: 
                    break
        if check_existence and not overwrite:
            logger.info('Found existing moment maps: "'+fname_dict['derived_broad_map']+'*". Will not overwrite.')
            return
        # 
        if not os.path.isfile(fname_dict['in_cube']):
            logger.warning('Input cube file not found: "'+fname_dict['in_cube']+'"')
            return
        # 
        logger.info('Generate hybridmask moment maps: "'+fname_dict['derived_broad_map']+'*"')
        lowres_mask = fits.getdata(lowres_fname_dict['signalmask'])
        mask, header = fits.getdata(fname_dict['signalmask'], header=True)
        hybridmask = np.logical_or(mask.astype(bool), lowres_mask.astype(bool))
        header['HISTORY'] = ''
        header['HISTORY'] = 'Hybridizing masks "%s" and "%s".'%(lowres_fname_dict['signalmask'], fname_dict['signalmask'])
        header['HISTORY'] = ''
        hybridmask_spectralcube = SpectralCube(data=mask.astype(int), wcs=WCS(header))
        hybridmask_spectralcube.write(fname_dict['hybridmask'], overwrite=True)
        moment_generator(fname_dict['in_cube'], 
                         root_name = fname_dict['derived_broad_map'], 
                         generate_mask = False, 
                         mask = fname_dict['hybridmask'], 
                         generate_noise = False, 
                         rms = fname_dict['derived_strict_map']+'_noise'+'.fits', 
                         )
                         # mask will have a file name: fname_dict['derived_strict_map']+'_signalmask.fits'
                         # noise will have a file name: fname_dict['derived_strict_map']+'_noise.fits'
        
        


        




