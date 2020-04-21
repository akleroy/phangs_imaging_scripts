"""imagingHandler

This module makes images and cubes out of the uv data products created by uvdataHandler. 

This code needs to be run inside CASA. 

Example:
    $ casa
    from phangsPipeline import handlerKeys as kh
    from phangsPipeline import handlerImaging as imh
    this_kh = kh.KeyHandler(master_key = 'config_keys/master_key.txt')
    this_imh = uvh.ImagingHandler(key_handler = this_kh)
    this_imh.set_targets(only = ['ngc3627'])
    this_imh.loop_imaging()

Args:
    make_dirty_image (boolean): If `True` then make dirty image cube from ms data. Default is `True`.
    revert_to_dirty (boolean): If `True` then reset current image cube by the dirty image cube (mainly for debugging). Default is `False`.
    read_in_clean_mask (boolean): If `True` then read in clean mask as defined in "cleanmask_keys.txt". Default is `True`.
    run_multiscale_clean (boolean): If `True` then run multiscale clean. Default is `True`.
    revert_to_multiscale (boolean): If `True` then reset current image cube by the multiscale cleaned image cube (mainly for debugging). Default is `False`.
    make_singlescale_mask (boolean): If `True` then make singlescale clean mask based on current image cube. Default is `True`.
    run_singlescale_clean (boolean): If `True` then run singlescale clean. Default is `True`.
    do_export_to_fits (boolean): If `True` then export ms data folders into fits-format image cube files. Default is `True`.

Notes:
    This code calls following functions:
        import casaImagingRoutines as imr
        import casaMaskingRoutines as msr
        if make_dirty_image: imr.make_dirty_map             
        if revert_to_dirty: imr.replace_cube_with_copy      
        if read_in_clean_mask: msr.import_and_align_mask    
        if run_multiscale_clean: imr.multiscale_loop        
        if revert_to_multiscale: imr.replace_cube_with_copy 
        if make_singlescale_mask: msr.signal_mask           
        if run_singlescale_clean: imr.singlescale_loop      
        if do_export_to_fits: imr.export_to_fits            

"""

#<DONE># 20200210 dzliu: self._kh._cleanmask_dict is always None. It is not yet implemented in "handlerKeys.py"!
#<TODO># 20200214 dzliu: will users want to do imaging for individual project instead of concatenated ms?
#<TODO># 20200214 dzliu: needing KeyHandler API: 
#<DONE># 20200214 dzliu:     self._kh._cleanmask_dict   --> self._kh.get_cleanmask() # input target name, output clean mask file
#<TODO># 20200214 dzliu:     self._kh._target_dict      --> self._kh.get_target_dict() # for rastring, decstring
#<TODO># 20200214 dzliu:     self._kh._override_dict    --> self._kh.get_overrides()
#<TODO># 20200214 dzliu:     self._kh._dir_keys         --> self._kh.get_target_name_for_multipart_name()
#<TODO># 20200214 dzliu:     self._kh._config_dict['interf_config']['clean_scales_arcsec'] # angular scales
#<TODO># 20200218 dzliu: CASA 5.4.0 works, but CASA 5.6.0 does not work!! -- now should work.
#<TODO># 20200218 dzliu: revert does not work! -- copy_imaging suffix bug fixed.
#<TODO># 20200218 dzliu: need to test 'cont'

import os, sys, re, shutil
import glob
import numpy as np

import logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

casa_enabled = ((sys.argv[0].endswith('start_casa.py'))
                or (sys.argv[0].endswith('casa')))

if casa_enabled:
    logger.debug('casa_enabled = True')
    import casaImagingRoutines as imr
    import casaMaskingRoutines as msr
    reload(imr)
    reload(msr)
else:
    logger.debug('casa_enabled = False')
    sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from clean_call import CleanCall, CleanCallFunctionDecorator

import utilsLines as lines
import handlerTemplate
import utilsFilenames

class ImagingHandler(handlerTemplate.HandlerTemplate):
    """
    Class to makes image cubes out of uv data from each spectral line and continuum of each galaxy. 
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
    
    ###############
    # _fname_dict #
    ###############
    
    def _fname_dict(
        self, 
        product, 
        imagename, 
        ):
        """
        Handles file names used in the imaging processes.
        """
        fname_dict = {}
        is_line_product = product in self._kh.get_line_products()
        if is_line_product:
            fname_dict['root'] = imagename
            fname_dict['suffix'] = ''
            fname_dict['image'] = imagename+'.image'
            fname_dict['model'] = imagename+'.model'
            fname_dict['residual'] = imagename+'.residual'
            fname_dict['mask'] = imagename+'.mask'
            fname_dict['pb'] = imagename+'.pb'
            fname_dict['psf'] = imagename+'.psf'
            fname_dict['weight'] = imagename+'.weight'
            fname_dict['sumwt'] = imagename+'.sumwt'
        else:
            fname_dict['root'] = imagename
            fname_dict['suffix'] = '.tt0'
            fname_dict['image'] = imagename+'.image.tt0'
            fname_dict['model'] = imagename+'.model.tt0'
            fname_dict['residual'] = imagename+'.residual.tt0'
            fname_dict['mask'] = imagename+'.mask.tt0'
            fname_dict['pb'] = imagename+'.pb.tt0'
            fname_dict['psf'] = imagename+'.psf.tt0'
            fname_dict['weight'] = imagename+'.weight.tt0'
            fname_dict['sumwt'] = imagename+'.sumwt.tt0'
            fname_dict['alpha'] = imagename+'.alpha'
            fname_dict['beta'] = imagename+'.beta'
        return fname_dict
    
    
    ################
    # loop_imaging #
    ################
    
    def loop_imaging(
        self, 
        do_all = False,
        do_dirty_image = False,
        do_revert_to_dirty = False,
        do_read_clean_mask = False, 
        do_multiscale_clean = False,
        do_revert_to_multiscale = False,
        do_singlescale_mask = False,
        do_singlescale_clean = False,
        do_revert_to_singlescale = False,
        do_export_to_fits = False, 
        extra_ext_in = None,
        suffix_in = None,
        extra_ext_out = None,        
        recipe = 'phangsalma',
        make_directories = True, 
        dynamic_sizing = True,
        force_square = False,
        overwrite = False, 
        ): 
        """
        Loops over the full set of targets, products, and
        configurations to do the imaging. Toggle the parts of the loop
        using the do_XXX booleans. Other choices affect algorithms
        used.
        """

        if do_all:
            do_dirty_image = True
            # debateable ...
            do_revert_to_dirty = True
            do_read_clean_mask = True
            do_multiscale_clean = True
            # debateable ...
            do_revert_to_multiscale = True
            do_singlescale_mask = True
            do_singlescale_clean = True
            # debateable ...
            do_revert_to_singlescale = True
            do_export_to_fits = True   

        if len(self.get_targets()) == 0:            
            logger.error("Need a target list.")
            return(None)
        
        if len(self.get_all_products()) == 0:            
            logger.error("Need a products list.")
            return(None)
        
        known_recipes = ['phangsalma']
        if recipe not in known_recipes:
            logger.error("Recipe not known "+recipe)
            return(None)

        if make_directories:
            self._kh.make_missing_directories(imaging=True)
        
        # Loop over targets, products, and interferometric configs

        for this_target, this_product, this_config in \
                self.looper(do_targets=True, do_products=True, do_configs=True, just_interf=True):

            # Change to the relevant directory

            this_imaging_dir = self._kh.get_imaging_dir_for_target(this_target, changeto=True)

            # print starting message
            logger.info("")
            logger.info("--------------------------------------------------------")
            logger.info("START: Imaging target "+this_target+" config "+this_config+" product "+this_product)
            logger.info("--------------------------------------------------------")
            logger.info('Imaging recipe: '+recipe)

            if recipe == 'phangsalma':
                self.recipe_phangsalma_imaging(
                    target = this_target, 
                    product = this_product, 
                    config = this_config, 
                    extra_ext_in = extra_ext_in, 
                    suffix_in = suffix_in, 
                    extra_ext_out = extra_ext_out, 
                    do_dirty_image = do_dirty_image, 
                    do_revert_to_dirty = do_revert_to_dirty, 
                    do_read_clean_mask = do_read_clean_mask, 
                    do_multiscale_clean = do_multiscale_clean, 
                    do_revert_to_multiscale = do_revert_to_multiscale, 
                    do_singlescale_mask = do_singlescale_mask, 
                    do_singlescale_clean = do_singlescale_clean, 
                    do_revert_to_singlescale = do_revert_to_singlescale, 
                    do_export_to_fits = do_export_to_fits, 
                    dynamic_sizing = dynamic_sizing, 
                    force_square = force_square, 
                    overwrite = overwrite, 
                    )
                
            # print ending message
            logger.info("--------------------------------------------------------")
            logger.info("END: Imaging target "+this_target+" config "+this_config+" product "+this_product)
            logger.info("--------------------------------------------------------")
        # 
        # end of for looper
    # 
    # end of loop_imaging()

    ####################
    # File name lookup #
    ####################
    
    ##################################################################
    # Tasks - discrete steps on target, product, config combinations #
    ##################################################################
    
    def task_initialize_clean_call(
        self,
        target = None, 
        config = None, 
        product = None, 
        extra_ext_in = None,
        suffix_in = None,
        extra_ext_out = None,
        stage = 'dirty',
        ):
        """
        Initialize a clean call object for a target, config, product
        combination and an imaging stage.
        """

        # check user input. 
        if target is None or config is None or product is None:
            logger.error('Please input target, config and product!')
            logger.error('e.g., target = "ngc3627_1", config = "12m+7m", product = "co21"')
            raise Exception('Please input target, config and product!')

        # Look up the recipes for this case
        recipe_list = self._kh.get_imaging_recipes(config=config, product=product)
        #logger.debug('self._kh._imaging_dict = ' + str(self._kh._imaging_dict))
        #print(self._kh._imaging_dict)
        if recipe_list is None:
            logger.error('Error! Could not get imaging recipe for config '+config+' product '+product+'. Please check your "imaging_recipes.txt".')
            raise Exception('Error! Could not get imaging recipe for config '+config+' product '+product+'. Please check your "imaging_recipes.txt".')

        # Initialize the clean call with the appropriate recipe list
        clean_call = CleanCall(recipe_list)

        # Get the visibility name
        vis_file = utilsFilenames.get_vis_filename(
            target=target, product=product, config=config,
            ext=extra_ext_in, suffix=suffix_in)

        # Test existence
        full_vis_file = self._kh.get_imaging_dir_for_target(target=target)+vis_file
        if not os.path.isdir(full_vis_file):
            logger.error('Visibility file not found: '+full_vis_file+" returning.")
            #raise Exception('Missing visibility needed for imaging.')
            return(None)

        # Set the visibility file (note that we assume we are in the working directory)

        clean_call.set_param('vis', vis_file, nowarning=True)

        # Set the output image file name (note no suffix for imaging root)
         
        image_root = utilsFilenames.get_cube_filename(
            target=target, product=product, config=config,
            ext=extra_ext_out, casa=True, casaext='')

        clean_call.set_param('imagename', image_root, nowarning=True)

        # Get the phase center associated with the target
        rastring, decstring = self._kh.get_phasecenter_for_target(target=target)
        phasecenter = 'J2000 '+rastring+' '+decstring

        clean_call.set_param('phasecenter', phasecenter)

        # Set the rest frequency or reference frequency

        is_line_product = product in self._kh.get_line_products()
        if is_line_product:
            if clean_call.get_param('specmode') != 'cube':
                logger.error('Line product detected but specmode is not cube.')
                raise Exception('Malformed clean call.')
            if stage == 'dirty':
                clean_call.set_param('deconvolver', 'hogbom')
            elif stage == 'multiscale':
                clean_call.set_param('deconvolver', 'multiscale')
            elif stage == 'singlescale':
                clean_call.set_param('deconvolver', 'hogbom')
        else:
            if clean_call.get_param('specmode') != 'mfs':
                logger.error('Continuum product detected but specmode is not msf.')
                raise Exception('Malformed clean call.')
            if stage == 'dirty':
                clean_call.set_param('deconvolver', 'mtmfs')
            elif stage == 'multiscale':
                clean_call.set_param('deconvolver', 'mtmfs')
            elif stage == 'singlescale':
                clean_call.set_param('deconvolver', 'mtmfs')

        if is_line_product:
            this_line_tag = self._kh.get_line_tag_for_line_product(product=product)
            if this_line_tag not in lines.line_list.keys():
                logger.error("Did not find line in line_list "+this_line_tag)
                raise Exception('Malformed clean call.')
            rest_freq_ghz = lines.line_list[this_line_tag]
            clean_call.set_restfreq_ghz(rest_freq_ghz)
        else:
            clean_call.set_reffreq_ghz(None)

        return(clean_call)
    
    @CleanCallFunctionDecorator
    def task_pick_cell_and_imsize(
        self, 
        clean_call = None,
        oversamp=5,
        force_square=False,
        check_files=True,
        ):
        """
        Pick a cell size and imsize for a clean call by inspecting the
        visibility data, then attach these values to the clean call.
        
        This will call casaImagingRoutines.estimate_cell_and_imsize()
        first, then apply custom overrides if they exist.
        """
                
        if clean_call is None:
            logger.warning("Require a clean_call object. Returning.")
            return()

        # Look up the file name
        vis_file = clean_call.get_param('vis')

        # Check visibility existence
        if check_files:
            if not (os.path.isdir(vis_file)):
                logger.warning("Missing visibility data "+vis_file)
                return()

        logger.info("")
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%")
        logger.info("Calculating target cell size and image size for:")
        logger.info(vis_file)
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%")
        logger.info("")
        
        # Call the estimating routine
        if (not self._dry_run) and casa_enabled:
            cell, imsize = \
                imr.estimate_cell_and_imsize(vis_file, 
                                             oversamp,
                                             force_square=force_square)
        else:
            cell, imsize = '0.1arcsec', [1000,1000]
            logger.info('DRY RUN skips calling imr.estimate_cell_and_imsize()')
          
        logger.info('cell='+cell+'arcsec, imsize='+str(imsize))

        clean_call.set_param('cell',cell,nowarning=True)
        clean_call.set_param('imsize',imsize,nowarning=True)

        # Return

        return(cell, imsize)
    
    @CleanCallFunctionDecorator
    def task_assign_multiscales(
        self, 
        clean_call = None,
        config = None, 
        ):
        """
        Look up the appropriate scales for multi-scale clean
        associated with this config.
        """
        
        if clean_call is None:
            logger.warning("Require a clean_call object. Returning.")
            return()

        scales_to_clean = self._kh.get_clean_scales_for_config(config=config)
        clean_call.set_multiscale_arcsec(scales=scales_to_clean)

        return()
    
    @CleanCallFunctionDecorator
    def task_make_dirty_image(
        self, 
        clean_call = None,
        backup = True,
        ):
        """
        Execute the attached clean call with zero iterations to create
        a dirty image.
        
        Set backup to True will make a copy of the dirty image as 
        {imagename}_dirty.image
        """

        logger.info("")
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%")
        logger.info("Making dirty image:")
        logger.info(str(clean_call.get_param('imagename')))
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%")
        logger.info("")
            
        if (not self._dry_run) and casa_enabled:
            imr.make_dirty_image(clean_call)
            if backup:
                imr.copy_imaging(
                    input_root=clean_call.get_param('imagename'),
                    output_root=clean_call.get_param('imagename')+'_dirty',
                    wipe_first=True)
                
        return()
    
    @CleanCallFunctionDecorator
    def task_revert_to_imaging(
        self, 
        clean_call = None, 
        tag = 'dirty', 
        ):
        """
        Reset the current imaging stack to an earlier imaging
        state. Requires a tag that will be used to define the old
        imaging moved into place. Wipes the current imaging.
        """

        logger.info("")
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%")
        logger.info("Resetting to "+tag+" imaging:")
        logger.info(str(clean_call.get_param('imagename')))
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%")
        logger.info("")
            
        if (not self._dry_run) and casa_enabled:
            imr.copy_imaging(
                input_root=clean_call.get_param('imagename')+'_'+tag,
                output_root=clean_call.get_param('imagename'),
                wipe_first=True)

        return()
    
    @CleanCallFunctionDecorator
    def task_read_clean_mask(
        self, 
        clean_call = None, 
        target = None, 
        config = None, 
        product = None, 
        ):
        """
        Identify the clean mask associated with the target and
        product, read it from disk (in FITS form) and align it to the
        astrometry and grid of the current imaging.
        """

        if target is None:
            logger.warning("Require a target. Returning.")
            return()

        #if config is None:
        #    logger.warning("Require a config. Returning.")
        #    return()

        if product is None:
            logger.warning("Require a product. Returning.")
            return()
        
        # Could add an error check here that the template imaging exists

        this_cleanmask = self._kh.get_cleanmask_filename(target = target, product = product)
        if this_cleanmask is None:
            logger.info("No clean mask found for target "+target+" product "+product)
            clean_call.set_param('usemask','pb')
            return()    

        logger.info("")
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%")
        logger.info("Reading and aligning clean mask.")
        logger.info('From target '+target+' product '+product)
        logger.info('To '+clean_call.get_param('imagename'))
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%")
        logger.info("")
        
        if self._dry_run:
            return()
        if not casa_enabled:
            return()

        # Get fname dict
        fname_dict = self._fname_dict(product = product, imagename = clean_call.get_param('imagename'))
        
        # import_and_align_mask
        msr.import_and_align_mask(in_file=this_cleanmask, \
                                  out_file=fname_dict['mask'], \
                                  template=fname_dict['image'])
        clean_call.set_param('usemask','user')

        return()
    
    @CleanCallFunctionDecorator
    def task_multiscale_clean(
        self, 
        clean_call = None, 
        backup = True, 
        ):
        """
        Run a multiscale clean loop to convergence. This task
        currently hardcodes some of the PHANGS-ALMA best choice
        parameters.
        
        Set backup to True will make a copy of the multiscale cleaned 
        image as {imagename}_multiscale.image
        """
        
        if clean_call.get_param('deconvolver') not in ['multiscale','mtmfs']:
            logger.warning("I expected a multiscale or mtmfs deconvolver but got "+str(clean_call.get_param('deconvolver'))+".")
            raise Exception("Incorrect clean call! Should have a multiscale or mtmfs deconvolver.")

        logger.info("")
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%")
        logger.info("Running clean call to convergence for:")
        logger.info(clean_call.get_param('imagename'))
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%")
        logger.info("")
            
        if self._dry_run:
            return()
        if not casa_enabled:
            return()

        imr.clean_loop(clean_call=clean_call, 
                       record_file=clean_call.get_param('imagename')+'_multiscale_record.txt',
                       niter_base_perchan = 10,
                       niter_growth_model = 'geometric', 
                       niter_growth_factor = 2.0, 
                       niter_saturation_perchan = 1000,    
                       niter_other_input = None,
                       cycleniter_base = 100,
                       cycleniter_growth_model='linear',
                       cycleniter_growth_factor = 1.0,
                       cycleniter_saturation_value = 1000,
                       cycleniter_other_input = None,
                       threshold_type = 'snr',
                       threshold_value = 4.0,
                       min_loops = 3,
                       max_loops = 20,
                       max_total_niter = None,
                       convergence_fracflux=0.01,
                       convergence_totalflux=None,
                       convergence_fluxperniter=None,
                       use_absolute_delta=True,
                       stop_at_negative=True,
                       remask_each_loop=False,
                       force_dirty_image=False,    
                       )
                       # log_ext='multiscale',
        if backup:
            imr.copy_imaging(
                input_root=clean_call.get_param('imagename'),
                output_root=clean_call.get_param('imagename')+'_multiscale',
                wipe_first=True)
        
        return()
    
    @CleanCallFunctionDecorator
    def task_singlescale_mask(
        self, 
        clean_call = None, 
        product = None, 
        ):
        """
        Create a signal-to-noise based mask within the existing clean
        mask for deep cleaning. Used before running a deep single
        scale clean.
        """
                
        imagename = clean_call.get_param('imagename')+'.image'
        if not os.path.isdir(imagename):
            logger.error("Image not found: "+imagename)
            return()

        logger.info("")
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%")
        logger.info("Creating signal-to-noise based clean mask for:")
        logger.info(str(clean_call.get_param('imagename')))
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%")
        logger.info("")
        
        if self._dry_run:
            return()
        if not casa_enabled:
            return()

        # Get fname dict
        fname_dict = self._fname_dict(product = product, imagename = clean_call.get_param('imagename'))

        # signal_mask
        msr.signal_mask(cube_root=fname_dict['root'],
                        out_file=fname_dict['mask'],
                        suffix_in=fname_dict['suffix'], 
                        suffix_out=fname_dict['suffix'], 
                        operation='AND',
                        high_snr=4.0,
                        low_snr=2.0,
                        absolute=False)

        clean_call.set_param('usemask','user')

        return()
    
    @CleanCallFunctionDecorator
    def task_singlescale_clean(
        self, 
        clean_call = None, 
        backup = True, 
        ):
        """
        Run a singlescale clean loop to convergence. This task
        currently hardcodes some of the PHANGS-ALMA best choice
        parameters.
        
        Set backup to True will make a copy of the singlescale cleaned 
        image as {imagename}_singlescale.image
        """
        
        if not (clean_call.get_param('deconvolver') in ['hogbom','mtmfs']):
            logger.warning("I expected a singlescale or mtmfs deconvolver but got "+str(clean_call.get_param('deconvolver'))+".")
            raise Exception("Incorrect clean call! Should have a hogbom or mtmfs deconvolver.")
            return()

        logger.info("")
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%")
        logger.info("Running clean call to convergence for:")
        logger.info(clean_call.get_param('imagename'))
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%")
        logger.info("")
            
        if self._dry_run:
            return()
        if not casa_enabled:
            return()

        imr.clean_loop(clean_call=clean_call, 
                       record_file=clean_call.get_param('imagename')+'_singlescale_record.txt',
                       niter_base_perchan = 10,
                       niter_growth_model = 'geometric', 
                       niter_growth_factor = 2.0, 
                       niter_saturation_perchan = 1000,    
                       niter_other_input = None,
                       cycleniter_base = 100,
                       cycleniter_growth_model='linear',
                       cycleniter_growth_factor = 1.0,
                       cycleniter_saturation_value = 1000,
                       cycleniter_other_input = None,
                       threshold_type = 'snr',
                       threshold_value = 1.0,
                       min_loops = 3,
                       max_loops = 20,
                       max_total_niter = None,
                       convergence_fracflux=0.01,
                       convergence_totalflux=None,
                       convergence_fluxperniter=None,
                       use_absolute_delta=True,
                       stop_at_negative=False,
                       remask_each_loop=False,
                       force_dirty_image=False,    
                       )
                       # log_ext='singlescale',
        if backup:
            imr.copy_imaging(
                input_root=clean_call.get_param('imagename'),
                output_root=clean_call.get_param('imagename')+'_singlescale',
                wipe_first=True)
        
        return()
    
    @CleanCallFunctionDecorator
    def task_export_to_fits(
        self, 
        clean_call = None, 
        tag = None, 
        ):
        """ 
        Export the results of a clean_call to FITS files. Optionally
        append _tag to the file.
        """

        image_root = str(clean_call.get_param('imagename'))
        if tag is not None:
            image_root += '_'+tag

        logger.info("")
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%")
        logger.info("Exporting to FITS:"+image_root)
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%")
        logger.info("")
                    
        if (not self._dry_run) and casa_enabled:
            imr.export_imaging_to_fits(image_root)
        
        return()
    
    #############################
    # recipe_imaging_one_target #
    #############################
    
    def recipe_phangsalma_imaging(
        self, 
        target=None,
        product=None,
        config=None,
        extra_ext_in = None,
        suffix_in = None,
        extra_ext_out = None,
        do_dirty_image = True,
        do_revert_to_dirty = True,
        do_read_clean_mask = True, 
        do_multiscale_clean = True,
        do_revert_to_multiscale = True,
        do_singlescale_mask = True,
        do_singlescale_clean = True,
        do_revert_to_singlescale = True,
        do_export_to_fits = True, 
        dynamic_sizing = True,
        force_square = False,
        export_dirty = False,
        export_multiscale = False,
        overwrite = False, 
        ): 
        """
        PHANGS-ALMA basic imaging recipe. 
        
        Major steps:
       
        (1) Dirty imaging
        (2) Align a broad clean mask (if present)
        (3) Lightly masked multiscale clean to S/N ~ 4
        (4) Heavily masked single scale clean until convergence
        (5) Export to FITS

        Operates by passing a "clean_call" object between steps. Can
        restart from any individual step.
        
        Input file names: {target}_{config}_{product}_{extra_ext_in}.ms{.suffix_in}
        
        Output file names: {target}_{config}_{product}_{extra_ext_out}.image
        """

        cell = None
        imsize = None

        # These calls instantiate a clean call object with default
        # clean parameters. For the next few steps, we use the clean
        # call appropriate for stage "dirty" imaging.

        clean_call = self.task_initialize_clean_call(
            target = target, config = config, product = product, 
            extra_ext_in = extra_ext_in, suffix_in = suffix_in,
            extra_ext_out = extra_ext_out,
            stage = 'dirty')

        if clean_call is None:
            logger.warning("I could not make a well-formed clean call.")
            return()
        
        if dynamic_sizing:
            if cell is None or imsize is None:
                cell, imsize = self.task_pick_cell_and_imsize(
                    clean_call=clean_call, force_square=force_square)
            else:                    
                clean_call.set_param('cell',cell,nowarning=True)
                clean_call.set_param('imsize',imsize,nowarning=True)                    

        # Make a dirty image (niter=0)

        if do_dirty_image: 
            
            self.task_make_dirty_image(clean_call=clean_call)

        if do_export_to_fits and export_dirty:

            self.task_export_to_fits(clean_call=clean_call, tag='dirty')

        # Reset the current imaging to the dirty image.

        if do_revert_to_dirty: 
            
            self.task_revert_to_imaging(clean_call=clean_call, tag='dirty')

        # Read and align the clean mask to the astrometry of the image.

        if do_read_clean_mask: 

            self.task_read_clean_mask(
                clean_call = clean_call,
                target = target, config = config, product = product)

        # For the next few steps, we use the clean call appropriate
        # for stage "multiscale" imaging.

        clean_call = self.task_initialize_clean_call(
            target = target, config = config, product = product, 
            extra_ext_in = extra_ext_in, suffix_in = suffix_in,
            extra_ext_out = extra_ext_out,
            stage = 'multiscale')
        
        if not do_read_clean_mask: 
            clean_call.set_param('usemask', 'pb')
            clean_call.set_param('pbmask', 0.2)
        
        if dynamic_sizing:
            clean_call.set_param('cell',cell,nowarning=True)
            clean_call.set_param('imsize',imsize,nowarning=True)                    

        # Look up the angular scales to clean for this config
        
        self.task_assign_multiscales(config=config, clean_call=clean_call)
            
        # Run a multiscale clean until it converges.

        if do_multiscale_clean:
            
            self.task_multiscale_clean(clean_call=clean_call)

        if do_export_to_fits and export_multiscale:

            self.task_export_to_fits(clean_call=clean_call, tag='multiscale')

        # Reset the current imaging to the results of the multiscale clean.

        if do_revert_to_multiscale:
                        
            self.task_revert_to_imaging(clean_call=clean_call, tag='multiscale')

        # For the next few steps, we use the clean call appropriate
        # for stage "singlescale" imaging.

        clean_call = self.task_initialize_clean_call(
            target = target, config = config, product = product, 
            extra_ext_in = extra_ext_in, suffix_in = suffix_in,
            extra_ext_out = extra_ext_out,
            stage = 'singlescale')

        if dynamic_sizing:
            clean_call.set_param('cell',cell,nowarning=True)
            clean_call.set_param('imsize',imsize,nowarning=True)                    

        # Make a signal-to-noise based mask for use in singlescale clean.

        if do_singlescale_mask:
            
            self.task_singlescale_mask(clean_call=clean_call, product = product)

        # Run a singlescale clean until it converges.

        if do_singlescale_clean:

            self.task_singlescale_clean(clean_call=clean_call)

        # Reset the current imaging to the results of the singlescale clean.

        if do_revert_to_singlescale:
                        
            self.task_revert_to_imaging(clean_call=clean_call, tag='singlescale')

        # Export the products of the current clean to FITS files.

        if do_export_to_fits:

            self.task_export_to_fits(clean_call=clean_call)
        
        # Return
        
        return
    
    # end of recipe_phangsalma_imaging()


    
    
    
    
    
    







