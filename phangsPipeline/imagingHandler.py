"""imagingHandler

This module makes image cubes out of the uv data from each spectral line and continuum of each galaxy. 

This code needs to be run inside CASA. 

Example:
    $ casa
    from phangsPipeline import keyHandler as kh
    from phangsPipeline import imagingHandler as imh
    this_kh = kh.KeyHandler(master_key = 'config_keys/master_key.txt')
    this_imh = uvh.ImagingHandler(key_handler = this_kh)
    this_imh.set_targets(only = ['ngc3627'])
    this_imh.loop_imaging()
    
"""

#<TODO># 20200210 dzliu: self._kh._cleanmask_dict is always None. It is not yet implemented in "keyHandler.py"!
#<TODO># 20200214 dzliu: will users want to do imaging for individual project instead of concatenated ms?
#<TODO># 20200214 dzliu: needing KeyHandler API: 
#<TODO># 20200214 dzliu:     self._kh._cleanmask_dict   --> self._kh.get_cleanmask() # input target name, output clean mask file
#<TODO># 20200214 dzliu:     self._kh._target_dict      --> self._kh.get_target_dict() # for rastring, decstring
#<TODO># 20200214 dzliu:     self._kh._override_dict    --> self._kh.get_overrides()
#<TODO># 20200214 dzliu:     self._kh._dir_keys         --> self._kh.get_target_name_for_multipart_name()
#<TODO># 20200214 dzliu:     self._kh._config_dict['interf_config']['clean_scales_arcsec'] # angular scales

import os, sys, re, shutil
import glob
import numpy as np

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
    from imr import CleanCall
else:
    logger.debug('casa_enabled = False')
    sys.path.append(os.path.dirname(os.path.abspath(__file__)))
    from clean_call import CleanCall

import utils
import line_list
import handlerTemplate


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
    
    
    
    ################
    # loop_imaging #
    ################
    
    def loop_imaging(
        self, 
        target = None, 
        project = None, 
        config = None, 
        product = None, 
        process_individual_mosaic = False, 
        process_individual_project = False, 
        do_chan0 = False, 
        make_dirty_image = True, 
        revert_to_dirty = False, 
        read_in_clean_mask = True, 
        run_multiscale_clean = True, 
        revert_to_multiscale = False, 
        make_singlescale_mask = True, 
        run_singlescale_clean = True, 
        do_export_to_fits = True, 
        forceSquare = False, 
        make_directories = True, 
        overwrite = False, 
        dry_run = False, 
        ): 
        """This function makes an image cube for each product, each config of each target. 
        
        This includes following tasks: 
            task_clean_call()
        
        Args:
            target (list or str): Galaxy name(s).
            config (list or str): Configuration(s) like '12m', '7m', '12m+7m'.
            product (list or str): Product name(s) like 'co21', 'c18o21', 'cont'.
            
        """
        
        if len(self.get_targets()) == 0:            
            logger.error("Need a target list.")
            return(None)
        
        if len(self.get_all_products()) == 0:            
            logger.error("Need a products list.")
            return(None)
        
        if make_directories:
            self._kh.make_missing_directories(imaging=True)
        
        # loop
        for this_target in self.get_targets():
            for this_product in self.get_all_products():
                for this_config in self.get_interf_configs():
                    # 
                    # get imaging dir
                    this_imaging_dir = self._kh.get_imaging_dir_for_target(this_target)
                    # 
                    # set project to empty <TODO>
                    this_project = ''
                    # 
                    # print starting message
                    logger.info("")
                    logger.info("--------------------------------------------------------")
                    logger.info("START: Imaging the data set.")
                    logger.info("--------------------------------------------------------")
                    logger.info('Target: '+this_target)
                    logger.info('Config: '+this_config)
                    logger.info('Product: '+this_product)
                    logger.info('Imaging dir: '+this_imaging_dir)
                    # 
                    # check config suffix if processing individual project ms data
                    if process_individual_project:
                        this_config_suffix = re.sub(r'^(.*?)_([0-9]+)$', r'\2', this_config)
                        this_config = re.sub(r'^(.*?)_([0-9]+)$', r'\1', this_config)
                    # 
                    # 
                    # change directory to imaging dir
                    current_dir = os.getcwd()
                    os.chdir(this_imaging_dir)
                    # 
                    # 
                    # build clean call
                    logger.info('Running recipe_build_clean_call')
                    clean_call = \
                    self.recipe_build_clean_call(
                        target = this_target, 
                        project = this_project, 
                        config = this_config, 
                        product = this_product, 
                        tag = '', 
                        forceSquare = forceSquare, 
                        overwrite = overwrite, 
                        )
                    # 
                    # 
                    # make imaging recipe
                    if clean_call is not None:
                        logger.info('Running recipe_imaging_one_target')
                        self.recipe_imaging_one_target(
                            clean_call = clean_call, 
                            make_dirty_image = make_dirty_image, 
                            revert_to_dirty = revert_to_dirty, 
                            read_in_clean_mask = read_in_clean_mask, 
                            run_multiscale_clean = run_multiscale_clean, 
                            revert_to_multiscale = revert_to_multiscale, 
                            make_singlescale_mask = make_singlescale_mask, 
                            run_singlescale_clean = run_singlescale_clean, 
                            do_export_to_fits = do_export_to_fits, 
                            overwrite = overwrite, 
                            )
                    # 
                    # image chan0 <TODO>
                    #if do_chan0:
                    #    this_product+'_chan0'
                    #    raise NotImplementedError()
                    # 
                    # 
                    # change dir back
                    os.chdir(current_dir)
                    # 
                    # 
                    # print ending message
                    logger.info("--------------------------------------------------------")
                    logger.info("END: Imaging the data set.")
                    logger.info("--------------------------------------------------------")
                # 
                # end of for configs
            # 
            # end of for line or cont products
        # 
        # end of for targets
    # 
    # end of loop_imaging()
    
    
    
    #############################
    # task_pick_cell_and_imsize #
    #############################
    
    def task_pick_cell_and_imsize(
        self, 
        in_file=None,
        oversamp=5,
        forceSquare=False
        ):
        """Pick a cell size and imsize for cleaning. 
        
        This will call casaImagingRoutines.estimate_cell_and_imsize() first, 
        then apply custom overrides if any.
        """
        
        if (not self._dry_run) and casa_enabled:
            cell_size_string, x_size_string, y_size_string = \
                imr.estimate_cell_and_imsize(in_file, 
                                             oversamp,
                                             forceSquare=forceSquare)
        else:
            cell_size_string, x_size_string, y_size_string = \
                '0.1', '1000', '1000'
            logger.info('DRY RUN skips calling imr.estimate_cell_and_imsize()')
                
        # Check for overrides
        #<TODO> should we remove the '.ms' suffix?
        logger.debug('Checking overrides for "'+in_file+'"')
        if self._kh.has_overrides_for_key(in_file):
            cell_size_string = self._kh.get_overrides(key = in_file, param = 'cell_size', default = cell_size_string)
            x_size_string = self._kh.get_overrides(key = in_file, param = 'x_size', default = x_size_string)
            y_size_string = self._kh.get_overrides(key = in_file, param = 'y_size', default = y_size_string)
        return cell_size_string, x_size_string, y_size_string
    
    
    
    ###########################
    # recipe_build_clean_call #
    ###########################
    
    def recipe_build_clean_call(
        self, 
        target = None, 
        project = None, 
        config = None, 
        product = None, 
        tag = '', 
        forceSquare = False, 
        overwrite = False, 
        dry_run = False, 
        ):
        """Build a clean call before running the clean task with the function phangsImagingRecipe().
        """
        
        # 
        # This code is adapted from the function buildPhangsCleanCall() in phangsPipeline.py / imagingPipeline.py.
        # 
        
        # 
        # check user input. 
        if target is None or config is None or product is None:
            logger.error('Please input target, config and product! (e.g., target = "ngc3627_1", config = "12m+7m", product = "co21")')
            raise Exception('Please input target, config and product!')
        # 
        # Initialize the call
        clean_call = CleanCall()
        # 
        # Set ms data file name
        if project is None or project == '':
            clean_call.vis = target+'_'+config+'_'+product+'.ms'
        else:
            clean_call.vis = target+'_'+project+'_'+config+'_'+product+'.ms'
        # 
        # append tag
        if tag == '':
            clean_call.image_root = target+'_'+config+'_'+product
        else:
            clean_call.image_root = target+'_'+config+'_'+product+'_'+tag
        # 
        # select antenna
        if config == '12m+7m':
            clean_call.antenna = ''
            #clean_call.antenna = select_12m7m
        if config == '12m':
            clean_call.antenna = ''
            #clean_call.antenna = select_12m
        if config == '7m':
            clean_call.antenna = ''
            #clean_call.antenna = select_7m
        # 
        # Look up the center and shape of the mosaic
        mosaic_key = self._kh._target_dict # read_mosaic_key()
        this_ra = mosaic_key[target]['rastring']
        this_dec = mosaic_key[target]['decstring']
        clean_call.phase_center = 'J2000 '+this_ra+' '+this_dec
        # 
        # check ms data file
        if not os.path.isdir(clean_call.vis):
            logger.warning("Visibility data "+'"'+os.getcwd()+os.sep+clean_call.vis+'"'+" not found. Returning empty.")
            return None
        # 
        # get cell size and imsize
        cell_size, x_size, y_size = \
            self.task_pick_cell_and_imsize(\
                clean_call.vis, 
                forceSquare=forceSquare) #<TODO># dzliu: this can be improved
        image_size = [int(x_size), int(y_size)]
        
        clean_call.cell_size = cell_size
        clean_call.image_size = image_size
        
        # Look up the line and data product
        
        if product == 'co21':
            clean_call.specmode = 'cube'
            clean_call.restfreq_ghz = line_list.line_list['co21']

        if product == 'co21_chan0':
            clean_call.specmode = 'mfs'
            clean_call.restfreq_ghz = line_list.line_list['co21']

        if product == 'c18o21':
            clean_call.specmode = 'cube'
            clean_call.restfreq_ghz = line_list.line_list['c18o21']

        if product == 'c18o21_chan0':
            clean_call.specmode = 'mfs'
            clean_call.restfreq_ghz = line_list.line_list['c18o21']

        if product == 'cont':
            clean_call.specmode = 'mfs'
            clean_call.restfreq_ghz = -1.0
        
        # Set angular scales to be used in multiscale clean
        #scales_for_clean = self._config_dict['interf_config'][this_config]['clean_scales_arcsec']
        
        if config == '7m':
            clean_call.pblimit = 0.25
            clean_call.smallscalebias = 0.6
            clean_call.scales_as_angle = [0, 5, 10]
        elif config == '12m':
            clean_call.smallscalebias = 0.6
            clean_call.scales_as_angle = [0, 1, 2.5, 5]
        elif config == '12m+7m':
            clean_call.smallscalebias = 0.8
            clean_call.scales_as_angle = [0, 1, 2.5, 5, 10]
        
        # Look up overrides in the imaging parameters
        override_dict = self._kh._override_dict # read_override_imaging_params()

        if clean_call.image_root in override_dict:
            this_override_dict = override_dict[clean_call.image_root]

            if 'smallscalebias' in this_override_dict:
                clean_call.smallscalebias = float(this_override_dict['smallscalebias'])
            if 'x_size' in this_override_dict:
                clean_call.image_size[0] = int(this_override_dict['x_size'])
            if 'y_size' in this_override_dict:
                clean_call.image_size[1] = int(this_override_dict['y_size'])
            if 'pblimit' in this_override_dict:
                clean_call.pblimit = float(this_override_dict['pblimit'])
            if 'scales_as_angle' in this_override_dict:
                scales_as_angle_string = this_override_dict['scales_as_angle']
                tokens = scales_as_angle_string.split(',')
                scales_as_angle = []
                for token in tokens:
                    if token == '':
                        continue
                    scales_as_angle.append(float(token))
                clean_call.scales_as_angle = scales_as_angle
        
        # Define the clean mask (note one mask per galaxy)
        
        dir_key = self._kh._dir_keys # read_dir_key() #<TODO># Need KeyHandler function get_target_name_by_multipart_name()
        if target in dir_key:
            this_gal = dir_key[target]
        else:
            this_gal = target
        
        cleanmask_dict = self._kh._cleanmask_dict
        if cleanmask_dict is None:
            cleanmask_dict = {}
        if this_gal in cleanmask_dict:
            this_cleanmask = cleanmask_dict[this_gal]
        else:
            logger.error('Error! Clean mask is not defined for target "'+this_gal+'" in "cleanmask_key.txt"! cleanmask_dict: '+str(cleanmask_dict.keys()))
            #raise Exception('Error! Clean mask is not defined for target "'+this_gal+'" in "cleanmask_key.txt"!')
            #this_cleanmask = '../clean_masks/'+this_gal+'_co21_clean_mask.fits' #<TODO># this needs discussion and improvement
            this_cleanmask = ''
            #<TODO># 20200210 dzliu: self._kh._cleanmask_dict is always None. It is not yet implemented in "keyHandler.py"!
        
        if os.path.isfile(this_cleanmask):
            clean_call.clean_mask_file = this_cleanmask
        else:
            clean_call.clean_mask_file = None
            logger.warning('Warning! Clean mask for target "'+this_gal+'" was not found: "'+this_cleanmask+'"')
        
        # 
        # Return
        return clean_call
    
    # end of recipe_build_clean_call()
    
    
    
    
    #############################
    # recipe_imaging_one_target #
    #############################
    
    def recipe_imaging_one_target(
        self, 
        clean_call = None, 
        make_dirty_image = True, 
        revert_to_dirty = False, 
        read_in_clean_mask = True, 
        run_multiscale_clean = True, 
        revert_to_multiscale = False, 
        make_singlescale_mask = True, 
        run_singlescale_clean = True, 
        do_export_to_fits = True, 
        overwrite = False, 
        ): 
        """PHANGS imaging recipe. 
        
        Including steps:
            Dirty image -> 
            mask alignment -> 
            lightly masked multiscale clean -> 
            heavily masked single scale clean -> 
            export.

        """
        
        # 
        # This code is adapted from the function phangsImagingRecipe() in phangsPipeline.py / imagingPipeline.py.
        # 
        
        if make_dirty_image:
            if (not self._dry_run) and casa_enabled:
                logger.info("")
                logger.info("MAKING THE DIRTY IMAGE.")
                logger.info("")
                imr.make_dirty_map(clean_call)
            else:
                logger.info('DRY RUN skips outputting to '+clean_call.image_root+'_dirty'+'.*')

        if revert_to_dirty:
            if (not self._dry_run) and casa_enabled:
                logger.info("")
                logger.info("RESETING THE IMAGING TO THE DIRTY IMAGE.")
                logger.info("")
                imr.replace_cube_with_copy(
                    to_root=clean_call.image_root,
                    from_root=clean_call.image_root+'_dirty')
            else:
                logger.info('DRY RUN skips outputting to '+clean_call.image_root+'.*')

        if read_in_clean_mask:
            if (not self._dry_run) and casa_enabled:
                logger.info("")
                logger.info("READING IN THE CLEAN MASK.")
                logger.info("")
                if clean_call.clean_mask_file is not None:
                    msr.import_and_align_mask(
                        in_file=clean_call.clean_mask_file,
                        out_file=clean_call.image_root+'.mask',
                        template=clean_call.image_root+'.image',
                        )
                    clean_call.usemask = 'user'
                else:
                    logger.info("No clean mask defined.")
                    clean_call.usemask = 'pb'
            else:
                logger.info("DRY RUN skips setting the clean mask")
                clean_call.usemask = 'pb'
        
        if run_multiscale_clean:
            if (not self._dry_run) and casa_enabled:
                logger.info("")
                logger.info("RUNNING THE MULTISCALE CLEAN.")
                logger.info("")
                imr.multiscale_loop(
                    clean_call = clean_call,
                    record_file = clean_call.image_root+'_multiscale_record.txt',
                    delta_flux_threshold=0.01,
                    absolute_threshold=None,
                    snr_threshold=4.0,
                    stop_at_negative=True,
                    max_loop = 20
                    )
            else:
                logger.info("DRY RUN skips running imr.multiscale_loop() and outputting to "+clean_call.image_root+'_multiscale'+'.*')

        if revert_to_multiscale:
            if (not self._dry_run) and casa_enabled:
                logger.info("")
                logger.info("RESETING THE IMAGING TO THE OUTPUT OF MULTISCALE CLEAN.")
                logger.info("")
                imr.replace_cube_with_copy(
                    to_root=clean_call.image_root,
                    from_root=clean_call.image_root+'_multiscale')
            else:
                logger.info("DRY RUN skips outputting to "+clean_call.image_root+'_multiscale'+'.*')

        if make_singlescale_mask:
            if (not self._dry_run) and casa_enabled:
                logger.info("")
                logger.info("MAKING THE MASK FOR SINGLE SCALE CLEAN.")
                logger.info("")
                msr.signal_mask(
                    cube_root=clean_call.image_root,
                    out_file=clean_call.image_root+'.mask',
                    operation='AND',
                    high_snr=4.0,
                    low_snr=2.0,
                    absolute=False)
                clean_call.usemask='user'
            else:
                logger.info("DRY RUN skips running msr.signal_mask() and setting the clean mask")
            
        if run_singlescale_clean:
            if (not self._dry_run) and casa_enabled:
                logger.info("")
                logger.info("RUNNING THE SINGLE SCALE CLEAN.")
                logger.info("")
                imr.singlescale_loop(
                    clean_call = clean_call,
                    record_file = clean_call.image_root+'_singlescale_record.txt',
                    delta_flux_threshold=0.01,
                    absolute_delta=True,
                    absolute_threshold=None,
                    snr_threshold=1.0,
                    stop_at_negative=False,
                    remask=False,
                    max_loop = 20
                    )
            else:
                logger.info("DRY RUN skips running imr.singlescale_loop() and outputting to "+clean_call.image_root+'_singlescale'+'.*')

        if do_export_to_fits:
            if (not self._dry_run) and casa_enabled:
                logger.info("")
                logger.info("EXPORTING PRODUCTS TO FITS.")
                logger.info("")
                imr.export_to_fits(clean_call.image_root)
                imr.export_to_fits(clean_call.image_root+'_dirty')
                imr.export_to_fits(clean_call.image_root+'_multiscale')
            else:
                logger.info("DRY RUN skips running imr.export_to_fits() and outputting to "+clean_call.image_root+'*.fits')
        
        return
    
    # end of recipe_imaging_one_target()
    
    


    
    
    
    
    
    







