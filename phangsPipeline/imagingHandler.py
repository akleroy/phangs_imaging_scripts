"""imagingHandler

This module makes image cubes out of the uv data from each spectral line and continuum of each galaxy. 

This code needs to be run inside CASA. 

Example:
    $ casa
    from phangsPipeline import keyHandler as kh
    from phangsPipeline import imagingHandler as imh
    this_kh = kh.KeyHandler(master_key = 'config_keys/master_key.txt')
    this_imh = uvh.ImagingHandler(key_handler = this_kh)
    this_imh.image_data(targets = ['ngc3627'])
    
"""

import os, sys, re, shutil
import glob
import numpy as np

import logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

import utils
import line_list

casa_enabled = True

if casa_enabled:
    import casaCubeRoutines as ccr
    import casaMosaicRoutines as cmr
    import casaFeatherRoutines as cfr
    import casaImagingRoutines as imr
    import casaMaskingRoutines as msr


class imagingHandler:
    """Class to makes image cubes out of uv data from each spectral line and continuum of each galaxy. 
    
    THIS IS JUST A SKETCH OF THIS CODE. FUNCTIONS ARE NOT IMPLEMENTED YET!
    """
    
    ############
    # __init__ #
    ############
    
    def __init__(
        self, 
        key_handler = None, 
        ): 
        self._kh = key_handler
        self._cc = None # clean_call
        pass
    
    ########################################
    # Set up loop and processing variables #
    ########################################

    def set_targets(self, ): pass
    def set_line_products(self, ): pass
    def set_cont_products(self, ): pass
    def set_interf_configs(self, ): pass
    def set_no_line(self, ): pass
    def set_no_cont(self, ): pass
    def set_dry_run(self, ): pass
    def set_key_handler(self, ): pass
    def _build_lists(self, ): pass
    def _all_products(self, ): pass

    ##############
    # image_data #
    ##############
    
    # was image_data.py
    def loop_image_data(self, **kwargs):
        """Loop and make image cubes for all uv data.
        """
        # 
        self._kh.set_looper(**kwargs) # this KeyHandler.set_looper() can be a new API in KeyHandler to generate the loop list of (target,product,config) tuple.
        # 
        for this_target, this_product, this_config in self._kh.get_looper(): # this KeyHandler.get_looper() can be a new API in KeyHandler corresponding to KeyHandler.set_looper().
            # 
            this_imaging_dir = self._kh.get_imaging_dir(this_target)
            current_dir = os.getcwd()
            # 
            os.chdir(this_imaging_dir)
            logger.info("START: Imaging the data set.")
            # 
            self.buildPhangsCleanCall() # build clean call. we can rename it to something like build_clean_call.
            # 
            TODO_decide_which = False
            if TODO_decide_which:
                # 
                # We can either call the old function phangsImagingRecipe() as in phangsPipeline.py, 
                # 
                self.phangsImagingRecipe()
            else:
                # 
                # Or discard that function and move its content directly here. 
                # Note that imr = casaImagingRoutines and msr = casaMaskingRoutines.
                # 
                imr.make_dirty_map(clean_call)
                # 
                imr.replace_cube_with_copy('', '_dirty')
                # 
                if read_in_clean_mask:
                    msr.import_and_align_mask()
                    clean_call.usemask = 'user'
                else:
                    clean_call.usemask = 'pb'
                # 
                if run_multiscale_clean:
                    imr.multiscale_loop()
                # 
                if revert_to_multiscale:
                    imr.replace_cube_with_copy('', '_multiscale')
                # 
                if make_singlescale_mask:
                    msr.signal_mask()
                    clean_call.usemask='user'
                # 
                if run_singlescale_clean:
                    imr.singlescale_loop()
                # 
                if do_export_to_fits:
                    imr.export_to_fits('')
                    imr.export_to_fits('_dirty')
                    imr.export_to_fits('_multiscale')
            # 
            logger.info("END: Imaging the data set.")
            os.chdir(current_dir)
        # 
        return()
    
    # was phangsPipeline.py pick_phangs_cell_and_imsize function
    def pick_phangs_cell_and_imsize(self, ):
        """Calculates cell and imsize by calling casaImagingRoutines.estimate_cell_and_imsize().
        """
        imr.estimate_cell_and_imsize() # call the casaImagingRoutines (imr) estimate_cell_and_imsize function.
        # and some overrides (same as in old phangsPipeline.py)
        if in_file in override_dict:
            if 'cell_size' in override_dict:
                cell_size_string = override_dict[this_vis]['cell_size']
            if 'x_size' in override_dict:
                x_size_string = override_dict[this_vis]['x_size']
            if 'y_size' in override_dict:
                y_size_string = override_dict[this_vis]['y_size']    
        pass
    
    # was phangsPipeline.py buildPhangsCleanCall function
    def buildPhangsCleanCall(self, ):
        """Build a clean call before running the clean task with the function phangsImagingRecipe().
        """
        # 
        # create clean call object
        clean_call = imr.cleanCall() # cleanCall is a class in casaImagingRoutines (imr)
        # calculates cell and imsize
        self.pick_phangs_cell_and_imsize()
        # set specmode and restfreq
        if product == ...:
            clean_call.specmode = ...
            clean_call.restfreq_ghz = line_list.line_list[...]
        # set angular scales to be used in multiscale clean
        if array == ...:
            clean_call.pblimit = ...
            clean_call.smallscalebias = ...
            clean_call.scales_as_angle = ...
        # set overrides
        clean_call.set_overrides(self._kh.get_overrides_for_taget()) # this can be a new function
        # set cleanmask
        clean_call.set_cleanmask(self._kh.get_cleanmask_for_taget()) # this can be a new function # it updates clean_call.clean_mask_file
        # return
        return clean_call
    
    # was phangsPipeline.py phangsImagingRecipe function
    def phangsImagingRecipe(self, ):
        """PHANGS imaging recipe. 
        Including steps:
            - make dirty image.
            - mask alignment.
            - lightly masked multiscale clean.
            - heavily masked single scale clean.
            - export fits.
        This calls the functions in casaImagingRoutines (imr) and casaMaskingRoutines (msr)
        """
        # 
        imr.make_dirty_map(clean_call)
        # 
        imr.replace_cube_with_copy('', '_dirty')
        # 
        if read_in_clean_mask:
            msr.import_and_align_mask()
            clean_call.usemask = 'user'
        else:
            clean_call.usemask = 'pb'
        # 
        if run_multiscale_clean:
            imr.multiscale_loop()
        # 
        if revert_to_multiscale:
            imr.replace_cube_with_copy('', '_multiscale')
        # 
        if make_singlescale_mask:
            msr.signal_mask()
            clean_call.usemask='user'
        # 
        if run_singlescale_clean:
            imr.singlescale_loop()
        # 
        if do_export_to_fits:
            imr.export_to_fits('')
            imr.export_to_fits('_dirty')
            imr.export_to_fits('_multiscale')
        # 
        return


    
    
    
    
    
    







