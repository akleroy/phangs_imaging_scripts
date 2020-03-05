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
        do_signalmask_moment_maps=True, 
        do_hybridmask_moment_maps=True, 
        make_directories=True, 
        extra_ext_in='', 
        extra_ext_out='', 
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
            self._kh.make_missing_directories(derived = True)

        # Loop over target, product, config combinations

        for this_target, this_product, this_config in \
            self.looper(do_targets=True,do_products=True,do_configs=True):
            
            # do signalmask moment maps for each resolution cube
            if do_signalmask_moment_maps:
                for this_res in self._kh.get_res_for_config(this_config):
                    self.task_generate_moment_maps(target=this_target, product=this_product, config=this_config, res=this_res, extra_ext_in=extra_ext_in, extra_ext_out=extra_ext_out)
            
            # do hybridmask moment maps for each resolution cube, using a cube close to 10.72 arcsec resolution
            if do_hybridmask_moment_maps:
                lowres, lowres_tag = self._find_lowest_res(target=this_target, config=this_config, product=this_product, closeto='10.72arcsec')
                for this_res in self._kh.get_res_for_config(this_config):
                    self.task_hybridize_masks(target=this_target, product=this_product, config=this_config, res=this_res, lowres=lowres, extra_ext_in=extra_ext_in, extra_ext_out=extra_ext_out)
            
            # end of loop


    ###########################################
    # Defined file names for various products #
    ###########################################
    
    def _fname_dict(
        self,
        target=None,
        config=None,
        product=None,
        extra_ext_in='',
        extra_ext_out='', 
        res=None,
        res_lowresmask='10p72', 
        ):
        """
        Function to define file names used in other functions.
        
        Output name dict contains following keys: 'in_cube', 'signalmask', 'hybridmask', 'derived_broad_map' and 'derived_strict_map'. 
        
        If input res == 'lowres', then we will find the lowest resolution (closest to res_lowresmask = 10.72 arcsec) 'in_cube'.
        """

        # Error checking
        if target is None:
            raise Exception("Need a target.")
        if product is None:
            raise Exception("Need a product.")
        if config is None:
            raise Exception("Need a config.")
        if res is None:
            raise Exception("Need a res.")
        
        # The output is a nested dictionary structure, for each cube
        # resolution (res_tag)
        fname_dict = {}
        indir = self._kh.get_postprocess_dir_for_target(target=target, changeto=False)
        indir = os.path.abspath(indir)
        #logger.debug('indir: '+indir)
        outdir = self._kh.get_derived_dir_for_target(target=target, changeto=False)
        outdir = os.path.abspath(outdir)
        #logger.debug('outdir: '+outdir)
        if not os.path.isdir(outdir):
            os.makedirs(outdir)
        
        res_tag = None
        if type(res) is str:
            if res == 'lowres':
                res, res_tag = self._find_lowest_res(target=target, config=config, product=product, closeto=res_lowresmask)
        if res_tag is None:
            res_tag = utilsResolutions.get_tag_for_res(res) # this will be like either '5p00' or '80pc'
        
        fname_dict['res'] = res
        fname_dict['res_tag'] = res_tag
        
        cube_filename = utilsFilenames.get_cube_filename(target = target, 
                                                         config = config, 
                                                         product = product,
                                                         ext = 'pbcorr_trimmed_k'+extra_ext_in+'_res'+res_tag,
                                                         casa = False)
        fname_dict['in_cube'] = os.path.join(indir, cube_filename)
        for tag in ['signalmask', 'hybridmask']:
            cube_filename = utilsFilenames.get_cube_filename(target = target, 
                                                             config = config, 
                                                             product = product,
                                                             ext = 'pbcorr_trimmed_k'+extra_ext_in+'_res'+res_tag+extra_ext_out+'_'+tag,
                                                             casa = False)
            fname_dict[tag] = os.path.join(outdir, cube_filename)
        for tag in ['broad', 'strict']:
            derived_name = utilsFilenames.get_derived_rootname(target=target,
                                                               config=config,
                                                               product=product,
                                                               ext='pbcorr_trimmed_k'+extra_ext_in,
                                                               res_tag=res_tag+extra_ext_out,
                                                               derived=tag)
            fname_dict['derived_%s_map'%(tag)] = os.path.join(outdir, derived_name)
        
        return fname_dict
    
    
    def _find_lowest_res(
        self, 
        target = None, 
        config = None, 
        product = None, 
        closeto = None, 
        extra_ext_in = '', 
        ):
        """
        Find the lowest resolution cube data, or closest to the 'closeto' resolution if set.
        """
        indir = self._kh.get_postprocess_dir_for_target(target=target, changeto=False)
        indir = os.path.abspath(indir)
        res_list = self._kh.get_res_for_config(config)
        if res_list is None:
            logger.error('No target resolutions found for target '+target+' and config'+config)
            raise Exception('No target resolutions found for target '+target+' and config'+config)
        # 
        if closeto is not None:
            closest_res_tag = '' # we will find the cloest-resolution cube for the given 'closeto' resolution.
            closest_res = None
        else:
            lowest_res_tag = '' # we will find the lowest-resolution cube for the given 'closeto' resolution.
            lowest_res = None
        # 
        if utilsResolutions.is_physical_resolution(this_res):
            distance = self._kh.get_distance_for_target(target=target)
        elif closest_res is not None and utilsResolutions.is_physical_resolution(closest_res):
            distance = self._kh.get_distance_for_target(target=target)
        else:
            distance = None
        # 
        for this_res in res_list:
            res_tag = utilsResolutions.get_tag_for_res(this_res) # this will be like either '5p00' or '80pc'
            cube_filename = utilsFilenames.get_cube_filename(target = target, 
                                                             config = config, 
                                                             product = product,
                                                             ext = 'pbcorr_trimmed_k'+extra_ext_in+'_res'+res_tag,
                                                             casa = False)
            if os.path.isfile(os.path.join(indir, cube_filename)):
                # 
                if closeto is not None:
                    if closest_res is None:
                        closest_res = this_res
                        closest_res_tag = res_tag
                    elif np.abs(utilsResolutions.get_angular_resolution_for_res(this_res, distance=distance) - utilsResolutions.get_angular_resolution_for_res(closeto, distance=distance)) < \
                         np.abs(utilsResolutions.get_angular_resolution_for_res(closest_res, distance=distance) - utilsResolutions.get_angular_resolution_for_res(closeto, distance=distance)):
                        closest_res = this_res
                        closest_res_tag = res_tag
                else:
                    if lowest_res is None:
                        lowest_res = this_res
                        lowest_res_tag = res_tag
                    elif utilsResolutions.get_angular_resolution_for_res(this_res) > utilsResolutions.get_angular_resolution_for_res(lowest_res):
                        lowest_res = this_res
                        lowest_res_tag = res_tag
            else:
                # file not found
                logger.warning('Cube with resolution tag '+res_tag+' was not found: "'+os.path.join(indir, cube_filename)+'"')
        
        if closeto is not None:
            return closest_res, closest_res_tag
        else:
            return lowest_res, lowest_res_tag
    
    
    ##################################################################
    # Tasks - discrete steps on target, product, config combinations #
    ##################################################################

    def task_generate_moment_maps(
        self,
        target = None, 
        config = None, 
        product = None, 
        res = None, 
        extra_ext_in = '', 
        extra_ext_out = '', 
        ):
        """
        Placeholder for a task to generate a noise cube.
        """
        fname_dict = self._fname_dict(target=target, config=config, product=product, res=res, extra_ext_in=extra_ext_in, extra_ext_out=extra_ext_out)
        if not os.path.isfile(fname_dict['in_cube']):
            logger.warning('Input cube file not found: "'+fname_dict['in_cube']+'"')
            return
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
        ):
        """
        Hybridize each res mask with lowres mask.
        """
        lowres_fname_dict = self._fname_dict(target=target, config=config, product=product, res=lowres, extra_ext_in=extra_ext_in, extra_ext_out=extra_ext_out)
        lowres_mask = fits.getdata(lowres_fname_dict['signalmask'])
        # 
        fname_dict = self._fname_dict(target=target, config=config, product=product, res=res, extra_ext_in=extra_ext_in, extra_ext_out=extra_ext_out)
        mask, header = fits.getdata(fname_dict['signalmask'], header=True)
        # 
        hybridmask = np.logical_or(mask.astype(bool), lowres_mask.astype(bool))
        header['HISTORY'] = ''
        header['HISTORY'] = 'Hybridizing masks "%s" and "%s".'%(lowres_fname_dict['signalmask'], fname_dict['signalmask'])
        header['HISTORY'] = ''
        hybridmask_spectralcube = SpectralCube(data=mask.astype(int), wcs=WCS(header))
        hybridmask_spectralcube.write(fname_dict['hybridmask'], overwrite=True)
        # 
        moment_generator(fname_dict['in_cube'], 
                         root_name = fname_dict['derived_broad_map'], 
                         generate_mask = False, 
                         mask = fname_dict['hybridmask'], 
                         generate_noise = True, 
                         )
                         # mask will have a file name: fname_dict['derived_strict_map']+'_signalmask.fits'
                         # noise will have a file name: fname_dict['derived_strict_map']+'_noise.fits'
        
        


        





