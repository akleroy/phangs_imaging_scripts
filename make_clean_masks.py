#!/usr/bin/env python
# 

from __future__ import print_function
import os, sys, re, shutil
import glob
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales, proj_plane_pixel_area
from astropy.convolution import convolve, Gaussian2DKernel, Box1DKernel
import astropy.units as u
import scipy.ndimage.morphology as morph
import scipy.ndimage as nd
import spectral_cube
from spectral_cube import SpectralCube
from radio_beam import Beam

import logging
logging.basicConfig()
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

casa_enabled = (sys.argv[0].endswith('start_casa.py'))

if casa_enabled:
    logger.debug('casa_enabled = True')
    reload(scmasking)
    reload(scproduct)
else:
    logger.debug('casa_enabled = False')
    sys.path.append(os.path.dirname(os.path.abspath(__file__)))
    sys.path.append(os.path.dirname(os.path.abspath(__file__))+os.sep+'phangsPipeline')

#print(sys.path)
from phangsPipeline import scMaskingRoutines as scmasking
from phangsPipeline import scDerivativeRoutines as scproduct
from phangsPipeline import handlerTemplate
from phangsPipeline import handlerKeys
from phangsPipeline import utilsFilenames


class MakeCleanMasksHandler(handlerTemplate.HandlerTemplate):
    """Make clean masks for next round of cleaning processes.
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
        handlerTemplate.HandlerTemplate.__init__(self, key_handler = key_handler, dry_run = dry_run)
        
        # constants
        self.lowres_kernel_size = {} # the 'kernel_lowres' in IDL code "make_clean_masks.pro"
        self.lowres_kernel_size['default'] = 33.0 * u.arcsec
        self.lowres_kernel_size['ngc3599'] = 20.0 * u.arcsec
        self.lowres_kernel_size['ngc4596'] = 20.0 * u.arcsec
        
        self.smoothing_channel_width = {} # the 'chan_lores' in IDL code "make_clean_masks.pro"
        self.smoothing_channel_width['default'] = 20.0
        self.smoothing_channel_width['ngc4596'] = 20.0
        
        self.using_edge_channels = {}
        self.using_edge_channels['default'] = -1
        use_edge_chan = ['ngc3489','ngc3599','ngc3489','ngc4459','ngc4476','ngc4477','ngc4596','ngc7743']
        for t in use_edge_chan:
            self.using_edge_channels[t] = 50
        
        self.masking_thresholds = {}
        self.masking_thresholds['default'] = (3.0, 10.0)
        is_super_faint = ['ngc0300', 'ngc3489', 'ngc3599', 'ngc3489', 'ngc4459', 'ngc4476', 'ngc4477', 'ngc4596', 'ngc7743']
        for t in is_super_faint:
            self.masking_thresholds[t] = (3.0, 5.0)
        
        self.unmasking_edges = {} # we will set these edge pixels to mask=False
        self.unmasking_edges['default'] = None
        self.unmasking_edges['ngc0300'] = {'left':24, 'right':150}
        
        self.masking_dilations = {}
        self.masking_dilations['default'] = {'xy':None, 'v':4}
        self.masking_dilations['ngc0300'] = {'xy':20, 'v':4}
        
        self.masking_bright_center = {}
        self.masking_bright_center['default'] = 0.0*u.arcsec
        has_bright_center = \
            ['circinus', 'ngc0253', 'ngc1097', 'ngc1300', 'ngc1317', 
             'ngc1365', 'ngc1433', 'ngc1512', 'ngc1566', 'ngc1637', 
             'ngc1672', 'ngc2566', 'ngc2903', 'ngc2997', 'ngc3351', 
             'ngc3507', 'ngc3626', 'ngc3627', 'ngc4293', 'ngc4303', 
             'ngc4321', 'ngc4457', 'ngc4535', 'ngc4536', 'ngc4548', 
             'ngc4569', 'ngc4579', 'ngc4826', 'ngc4941', 'ngc4951', 
             'ngc5128', 'ngc5248', 'ngc5643', 'ngc6300', 'ngc7496' ]
        for t in has_bright_center:
            self.masking_bright_center[t] = 20.*u.arcsec
    
    
    def main_loop(
        self, 
        make_directories = True, 
        ):
        
        if make_directories:
            self._kh.make_missing_directories(release = True)

        for this_target, this_product in \
            self.looper(do_targets=True, do_products=True, do_configs=False):
            # 
            self.recipe_make_clean_mask(target=this_target, product=this_product)
    
    
    def task_read_cube(\
        self, 
        target = None, 
        config = None, 
        product = None, 
        history = None, 
        ):
        """
        We will loop over all configs and read in a cube with ideal resolution.
        """
        
        if target is None:
            raise Exception('Please input a target.')
        if product is None:
            raise Exception('Please input a product.')
        
        # record masking history
        masking_history = [] if history is None else history
        
        cube_root = self._kh.get_derived_dir_for_target(target = target) # get_product_dir_for_target
        cube_root = os.path.abspath(cube_root)+os.sep
        
        # loop all configs and find the most ideal one
        cube = None
        #interf_config = self._kh.get_interf_configs()
        #feather_config = self._kh.get_feather_configs()
        #all_configs = []
        #if interf_config is not None:
        #    all_configs.extend(interf_config)
        #if feather_config is not None:
        #    all_configs.extend(feather_config)
        #all_config_resolutions = [utilsResolutions.get_angular_resolution_for_config(t) for t in all_configs]
        #loop_config_indices = np.array(all_config_resolutions).argsort().tolist()
        all_configs = ['7m+tp','7m','12m','12m+7m+tp']
        for config in all_configs:
            cube_fname = utilsFilenames.get_cube_filename(\
                            target = target, 
                            product = product, 
                            config = config, 
                            ext = 'pbcorr_round_k', 
                            casa = False )
                            # ext = 'flat_round_k', -- <TODO>
            
            #cube_data, cube_header = fits.getdata(\
            #        cube_fname, 
            #        header = True )
            #
            #cube_wcs = WCS(cube_header)
            #
            #pixscale = proj_plane_pixel_scales(cube_wcs)[-1]*u.deg
            # 
            #cube = SpectralCube(data=cube_data, wcs=cube_wcs, header=cube_header)
            
            if os.path.isfile(cube_root+cube_fname):
                logger.info('Reading spectral cube "%s".'%(cube_root+cube_fname))
                cube = SpectralCube.read(cube_root+cube_fname)
                masking_history.append('Reading spectral cube "%s".'%(cube_root+cube_fname))
                break # do not need to read next config
            else:
                logger.info('Spectral cube "%s" not found. Try next.'%(cube_root+cube_fname))
                cube = None
        
        if cube is None:
            raise Exception()
        
        # return the SpectralCube object
        return cube, masking_history
    
    
    def task_make_convolution(\
        self, 
        target = None, 
        cube = None, 
        header = None, 
        pixscale = None, 
        history = None, 
        ):
        
        if target is None:
            raise Exception('Please input a target.')
        if cube is None:
            raise Exception('Please input a cube.')
        
        # record masking history
        masking_history = [] if history is None else history
        
        # target specification
        if target in self.lowres_kernel_size:
            kernel_lores = self.lowres_kernel_size[target]
        else:
            kernel_lores = self.lowres_kernel_size['default']
        
        logger.info('Convolving cube with kernel %s arcsec.'%(kernel_lores.to(u.arcsec).value))
        
        if isinstance(cube, np.ndarray):
            if pixscale is None and header is not None:
                if isinstance(header, fits.Header):
                    pixscale = proj_plane_pixel_scales(WCS(header))[-1] * u.deg
            if pixscale is None:
                raise Exception('Please input a pixscale or fits header.')
            kernel = Gaussian2DKernel((kernel_lores/pixscale).to(u.one).value / np.sqrt(8 * np.log(2)))
            lowres_cube = convolve(cube, kernel)
        elif isinstance(cube, SpectralCube):
            lowres_cube = cube.convolve_to(Beam(major=kernel_lores, minor=kernel_lores, pa=0.0*u.deg))
        else:
            raise Exception('The input cube is not a SpectralCube! Type is %s'%(str(type(cube))))
        
        masking_history.append('Convolving spectral cube with kernel %s arcsec.'%(kernel_lores.to(u.arcsec).value))
        
        return lowres_cube, masking_history
    
    
    def task_make_smoothing(\
        self, 
        target = None, 
        cube = None, 
        history = None, 
        ):
        
        if target is None:
            raise Exception('Please input a target.')
        if cube is None:
            raise Exception('Please input a cube.')
        
        # record masking history
        masking_history = [] if history is None else history
        
        # target specification
        if target in self.smoothing_channel_width:
            chan_lores = self.smoothing_channel_width[target]
        else:
            chan_lores = self.smoothing_channel_width['default']
        
        logger.info('Smoothing spectral cube with kernel %s channels.'%(chan_lores)) #<TODO># km/s?
        
        if isinstance(cube, np.ndarray):
            #smoothed_cube = convolve(cube, Box1DKernel(chan_lores))
            raise Exception('The input cube is not a SpectralCube! Type is %s'%(str(type(cube))))
            raise NotImplementedError() # should convolve only along spectral axis
        elif isinstance(cube, SpectralCube):
            smoothed_cube = cube.spectral_smooth(Box1DKernel(chan_lores))
        else:
            raise Exception('The input cube is not a SpectralCube! Type is %s'%(str(type(cube))))
        
        masking_history.append('Smoothing spectral cube with kernel %s channels.'%(chan_lores))
        
        return smoothed_cube, masking_history
    
    
    def task_calculate_rms(\
        self, 
        target = None, 
        cube = None, 
        history = None, 
        ):
        
        if target is None:
            raise Exception('Please input a target.')
        if cube is None:
            raise Exception('Please input a cube.')
        
        # record masking history
        masking_history = [] if history is None else history
        noise_params = {\
            'spec_box': 5,
            'box':5,
            'iterations':3
            }
        # target specification
        if target in self.using_edge_channels:
            nchan_edge = self.using_edge_channels[target]
        else:
            nchan_edge = self.using_edge_channels['default']
        
        if isinstance(cube, np.ndarray):
            cube_data = cube
        elif isinstance(cube, SpectralCube):
            cube_data = cube._data # cube.unmasked_data() # cube._data
        else:
            raise Exception('The input cube is not a SpectralCube! Type is %s'%(str(type(cube))))
        
        if nchan_edge > 0:
            logger.info('Calculating rms noise cube using only %d edge channels.'%(nchan_edge))
            rms_cube_1 = scmasking.noise_cube(cube_data[0:nchan_edge, :, :], **noise_params)
            rms_cube_2 = scmasking.noise_cube(
                cube_data[cube.shape[0]-nchan_edge:cube.shape[0], :, :], **noise_params)
            rms_cube = 0.5 * (rms_cube_1 + rms_cube_2)
            masking_history.append('Calculating rms noise cube using only %d edge channels.'%(nchan_edge))
        else:
            logger.info('Calculating rms noise cube.')
            rms_cube = scmasking.noise_cube(cube_data, **noise_params)
            masking_history.append('Calculating rms noise cube.')
        
        return rms_cube, masking_history
    
    
    def task_calculate_ppm(\
        self, 
        cube = None, 
        header = None, 
        history = None, 
        ):
        
        if cube is None:
            raise Exception('Please input a cube.')
        
        # record masking history
        masking_history = [] if history is None else history
        
        if isinstance(cube, np.ndarray):
            if not isinstance(header, fits.Header):
                raise Exception('Please input a fits header.')
            if 'BMAJ' in header and 'BMIN' in header and 'BPA' in header:
                #ppbeam = ( Beam(major=header['BMAJ']*u.deg, minor=header['BMAJ']*u.deg, pa=header['BPA']*u.deg).sr.value \
                #           / proj_plane_pixel_area(WCS(header)) \
                #         ).to(u.one).value # calculate pixel per beam #<TODO># NOT TESTED! # https://github.com/radio-astro-tools/spectral-cube/blob/master/spectral_cube/base_class.py
                ppbeam = np.pi/(4.0*np.log(2.0))*((header['BMAJ'])*(header['BMIN'])) / proj_plane_pixel_area(WCS(header))
            else:
                raise Exception('The input fits header does not contain BMAJ, BMIN or BPA!')
        elif isinstance(cube, SpectralCube):
            ppbeam = cube.pixels_per_beam # calculate pixel per beam 
            # tested: returns same as above.
        else:
            raise Exception('The input cube is not a SpectralCube! Type is %s'%(str(type(cube))))
        
        logger.info('Calculating pixels per beam: %s'%(ppbeam))
        masking_history.append('Calculating pixels per beam: %s'%(ppbeam))
        
        return ppbeam, masking_history
    
    
    def task_make_mask(\
        self, 
        target = None, 
        cube = None, 
        rms = None, 
        header = None, 
        history = None, 
        ):
        
        if target is None:
            raise Exception('Please input a target.')
        if cube is None:
            raise Exception('Please input a cube.')
        if rms is None:
            raise Exception('Please input a rms cube.')
        
        # record masking history
        masking_history = [] if history is None else history
        
        if target in self.masking_thresholds:
            lo_thresh_lowres, hi_thresh_lowres = self.masking_thresholds[target]
        else:
            lo_thresh_lowres, hi_thresh_lowres = self.masking_thresholds['default']
        
        #if target in self.masking_dilations:
        #    grow_xy, grow_v = self.masking_dilations[target]
        #else:
        #    grow_xy, grow_v = self.masking_dilations['default']
        
        if isinstance(cube, np.ndarray):
            cube_data = cube
            lowres_ppbeam, masking_history = self.task_calculate_ppm(cube_data, header = header, history = masking_history) # calculate pixel per beam #<TODO># NOT TESTED!
        elif isinstance(cube, SpectralCube):
            cube_data = cube._data
            lowres_ppbeam, masking_history = self.task_calculate_ppm(cube_data, header = cube.header, history = masking_history) # calculate pixel per beam #<TODO># NOT TESTED!
        else:
            raise Exception('The input cube is not a SpectralCube! Type is %s'%(str(type(cube))))
        
        # make cprops mask
        masking_params = {\
            'hi_thresh': hi_thresh_lowres, 'hi_nchan': 2,
            'lo_thresh': lo_thresh_lowres, 'lo_nchan': 2,
            'min_pix': 4.0*lowres_ppbeam, 
            'min_area': 2.0*lowres_ppbeam,
            }
            # 'grow_xy': grow_xy, 
            # 'grow_v': grow_v, 
            #--> we will do dilation afterwards because some galaxies need unmasking edge before dilation.
        masking_params_str = '('+', '.join("{!s}={!r}".format(k, masking_params[k]) for k in masking_params.keys())+')'
        logger.info('Making mask with params: '+masking_params_str)
        
        mask = scmasking.simple_mask(cube_data, rms, **masking_params)
        
        masking_history.append('Masking with scMaskingRoutines simple_mask function with params:')
        masking_history.extend(["... {!s} = {!r}".format(k, masking_params[k]) for k in masking_params.keys()])
        
        return mask, masking_history
    
    
    def task_grow_mask(
        self, 
        target = None, 
        mask = None, 
        grow_xy = None, 
        grow_v = None, 
        header = None, 
        pixscale = None, 
        history = None, 
        ):
        
        if target is None:
            raise Exception('Please input a target.')
        if mask is None:
            raise Exception('Please input a mask.')
        
        if isinstance(mask, np.ndarray):
            pass
        elif isinstance(mask, SpectralCube):
            pass
        else:
            raise Exception('The input mask is not a numpy ndarray! Type is %s'%(str(type(mask))))
        
        # record masking history
        masking_history = [] if history is None else history
        
        # unmasking image edges
        if target in self.unmasking_edges:
            if 'left' in self.unmasking_edges[target]:
                nchan_edge = self.unmasking_edges[target]['left']
                mask[:, :, 0:nchan_edge+1] = False
                masking_history.append('Unmasking left edge by %d.'%(nchan_edge))
            if 'right' in self.unmasking_edges[target]:
                nchan_edge = self.unmasking_edges[target]['right']
                mask[:, :, nchan_edge:] = False
                masking_history.append('Unmasking right edge by %d.'%(nchan_edge))
            if 'top' in self.unmasking_edges[target]:
                nchan_edge = self.unmasking_edges[target]['top']
                mask[:, 0:nchan_edge+1, :] = False
                masking_history.append('Unmasking top edge by %d.'%(nchan_edge))
            if 'bottom' in self.unmasking_edges[target]:
                nchan_edge = self.unmasking_edges[target]['bottom']
                mask[:, nchan_edge:, :] = False
                masking_history.append('Unmasking bottom edge by %d.'%(nchan_edge))
        
        # dilation
        if target in self.masking_dilations:
            if 'xy' in self.masking_dilations[target]:
                grow_xy = self.masking_dilations[target]['xy']
            else:
                grow_xy = None
            if 'v' in self.masking_dilations[target]:
                grow_v = self.masking_dilations[target]['v']
            else:
                grow_v = None
        else:
            if 'xy' in self.masking_dilations['default']:
                grow_xy = self.masking_dilations['default']['xy']
            else:
                grow_xy = None
            if 'v' in self.masking_dilations['default']:
                grow_v = self.masking_dilations['default']['v']
            else:
                grow_v = None
        
        logger.info('Growing mask.')
        
        # grow xy --> see inside scmasking.simple_mask
        if grow_xy is not None:
            struct = morph.iterate_structure(morph.generate_binary_structure(2, 1), grow_xy)
            mask = morph.binary_dilation(mask, struct[np.newaxis, :, :])
            masking_history.append('Growing mask by %d in spatial axes.'%(grow_xy))
        
        # grow v --> see inside scmasking.simple_mask
        if grow_v is not None:
            struct = np.ones(grow_v, dtype=np.bool)
            mask = morph.binary_dilation(mask, struct[:, np.newaxis, np.newaxis])
            masking_history.append('Growing mask by %d in spectral axis.'%(grow_v))
        
        # mask bright center
        if target in self.masking_bright_center:
            center_radius = self.masking_bright_center[target]
            if center_radius > 0.0*u.arcsec:
                if pixscale is None and header is not None:
                    if isinstance(header, fits.Header):
                        pixscale = proj_plane_pixel_scales(WCS(header))[-1] * u.deg
                if pixscale is None:
                    raise Exception('Please input a pixscale or header.')
                gridy, gridx = np.mgrid[0:mask.shape[1], 0:mask.shape[2]]
                gridy = gridy - ((mask.shape[1]-1)/2.0)
                gridx = gridx - ((mask.shape[2]-1)/2.0)
                rad = np.sqrt(gridx*gridx+gridy*gridy) * pixscale.to(u.arcsec)
                center_mask = (rad <= center_radius) # arcsec
                mask = np.logical_or(mask, center_mask)
                masking_history.append('Masking bright center within %s arcsec.'%(str(center_radius.to(u.arcsec).value)))
        
        return mask, masking_history
        
    
    def task_write_fits(
        self, 
        mask = None, 
        history = None, 
        header = None, 
        wcs = None, 
        outfile = None, 
        ):
        """
        Write a cube mask 3D array to disk.
        """
        if header is not None and wcs is None:
            wcs = WCS(header)
        if mask is not None and wcs is not None and outfile is not None:
            # check existing output file
            if os.path.isfile(outfile):
                os.remove(outfile)
            elif not os.path.isdir(os.path.dirname(outfile)):
                os.makedirs(os.path.dirname(outfile))
            # prepare fits header
            header_to_write = wcs.to_header().copy()
            for key in ['MJD-OBS']:
                if key in header_to_write:
                    del header_to_write[key] #<TODO># there are two 'MJD-OBS' in header?
            if header is not None:
                for key in ['BMAJ', 'BMIN', 'BPA', 'MJD-OBS']:
                #for key in header.keys():
                    if key in header and key not in header_to_write:
                        header_to_write[key] = header[key] # copy input header keywords into output header_to_write
            # write masking history into fits header
            if history is not None:
                if np.isscalar(history):
                    history = [history]
                header_to_write['HISTORY'] = ''
                for item in history:
                    header_to_write['HISTORY'] = str(item).strip()
                header_to_write['HISTORY'] = ''
            #print(header_to_write)
            mask_spectralcube = SpectralCube(data=mask.astype(int), wcs=wcs, header=header_to_write)
            if os.path.isfile(outfile):
                os.remove(outfile)
            mask_spectralcube.write(outfile, format="fits")
            logger.info('Output to "%s".'%(outfile))
    
    
    def recipe_make_clean_mask(
        self, 
        target = None, 
        product = None, 
        config = None,
        ):
        
        # read cube
        cube, masking_history = self.task_read_cube(target=target,
                                                    product=product,
                                                    config=config)
        
        # make convolution
        lowres_cube, masking_history = self.task_make_convolution(target=target, cube=cube, history=masking_history)
        
        # boxcar smoothing along spectral axis
        lowres_cube, masking_history = self.task_make_smoothing(target=target, cube=lowres_cube, history=masking_history)
        
        # calculate rms noise cube
        rms_lowres_cube, masking_history = self.task_calculate_rms(target=target, cube=lowres_cube, history=masking_history)
        
        # make cprops mask
        mask, masking_history = self.task_make_mask(target=target, cube=lowres_cube, rms=rms_lowres_cube, history=masking_history)
        
        # grow mask and deal with edge and bright center
        mask, masking_history = self.task_grow_mask(target=target, mask=mask, history=masking_history)
        
        # write fits file
        #cube_file_name = re.sub(r'^Reading spectral cube \"(.*)\"\.$', r'\1', masking_history[0])
        #cube_header = fits.getheader(cube_file_name)
        #clean_mask_fname = self._kh.get_cleanmask_filename(target=target, product=product) # target+'_'+product+'_clean_mask.fits'
        clean_mask_fname = (self._kh.get_cleanmask_dir_for_target(target)
                            + target+'_'+product+'_clean_mask.fits')
        self.task_write_fits(mask=mask, wcs=lowres_cube.wcs, header=cube.header, outfile=clean_mask_fname, history=masking_history)



# 
# main
# # 
# if __name__ == '__main__':
#     this_key_handler = handlerKeys.KeyHandler()
#     this_handler = MakeCleanMasksHandler(key_handler = this_key_handler)
#     this_handler.set_no_interf_configs(True)
#     this_handler.set_line_products(only=['co21'])
#     this_handler.set_no_cont_products(True)
#     this_handler.main_loop()








