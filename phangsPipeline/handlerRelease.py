"""releaseHandler
TBE.
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

import logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

casa_enabled = (sys.argv[0].endswith('start_casa.py'))

if casa_enabled:
    logger.debug('casa_enabled = True')
else:
    logger.debug('casa_enabled = False')
    sys.path.append(os.path.dirname(os.path.abspath(__file__)))

import utils
import line_list
import handlerTemplate


class ReleaseHandler(handlerTemplate.HandlerTemplate):
	"""
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



	###########################################
	# Defined file names for various products #
	###########################################

	def _fname_dict(
	    self,
	    target=None,
	    config=None,
	    product=None,
	    ):

	    if target is None:
	        logger.error("Need a target.")
	        return()

	    if product is None:
	        logger.error("Need a product.")
	        return()

	    if config is None:
	        logger.error("Need a config.")
	        return()
	    
	    fname_dict = {} # here we will create a subdict for each cube resolution (res_tag)
	    
	    indir_postprocess = self._kh.get_postprocess_dir_for_target(target=target, changeto=False)
	    indir_postprocess = os.path.abspath(indir_postprocess)
	    logger.debug('indir_postprocess: '+indir_postprocess)
	    
	    indir_product = self._kh.get_product_dir_for_target(target=target, changeto=False)
	    indir_product = os.path.abspath(indir_product)
	    logger.debug('indir_product: '+indir_product)

	    outdir = self._kh.get_release_dir_for_target(target=target, changeto=False)
	    outdir = os.path.abspath(outdir)
	    logger.debug('outdir: '+outdir)
	    
	    if not os.path.isdir(outdir):
	        os.makedirs(outdir)
	    
	    res_list = self._kh.get_res_for_config(config)
	    if res_list is None:
	        logger.error('No target resolutions found for target '+target+' and config'+config)
	        raise Exception('No target resolutions found for target '+target+' and config'+config)
	    
	    cube_filename = self._kh.get_cube_filename(target = target, 
	                                               config = config, 
	                                               product = product,
	                                               ext = 'pbcorr_trimmed_k',
	                                               casa = False)

	    fname_dict['cube'] = {}
	    if os.path.isfile(os.path.join(indir_postprocess, cube_filename)):
	        fname_dict['cube']['pbcorr_trimmed_k'] = os.path.join(indir_postprocess, cube_filename)
	        hybridmask_filename = self._kh.get_cube_filename(target = target, config = config, product = product,ext = 'hybridmask', casa = False)
	        if os.path.isfile(os.path.join(indir_product, hybridmask_filename)):
	            fname_dict['cube']['hybridmask'] = os.path.join(indir_product, hybridmask_filename)

	        signalmask_filename = self._kh.get_cube_filename(target = target, config = config, product = product, ext = 'signalmask', casa = False)

	        if os.path.isfile(os.path.join(indir_product, signalmask_filename)):
	            fname_dict['cube']['signalmask'] = os.path.join(indir_product, signalmask_filename)

	    for this_res in res_list:
	        res_tag = self._kh.get_tag_for_res(this_res)
	        for this_mask in ["broad", "strict"]:
	            for this_mom in ["mom0", "mom1", "mom2", "eq", "tmax", "vmax", "vquad"]:
	                cube_filename = self._kh.get_cube_filename(target = target, config = config, product = product,
	                                                           ext = this_mask+"_"+this_mom+"_res"+res_tag,
	                                                           casa = False)

                        #if os.path.isfile(os.path.join(indir, cube_filename)):

            return fname_dict

	    """
	    for this_res in res_list:
	        
	        res_tag = self._kh.get_tag_for_res(this_res)
	        cube_filename = self._kh.get_cube_filename(target = target, 
	                                                   config = config, 
	                                                   product = product,
	                                                   ext = 'pbcorr_trimmed_k'+'_res'+res_tag,
	                                                   casa = False)
	        
	        if os.path.isfile(os.path.join(indir, cube_filename)):
	            fname_dict[res_tag] = {}
	            fname_dict[res_tag]['pbcorr_trimmed_k'] = os.path.join(indir, cube_filename)
	            # 
	            image_basename = re.sub('_pbcorr_trimmed_k_res'+res_tag+r'\.fits$', '', os.path.basename(cube_filename)) # remove suffix
	            for tag in ['hybridmask', 'signalmask']:
	                fname_dict[res_tag][tag] = os.path.join(outdir, image_basename+'_'+tag+'_res'+res_tag+'.fits')
	            for tag in ['broad', 'strict']:
	                fname_dict[res_tag][tag] = os.path.join(outdir, image_basename+'_'+tag) # we will append mom0 mom1 then res_tag
	        else:
	            # file not found
	            logger.warning('Cube with tag '+res_tag+' at '+str(np.round(this_res,2))+' arcsec resolution was not found: "'+os.path.join(indir, cube_filename)+'"')
	    """


	########################
	# loop_products_making #
	########################

	def loop_build_release(
	    self,
	    make_directories=True,
	    ):
	    """
	    """
	    
	    # Replacing build_release.pro
	    
	    if len(self.get_targets()) == 0:            
	        logger.error("Need a target list.")
	        return(None)

	    if len(self.get_all_products()) == 0:            
	        logger.error("Need a products list.")
	        return(None)

	    if make_directories:
	        self._kh.make_missing_directories(release = True)

	    for this_target, this_product, this_config in \
	        self.looper(do_targets=True,do_products=True,do_configs=True):
	        fname_dict = self._fname_dict(target = this_target, 
	                                      config = this_config, 
	                                      product = this_product)

	    return fname_dict

	    """
	        # 
	        # convolve to specific resolution (this is now done by postprocessHandler.py)
	        # 
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
	    """

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
	        res_tag = self._kh.get_tag_for_res(this_res)
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
	                self.task_writing_mask(\
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
	                self.task_writing_mask(
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
	                self.task_writting_products(
	                    cube = broadcube_data,
	                    rms = broadcube_noise,
	                    outfile = fname_dict[res_tag]['broad'], 
	                    res_tag = res_tag, 
	                    )
	                # 
	                # make moment maps and other products with the signal (strict) mask
	                self.task_writting_products(
	                    cube = strictcube_data,
	                    rms = strictcube_noise,
	                    outfile = fname_dict[res_tag]['strict'], 
	                    res_tag = res_tag, 
	                    )
	        # 



	#####################
	# task_writing_mask #
	#####################

	def task_writing_mask(
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


	##########################
	# task_writting_products #
	##########################

	def task_writting_products(
	    self,
	    cube,
	    rms,
	    outfile, 
	    res_tag=None,
	    do_moment0 = True,
	    do_moment1 = True,
	    do_moment2 = True,
	    do_ew = True,
	    do_tmax = True,
	    do_vmax = True,
	    do_vquad = True,
	    do_errormaps = True,
	    ):
	    """Collapse masked cube as our products and write to disk. 
	    """
	    # 
	    if re.match(r'\.fits$', outfile):
	        outfile = re.sub(r'\.fits$', r'', outfile) # remove suffix
	    # 
	    if res_tag is not None:
	        res_tag_ext = '_res'+res_tag
	    else:
	        res_tag_ext = ''
	    # 
	    process_list = []
	    if do_moment0:   
	        process_list.append({'outfile': outfile+'_mom0'+res_tag_ext+'.fits',  
	                             'errorfile': outfile+'_emom0'+res_tag_ext+'.fits',  
	                             'func': scproduct.write_moment0, 
	                             'unit': u.K * u.km/u.s })
	    if do_moment1:   
	        process_list.append({'outfile': outfile+'_mom1'+res_tag_ext+'.fits',  
	                             'errorfile': outfile+'_emom1'+res_tag_ext+'.fits',  
	                             'func': scproduct.write_moment1, 
	                             'unit': u.km/u.s })
	    if do_moment2:   
	        process_list.append({'outfile': outfile+'_mom2'+res_tag_ext+'.fits',  
	                             'errorfile': outfile+'_emom2'+res_tag_ext+'.fits',  
	                             'func': scproduct.write_moment2, 
	                             'unit': u.km/u.s })
	    if do_ew:        
	        process_list.append({'outfile': outfile+'_ew'+res_tag_ext+'.fits',    
	                             'errorfile': outfile+'_eew'+res_tag_ext+'.fits',    
	                             'func': scproduct.write_ew,      
	                             'unit': u.km/u.s })
	    if do_tmax:      
	        process_list.append({'outfile': outfile+'_tmax'+res_tag_ext+'.fits',  
	                             'errorfile': outfile+'_etmax'+res_tag_ext+'.fits',  
	                             'func': scproduct.write_tmax,    
	                             'unit': u.K })
	    if do_vmax:      
	        process_list.append({'outfile': outfile+'_vmax'+res_tag_ext+'.fits',  
	                             'errorfile': outfile+'_evmax'+res_tag_ext+'.fits',  
	                             'func': scproduct.write_vmax,    
	                             'unit': u.km/u.s })
	    if do_vquad:     
	        process_list.append({'outfile': outfile+'_vquad'+res_tag_ext+'.fits', 
	                             'errorfile': outfile+'_evquad'+res_tag_ext+'.fits', 
	                             'func': scproduct.write_vquad,   
	                             'unit': u.km/u.s })
	    # 
	    # delete old files
	    for process_dict in process_list:
	        outfile = process_dict['outfile']
	        errorfile = process_dict['errorfile'] if do_errormaps else None
	        if os.path.isfile(outfile):
	            os.remove(outfile)
	            logger.debug('Deleting old file "'+outfile+'"')
	        if os.path.isfile(errorfile):
	            os.remove(errorfile)
	            logger.debug('Deleting old file "'+errorfile+'"')
	    # 
	    # make moment maps and other products using scProductRoutines functions.
	    for process_dict in process_list:
	        processfunction = process_dict['func']
	        outfile = process_dict['outfile']
	        errorfile = process_dict['errorfile'] if do_errormaps else None
	        logger.info('Producing "'+outfile+'"')
	        processfunction(cube,
	                        rms = rms, 
	                        outfile = outfile, 
	                        errorfile = errorfile)
	    
