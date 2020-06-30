#!/usr/bin/env python
# 

from __future__ import print_function
import os, sys, re, shutil
import glob
import numpy as np
import scipy.ndimage.morphology as morph
import scipy.ndimage as nd

import logging
logging.basicConfig()
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

casa_enabled = (sys.argv[0].endswith('start_casa.py'))

if casa_enabled:
    logger.debug('casa_enabled = True')
else:
    logger.debug('casa_enabled = False')
    sys.path.append(os.path.dirname(os.path.abspath(__file__)))
    sys.path.append(os.path.dirname(os.path.abspath(__file__))+os.sep+'phangsPipeline')

#print(sys.path)
from phangsPipeline import handlerTemplate
from phangsPipeline import handlerKeys
from phangsPipeline import utilsFilenames


class MakeCleanResidualRegressionHandler(handlerTemplate.HandlerTemplate):
    """Make clean residual regression
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
    
    
    def go(
        self, 
        ):
        """main function to execute.
        """
        for this_target, this_product, this_config in \
            self.looper(do_targets=True, do_products=True, do_configs=True, just_interf=True):
            # 
            self.recipe_make_clean_residual_regression(\
                target=this_target, 
                product=this_product, 
                config=this_config,
                casaext='_multiscale', 
                )
            # 
            self.recipe_make_clean_residual_regression(\
                target=this_target, 
                product=this_product, 
                config=this_config,
                casaext='_singlescale', 
                )
    
    def _fname_dict(
        self,
        product, 
        imagename, 
        ):
        """define file name dict for analysis in this code.
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
    
    
    def recipe_make_clean_residual_regression(
        self,
        target=None, 
        product=None, 
        config=None, 
        casaext='', 
        ):
        """recipe to make clean residual regression for one target.
        """
        if target is None:
            raise Exception('Please input a target!')
        if product is None:
            raise Exception('Please input a product!')
        if config is None:
            raise Exception('Please input a config!')
        # 
        logger.info('------------------------------------------------------------')
        logger.info('Processing target %s product %s config %s'%(target, product, config))
        logger.info('------------------------------------------------------------')
        # 
        image_root = utilsFilenames.get_cube_filename(\
            target=target, product=product, config=config,
            ext=None, casa=True, casaext=casaext)
        # 
        fname_dict = self._fname_dict(\
            product = product, 
            imagename = image_root)
        # 
        image_file = fname_dict['image']
        residual_file = fname_dict['residual']
        mask_file = fname_dict['mask']
        # 
        target_whole_name = self._kh.get_mosaic_target_for_parts(target)
        if target_whole_name is None:
            target_whole_name = target
        # 
        imaging_dir = self._kh.get_imaging_dir_for_target(target_whole_name, changeto=False)
        # 
        image_file = os.path.join(imaging_dir, image_file)
        residual_file = os.path.join(imaging_dir, residual_file)
        mask_file = os.path.join(imaging_dir, mask_file)
        # 
        self.task_make_clean_residual_regression(\
            image_file=image_file, 
            residual_file=residual_file, 
            mask_file=mask_file)
        

    def task_make_clean_residual_regression(
        self, 
        image_file, 
        residual_file, 
        mask_file, 
        ):
        """calculates the clean residual statistics
        """
        logger.info('image    : %r'%(image_file))
        logger.info('residual : %r'%(residual_file))
        logger.info('mask     : %r'%(mask_file))
        
        # TODO
        # exportfits?
        # analyze residual?


########
# main #
########
if __name__ == '__main__':
    
    master_key = 'keys/master_key.txt'
    #master_key = '/Users/dzliu/Work/AlmaPhangs/Works/20200630_PHANGS_ALMA_clean_records/test_phangs_working_dir/keys/master_key.txt'
    if not os.path.isfile(master_key):
        if sys.version_info.major <= 2:
            master_key = raw_input("Please input your master key file path: ")
        else:
            master_key = input("Please input your master key file path: ")
    
    if master_key.find("'") >= 0:
        master_key = master_key.replace("'", "")
    if master_key.find('"') >= 0:
        master_key = master_key.replace('"', '')
    
    if not os.path.isfile(master_key):
        raise Exception('Error! Input master key file does not exist: %r'%(master_key))
    
    this_key_handler = handlerKeys.KeyHandler(master_key = master_key)
    this_handler = MakeCleanResidualRegressionHandler(key_handler = this_key_handler)
    this_handler.set_line_products(only=['co21'])
    this_handler.set_no_cont_products(True)
    this_handler.go()








