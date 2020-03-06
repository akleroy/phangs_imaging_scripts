
# Test under:
#   CASA 5.4.0
# 
# How to run:
#   casapy-5.4.0 --nogui --log2term -c "exec(open('test_scripts/unittest_handlerImaging.py','r').read())"
# 
# TODO: 
# 
# DONE: 
#   Needs to add 'do_revert_to_singlescale' argument to the 'handlerImaging.loop_imaging' function.
#   Ran 14 tests in 33836.663s
# 

from __future__ import print_function
import os, sys, shutil, logging, unittest
logging.basicConfig()
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
sys.path.append('.')
from phangsPipeline import handlerTemplate
from phangsPipeline import handlerKeys
from phangsPipeline import handlerImaging
from phangsPipeline import utilsFilenames
if sys.version_info.major <= 2:
    reload(handlerTemplate)
    reload(handlerKeys)
    reload(handlerImaging)
    reload(utilsFilenames)



class test_imaging_handler(unittest.TestCase):
    
    key_handler = None
    imaging_handler = None
    current_dir = None
    passed_steps = None
    
    def setUp(self):
        self.current_dir = os.getcwd()
    
    @classmethod
    def setUpClass(cls):
        # requires Python >= 2.7
        cls.key_handler = handlerKeys.KeyHandler(master_key = 'test_keys/master_key.txt', dochecks = False)
        cls.imaging_handler = handlerImaging.ImagingHandler(key_handler = cls.key_handler)
        cls.key_handler.get_vis_filename = utilsFilenames.get_vis_filename
        cls.key_handler.get_cube_filename = utilsFilenames.get_cube_filename
        #self.key_handler._key_dir = os.getcwd()+os.sep+'key_templates'+os.sep
        #self.key_handler._imaging_root = os.getcwd()+os.sep+'test'+os.sep+'imaging'+os.sep
        #self.key_handler._postprocess_root = os.getcwd()+os.sep+'test'+os.sep+'postprocess'+os.sep
        #self.key_handler._product_root = os.getcwd()+os.sep+'test'+os.sep+'product'+os.sep
        #if not os.path.isdir(self.key_handler._imaging_root):
        #    os.makedirs(self.key_handler._imaging_root)
        cls.passed_steps = []
        cls.passed_steps = ['dirty','cleanmask','multiscale','singlescale']
        
    def tearDown(self):
        #if self.key_handler is not None:
        #    del self.key_handler
        #if self.imaging_handler is not None:
        #    del self.imaging_handler
        if self.current_dir is not None:
            os.chdir(self.current_dir)
        #if os.path.isdir('test'):
        #    shutil.rmtree('test')
        pass
    
    def test_initialize_key_handler(self):
        assert self.key_handler is not None
    
    def test_initialize_imaging_handler(self):
        assert self.imaging_handler is not None
    
    def test_do_dirty_imaging(self):
        if 'dirty' in self.passed_steps:
            return
        target = self.key_handler.get_targets()[0]
        product = 'co21' # self.key_handler.get_line_products()[0]
        config = '7m' # self.key_handler.get_interf_configs()[0]
        self.imaging_handler.set_targets(only=[target])
        self.imaging_handler.set_line_products(only=[product])
        self.imaging_handler.set_interf_configs(only=['7m'])
        self.imaging_handler.set_no_feather_configs(True)
        self.imaging_handler.set_no_cont_products(True)
        # 
        self.imaging_handler.loop_imaging(make_directories = True, \
                                          do_dirty_image = True,
                                          do_revert_to_dirty = False,
                                          do_read_clean_mask = False, 
                                          do_multiscale_clean = False,
                                          do_revert_to_multiscale = False,
                                          do_singlescale_mask = False,
                                          do_singlescale_clean = False,
                                          do_revert_to_singlescale = False,
                                          do_export_to_fits = False, 
                                          overwrite = False)
        # 
        # check the existence of dirty image
        #assert utilsFilenames.get_cube_filename(target = target, config = config, product = product, casa = True, ext='dirty') is not None
        #assert os.path.isdir(utilsFilenames.get_cube_filename(target = target, config = config, product = product, casa = True, ext='dirty'))
        assert self.key_handler.get_cube_filename(target = target, config = config, product = product, casa = True, ext='dirty') is not None
        this_cube_filename = self.key_handler.get_cube_filename(target = target, config = config, product = product, casa = True, ext='dirty')
        logger.debug('%r %s %s'%(this_cube_filename, 'isdir?', os.path.isdir(this_cube_filename) ) )
        assert os.path.isdir(this_cube_filename)
    
    def test_do_read_clean_mask(self):
        if 'cleanmask' in self.passed_steps:
            return
        target = self.key_handler.get_targets()[0]
        product = 'co21' # self.key_handler.get_line_products()[0]
        config = '7m' # self.key_handler.get_interf_configs()[0]
        self.imaging_handler.set_targets(only=[target])
        self.imaging_handler.set_line_products(only=[product])
        self.imaging_handler.set_interf_configs(only=['7m'])
        self.imaging_handler.set_no_feather_configs(True)
        self.imaging_handler.set_no_cont_products(True)
        # 
        # before reading clean mask, the mask data should not exist.
        #assert utilsFilenames.get_cube_filename(target = target, config = config, product = product, casa = True, casaext = '.mask') is not None
        #assert os.path.isdir(utilsFilenames.get_cube_filename(target = target, config = config, product = product, casa = True, casaext = '.mask'))
        assert self.key_handler.get_cube_filename(target = target, config = config, product = product, casa = True, casaext = '.mask') is not None
        this_cube_filename = self.key_handler.get_cube_filename(target = target, config = config, product = product, casa = True, casaext = '.mask')
        logger.debug('%r %s %s'%(this_cube_filename, 'isdir?', os.path.isdir(this_cube_filename) ) )
        assert not os.path.isdir(this_cube_filename)
        # 
        self.imaging_handler.loop_imaging(make_directories = True, \
                                          do_dirty_image = False,
                                          do_revert_to_dirty = True,
                                          do_read_clean_mask = True, 
                                          do_multiscale_clean = False,
                                          do_revert_to_multiscale = False,
                                          do_singlescale_mask = False,
                                          do_singlescale_clean = False,
                                          do_revert_to_singlescale = False,
                                          do_export_to_fits = False, 
                                          overwrite = False)
        # 
        # after reading clean mask, the mask data should exist.
        #assert utilsFilenames.get_cube_filename(target = target, config = config, product = product, casa = True, casaext = '.mask') is not None
        #assert os.path.isdir(utilsFilenames.get_cube_filename(target = target, config = config, product = product, casa = True, casaext = '.mask'))
        assert self.key_handler.get_cube_filename(target = target, config = config, product = product, casa = True, casaext = '.mask') is not None
        this_cube_filename = self.key_handler.get_cube_filename(target = target, config = config, product = product, casa = True, casaext = '.mask')
        logger.debug('%r %s %s'%(this_cube_filename, 'isdir?', os.path.isdir(this_cube_filename) ) )
        assert os.path.isdir(this_cube_filename)
    
    def test_do_multiscale_imaging(self):
        if 'multiscale' in self.passed_steps:
            return
        target = self.key_handler.get_targets()[0]
        product = 'co21' # self.key_handler.get_line_products()[0]
        config = '7m' # self.key_handler.get_interf_configs()[0]
        self.imaging_handler.set_targets(only=[target])
        self.imaging_handler.set_line_products(only=[product])
        self.imaging_handler.set_interf_configs(only=['7m'])
        self.imaging_handler.set_no_feather_configs(True)
        self.imaging_handler.set_no_cont_products(True)
        # 
        self.imaging_handler.loop_imaging(make_directories = True, \
                                          do_dirty_image = False,
                                          do_revert_to_dirty = True,
                                          do_read_clean_mask = True, 
                                          do_multiscale_clean = True,
                                          do_revert_to_multiscale = False,
                                          do_singlescale_mask = False,
                                          do_singlescale_clean = False,
                                          do_revert_to_singlescale = False,
                                          do_export_to_fits = False, 
                                          overwrite = False)
        # 
        # check the existence of multiscale-cleaned cube image
        #assert utilsFilenames.get_cube_filename(target = target, config = config, product = product, casa = True, ext = 'multiscale') is not None
        #assert os.path.isdir(utilsFilenames.get_cube_filename(target = target, config = config, product = product, casa = True, ext = 'multiscale'))
        assert self.key_handler.get_cube_filename(target = target, config = config, product = product, casa = True, ext = 'multiscale') is not None
        this_cube_filename = self.key_handler.get_cube_filename(target = target, config = config, product = product, casa = True, ext = 'multiscale')
        logger.debug('%r %s %s'%(this_cube_filename, 'isdir?', os.path.isdir(this_cube_filename) ) )
        assert os.path.isdir(this_cube_filename)
    
    def test_do_singlescale_imaging(self):
        if 'singlescale' in self.passed_steps:
            return
        target = self.key_handler.get_targets()[0]
        product = 'co21' # self.key_handler.get_line_products()[0]
        config = '7m' # self.key_handler.get_interf_configs()[0]
        self.imaging_handler.set_targets(only=[target])
        self.imaging_handler.set_line_products(only=[product])
        self.imaging_handler.set_interf_configs(only=['7m'])
        self.imaging_handler.set_no_feather_configs(True)
        self.imaging_handler.set_no_cont_products(True)
        # 
        self.imaging_handler.loop_imaging(make_directories = True, \
                                          do_dirty_image = False,
                                          do_revert_to_dirty = False,
                                          do_read_clean_mask = False, 
                                          do_multiscale_clean = False,
                                          do_revert_to_multiscale = True,
                                          do_singlescale_mask = True,
                                          do_singlescale_clean = True,
                                          do_revert_to_singlescale = False,
                                          do_export_to_fits = False, 
                                          overwrite = False)
        # 
        # check the existence of singlescale-cleaned cube image
        #assert utilsFilenames.get_cube_filename(target = target, config = config, product = product, casa = True, ext = 'singlescale') is not None
        #assert os.path.isdir(utilsFilenames.get_cube_filename(target = target, config = config, product = product, casa = True, ext = 'singlescale'))
        assert self.key_handler.get_cube_filename(target = target, config = config, product = product, casa = True, ext = 'singlescale') is not None
        this_cube_filename = self.key_handler.get_cube_filename(target = target, config = config, product = product, casa = True, ext = 'singlescale')
        logger.debug('%r %s %s'%(this_cube_filename, 'isdir?', os.path.isdir(this_cube_filename) ) )
        assert os.path.isdir(this_cube_filename)
    
    def test_do_export_to_fits(self):
        if 'exportfits' in self.passed_steps:
            return
        target = self.key_handler.get_targets()[0]
        product = 'co21' # self.key_handler.get_line_products()[0]
        config = '7m' # self.key_handler.get_interf_configs()[0]
        self.imaging_handler.set_targets(only=[target])
        self.imaging_handler.set_line_products(only=[product])
        self.imaging_handler.set_interf_configs(only=['7m'])
        self.imaging_handler.set_no_feather_configs(True)
        self.imaging_handler.set_no_cont_products(True)
        # 
        self.imaging_handler.loop_imaging(make_directories = True, \
                                          do_dirty_image = False,
                                          do_revert_to_dirty = False,
                                          do_read_clean_mask = False, 
                                          do_multiscale_clean = False,
                                          do_revert_to_multiscale = False,
                                          do_singlescale_mask = False,
                                          do_singlescale_clean = False,
                                          do_revert_to_singlescale = True,
                                          do_export_to_fits = True, 
                                          overwrite = False)
        # 
        # check the existence of exported fits files
        for ext in [None, 'model', 'mask', 'pb', 'psf']:
            #assert utilsFilenames.get_cube_filename(target = target, config = config, product = product, casa = False, ext = ext) is not None
            #assert os.path.isdir(utilsFilenames.get_cube_filename(target = target, config = config, product = product, casa = False, ext = ext))
            assert self.key_handler.get_cube_filename(target = target, config = config, product = product, casa = False, ext = ext) is not None
            this_cube_filename = self.key_handler.get_cube_filename(target = target, config = config, product = product, casa = False, ext = ext)
            logger.debug('%r %s %s'%(this_cube_filename, 'isfile?', os.path.isfile(this_cube_filename) ) )
            assert os.path.isfile(this_cube_filename)
    


if __name__ == '__main__':
    unittest.main(exit=False)
else:
    unittest.main()





