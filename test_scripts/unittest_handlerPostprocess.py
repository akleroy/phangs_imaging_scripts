
# Test under:
#   CASA 5.4.0
# 
# How to run:
#   casapy-5.4.0 -c "exec(open('test_scripts/unittest_handlerPostprocess.py','r').read())"
# 
# TODO: 
# 
# DONE: 
#   Needs to add 'do_revert_to_singlescale' argument to the 'handlerPostprocess.loop_postprocess' function.
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
from phangsPipeline import handlerPostprocess
from phangsPipeline import utilsFilenames
if sys.version_info.major <= 2:
    reload(handlerTemplate)
    reload(handlerKeys)
    reload(handlerPostprocess)
    reload(utilsFilenames)



class test_postprocess_handler(unittest.TestCase):
    
    key_handler = None
    postprocess_handler = None
    current_dir = None
    passed_steps = None
    
    def setUp(self):
        self.current_dir = os.getcwd()
    
    @classmethod
    def setUpClass(cls):
        # requires Python >= 2.7
        cls.key_handler = handlerKeys.KeyHandler(master_key = 'test_keys/master_key.txt', dochecks = False)
        cls.postprocess_handler = handlerPostprocess.PostProcessHandler(key_handler = cls.key_handler)
        #cls.key_handler.get_vis_filename = utilsFilenames.get_vis_filename
        #cls.key_handler.get_cube_filename = utilsFilenames.get_cube_filename
        #self.key_handler._key_dir = os.getcwd()+os.sep+'key_templates'+os.sep
        #self.key_handler._imaging_root = os.getcwd()+os.sep+'test'+os.sep+'imaging'+os.sep
        #self.key_handler._postprocess_root = os.getcwd()+os.sep+'test'+os.sep+'postprocess'+os.sep
        #self.key_handler._product_root = os.getcwd()+os.sep+'test'+os.sep+'product'+os.sep
        #if not os.path.isdir(self.key_handler._postprocess_root):
        #    os.makedirs(self.key_handler._postprocess_root)
        cls.passed_steps = []
        #cls.passed_steps = ['dirty','cleanmask','multiscale','singlescale']
        
    def tearDown(self):
        #if self.key_handler is not None:
        #    del self.key_handler
        #if self.postprocess_handler is not None:
        #    del self.postprocess_handler
        if self.current_dir is not None:
            os.chdir(self.current_dir)
        #if os.path.isdir('test'):
        #    shutil.rmtree('test')
        pass
    
    def test_initialize_key_handler(self):
        assert self.key_handler is not None
    
    def test_initialize_postprocess_handler(self):
        assert self.postprocess_handler is not None
    
    def test_do_prep(self):
        if 'prep' in self.passed_steps:
            return
        #target = self.key_handler.get_targets()[0]
        product = 'co21' # self.key_handler.get_line_products()[0]
        #self.postprocess_handler.set_targets(only=[target])
        self.postprocess_handler.set_line_products(only=[product])
        self.postprocess_handler.set_no_cont_products(True)
        self.postprocess_handler.set_interf_configs(only=['7m','12m+7m'])
        #self.postprocess_handler.set_no_feather_configs(True)
        # 
        self.postprocess_handler.loop_postprocess(\
            do_prep=True,
            do_feather=False, 
            do_mosaic=False,
            do_cleanup=False,
            do_convolve=False,
            feather_apod=False, 
            feather_noapod=False,
            feather_before_mosaic=False,
            feather_after_mosaic=False,
            )
        # 
        for config in self.key_handler.get_interf_configs():
            assert utilsFilenames.get_cube_filename(target = target, config = config, product = product, ext = 'trimmed', casa = True, casaext = '.image') is not None
            this_cube_filename = utilsFilenames.get_cube_filename(target = target, config = config, product = product, ext = 'trimmed', casa = True, casaext = '.image')
            logger.debug('%r %s %s'%(this_cube_filename, 'isdir?', os.path.isdir(this_cube_filename) ) )
            assert os.path.isdir(this_cube_filename)
    


if __name__ == '__main__':
    unittest.main(exit=False)
else:
    unittest.main()





