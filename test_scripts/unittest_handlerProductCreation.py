
# Test under:
#   python2.7
# 
# How to run:
#   python2.7 test_scripts/unittest_handlerProductCreation.py
# 
# TODO: 
# 
# DONE: 
# 

from __future__ import print_function
import os, sys, shutil, logging, unittest, inspect
logging.basicConfig()
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
sys.path.append('.')
from phangsPipeline import handlerTemplate
from phangsPipeline import handlerKeys
from phangsPipeline import handlerDerived
from phangsPipeline import utilsFilenames
from phangsPipeline import utilsResolutions
if sys.version_info.major <= 2:
    reload(handlerTemplate)
    reload(handlerKeys)
    reload(handlerDerived)
    reload(utilsFilenames)
    reload(utilsResolutions)



class test_productcreation_handler(unittest.TestCase):
    
    key_handler = None
    productcreation_handler = None
    current_dir = None
    passed_steps = None
    
    def setUp(self):
        self.current_dir = os.getcwd()
    
    @classmethod
    def setUpClass(cls):
        # requires Python >= 2.7
        cls.key_handler = handlerKeys.KeyHandler(master_key = 'test_keys/master_key.txt', dochecks = False)
        cls.productcreation_handler = handlerDerived.DerivedHandler(key_handler = cls.key_handler)
        #for name, obj in inspect.getmembers(handlerDerived):
        #    if inspect.isclass(obj):
        #        cls.productcreation_handler = obj(key_handler = cls.key_handler)
        #        break
        cls.passed_steps = []
        #cls.passed_steps = ['step1','step2']
        
    def tearDown(self):
        #if self.key_handler is not None:
        #    del self.key_handler
        #if self.productcreation_handler is not None:
        #    del self.productcreation_handler
        if self.current_dir is not None:
            os.chdir(self.current_dir)
        #if os.path.isdir('test'):
        #    shutil.rmtree('test')
        pass
    
    def test_initialize_key_handler(self):
        assert self.key_handler is not None
    
    def test_initialize_productcreation_handler(self):
        assert self.productcreation_handler is not None
    
    def test_1(self):
        if 'prep' in self.passed_steps:
            return
        target = self.key_handler.get_targets()[0]
        product = 'co21' # self.key_handler.get_line_products()[0]
        #self.productcreation_handler.set_targets(only=[target])
        self.productcreation_handler.set_line_products(only=[product])
        self.productcreation_handler.set_no_cont_products(True)
        self.productcreation_handler.set_interf_configs(only=['7m'])
        self.productcreation_handler.set_no_feather_configs(True)
        # 
        for config in self.productcreation_handler.get_all_configs():
            for res in self.key_handler.get_res_for_config(config):
                res_tag = utilsResolutions.get_tag_for_res(res)
                postprocess_root = self.key_handler.get_postprocess_dir_for_target(target=target, changeto=False)
                assert utilsFilenames.get_cube_filename(target = target, config = config, product = product, ext = 'pbcorr_trimmed_k_res'+res_tag, casa = False) is not None
                this_cube_filename = utilsFilenames.get_cube_filename(target = target, config = config, product = product, ext = 'pbcorr_trimmed_k_res'+res_tag, casa = False)
                logger.debug('%r %s %s'%(postprocess_root+this_cube_filename, 'isfile?', os.path.isfile(postprocess_root+this_cube_filename) ) )
                #assert os.path.isfile(postprocess_root+this_cube_filename)
        # 
        self.productcreation_handler.loop_make_products(\
            do_signalmask_moment_maps=True,
            do_hybridmask_moment_maps=True, 
            extra_ext_out = '_for_unittest', 
            )
        # 
        for config in self.productcreation_handler.get_all_configs():
            for res in self.key_handler.get_res_for_config(config):
                res_tag = utilsResolutions.get_tag_for_res(res)
                productcreation_root = self.key_handler.get_derived_dir_for_target(target=target, changeto=False)
                for tag in ['mom0']:
                    assert utilsFilenames.get_cube_filename(target = target, config = config, product = product, ext = 'pbcorr_trimmed_k_res'+res_tag+'_for_unittest_broad_'+tag, casa = False) is not None
                    this_cube_filename = utilsFilenames.get_cube_filename(target = target, config = config, product = product, ext = 'pbcorr_trimmed_k_res'+res_tag+'_for_unittest_broad_'+tag, casa = False)
                    logger.debug('%r %s %s'%(productcreation_root+this_cube_filename, 'isfile?', os.path.isfile(productcreation_root+this_cube_filename) ) )
                    #assert os.path.isfile(productcreation_root+this_cube_filename)
    


if __name__ == '__main__':
    unittest.main(exit=False)
else:
    unittest.main()





