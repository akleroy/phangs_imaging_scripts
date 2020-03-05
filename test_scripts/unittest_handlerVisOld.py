
# Test under:
#   CASA 5.4.0
# 
# How to run:
#   casapy-5.4.0 --nogui --log2term -c "exec(open('test_scripts/unittest_handlerVisOld.py','r').read())"
#   casapy-5.4.0 -c "execfile('test_scripts/unittest_handlerVisOld.py')"
# 
# TODO: 
#   waiting for new uvdata handler
# 

from __future__ import print_function
import os, sys, shutil, logging, unittest
logging.basicConfig()
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
sys.path.append('.')
from phangsPipeline import handlerTemplate
from phangsPipeline import handlerKeys
from phangsPipeline import handlerVis
from phangsPipeline import utilsFilenames
if sys.version_info.major <= 2:
    reload(handlerTemplate)
    reload(handlerKeys)
    reload(handlerVis)
    reload(utilsFilenames)



class test_uvdata_handler(unittest.TestCase):
    
    key_handler = None
    uvdata_handler = None
    current_dir = None
    passed_steps = None
    
    def setUp(self):
        self.current_dir = os.getcwd()
    
    @classmethod
    def setUpClass(cls):
        # requires Python >= 2.7
        cls.key_handler = handlerKeys.KeyHandler(master_key = 'test_keys/master_key.txt', dochecks = False)
        cls.uvdata_handler = handlerVis.VisHandler(key_handler = cls.key_handler)
        #cls.key_handler.get_vis_filename = utilsFilenames.get_vis_filename
        #cls.key_handler.get_cube_filename = utilsFilenames.get_cube_filename
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
        #if self.uvdata_handler is not None:
        #    del self.uvdata_handler
        if self.current_dir is not None:
            os.chdir(self.current_dir)
        #if os.path.isdir('test'):
        #    shutil.rmtree('test')
        pass
    
    def test_initialize_key_handler(self):
        assert self.key_handler is not None
    
    def test_initialize_uvdata_handler(self):
        assert self.uvdata_handler is not None
    
    def test_do_copy(self):
        if 'copy' in self.passed_steps:
            return
        target = self.key_handler.get_targets()[0]
        product = 'co21' # self.key_handler.get_line_products()[0]
        config = '7m' # self.key_handler.get_interf_configs()[0]
        self.uvdata_handler.set_targets(only=[target])
        self.uvdata_handler.set_line_products(only=[product])
        self.uvdata_handler.set_interf_configs(only=['7m'])
        self.uvdata_handler.set_no_feather_configs(True)
        self.uvdata_handler.set_no_cont_products(True)
        self.uvdata_handler.loop_stage_uvdata(make_directories = True, do_copy = True, extra_ext = 'for_unit_test', overwrite = False)
        assert self.key_handler.get_all_input_ms(target = target, config = config) is not None
        this_ms_filenames, this_ms_filepaths = self.key_handler.get_all_input_ms(target = target, config = config)
        for i in range(len(this_ms_filenames)):
            this_ms_filename = this_ms_filenames[i]+'.ms'
            logger.debug('%r %s %s'%(this_ms_filename, 'isdir?', os.path.isdir(this_ms_filename) ) )
            assert os.path.isdir(this_ms_filename)
    
    def test_do_extract_line(self):
        if 'extract_line' in self.passed_steps:
            return
        target = self.key_handler.get_targets()[0]
        product = 'co21' # self.key_handler.get_line_products()[0]
        config = '7m' # self.key_handler.get_interf_configs()[0]
        self.uvdata_handler.set_targets(only=[target])
        self.uvdata_handler.set_line_products(only=[product])
        self.uvdata_handler.set_interf_configs(only=['7m'])
        self.uvdata_handler.set_no_feather_configs(True)
        self.uvdata_handler.set_no_cont_products(True)
        self.uvdata_handler.loop_stage_uvdata(make_directories = True, do_copy = False, do_extract_line = True, do_concat_line = True, extra_ext = 'for_unit_test', overwrite = False)
        assert self.key_handler.get_vis_filename(target = target, config = config, product = product, ext = 'for_unit_test') is not None
        this_dir = self.key_handler.get_imaging_dir_for_target(target = target)
        this_ms_filename = self.key_handler.get_vis_filename(target = target, config = config, product = product, ext = 'for_unit_test')
        logger.debug('%r %s %s'%(this_dir+os.sep+this_ms_filename, 'isdir?', os.path.isdir(this_dir+os.sep+this_ms_filename) ) )
        assert os.path.isdir(this_dir+os.sep+this_ms_filename)
    
    def test_do_extract_cont(self):
        if 'extract_cont' in self.passed_steps:
            return
        target = self.key_handler.get_targets()[0]
        product = 'cont' # self.key_handler.get_line_products()[0]
        config = '7m' # self.key_handler.get_interf_configs()[0]
        self.uvdata_handler.set_targets(only=[target])
        self.uvdata_handler.set_interf_configs(only=['7m'])
        self.uvdata_handler.set_no_feather_configs(True)
        self.uvdata_handler.set_no_line_products(True)
        self.uvdata_handler.set_no_cont_products(False)
        self.uvdata_handler.loop_stage_uvdata(make_directories = True, do_copy = False, do_extract_line = False, do_concat_line = False, do_extract_cont = True, do_concat_cont = True, extra_ext = 'for_unit_test', overwrite = False)
        assert self.key_handler.get_vis_filename(target = target, config = config, product = product, ext = 'for_unit_test') is not None
        this_dir = self.key_handler.get_imaging_dir_for_target(target = target)
        this_ms_filename = self.key_handler.get_vis_filename(target = target, config = config, product = product, ext = 'for_unit_test')
        logger.debug('%r %s %s'%(this_dir+os.sep+this_ms_filename, 'isdir?', os.path.isdir(this_dir+os.sep+this_ms_filename) ) )
        assert os.path.isdir(this_dir+os.sep+this_ms_filename)
    


if __name__ == '__main__':
    unittest.main(exit=False)
else:
    unittest.main()





