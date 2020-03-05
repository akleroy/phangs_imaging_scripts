
# Test under:
#   CASA 5.4.0
# 
# How to run:
#   casapy-5.4.0 --nogui --log2term -c "exec(open('test_scripts/unittest_handlerKeys.py','r').read())"
#   python2.7 test_scripts/unittest_handlerKeys.py
#   

from __future__ import print_function
import os, sys, shutil, logging, unittest
logging.basicConfig()
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
sys.path.append('.')
from phangsPipeline import handlerKeys
from phangsPipeline import utilsFilenames
from phangsPipeline import utilsResolutions
if sys.version_info.major <= 2:
    reload(handlerKeys)
    reload(utilsFilenames)
    reload(utilsResolutions)



class test_handlerKeys(unittest.TestCase):
    
    key_handler = None
    
    def setUp(self):
        self.current_dir = os.getcwd()
    
    @classmethod
    def setUpClass(cls):
        # requires Python >= 2.7
        cls.key_handler = handlerKeys.KeyHandler(master_key = 'test_keys/master_key.txt', dochecks = False)
        cls.key_handler._key_dir = os.getcwd()+os.sep+'key_templates'+os.sep # do not forget +os.sep
        cls.key_handler.check_key_existence()
        
    def tearDown(self):
        if self.current_dir is not None:
            os.chdir(self.current_dir)
    
    def test_initialize_key_handler(self):
        assert self.key_handler is not None
    
    def test_get_targets(self):
        assert self.key_handler.get_targets() is not None
        assert self.key_handler.get_targets() != []
    
    def test_get_configs(self):
        assert self.key_handler.get_interf_configs() is not None
        assert self.key_handler.get_interf_configs() != []
        assert self.key_handler.get_feather_configs() is not None
        assert self.key_handler.get_feather_configs() != []
    
    def test_get_products(self):
        assert self.key_handler.get_line_products() is not None
        assert self.key_handler.get_line_products() != []
        assert self.key_handler.get_continuum_products() is not None
        assert self.key_handler.get_continuum_products() != []
    
    def test_get_imaging_dir(self):
        assert self.key_handler.get_imaging_dir_for_target(self.key_handler.get_targets()[0]) is not None
        assert self.key_handler.get_imaging_dir_for_target(self.key_handler.get_targets()[0]) != ''
    
    def test_get_imaging_recipes(self):
        assert self.key_handler.get_imaging_recipes(config=self.key_handler.get_interf_configs()[0], product=self.key_handler.get_line_products()[0], stage=handlerKeys.VALID_IMAGING_STAGES[0]) is not None
        assert self.key_handler.get_imaging_recipes(config=self.key_handler.get_interf_configs()[0], product=self.key_handler.get_line_products()[0], stage=handlerKeys.VALID_IMAGING_STAGES[0]) != []
    
    def test_get_system_velocity_and_velocity_width_for_target(self):
        assert self.key_handler.get_system_velocity_and_velocity_width_for_target(self.key_handler.get_targets()[0]) is not None
        #assert self.key_handler.get_system_velocity_and_velocity_width_for_target('nono') is not None
    
    def test_get_distance(self):
        config = '7m'
        assert self.key_handler.get_res_for_config(config) is not None
        res_list = self.key_handler.get_res_for_config(config)
        for this_res in res_list:
            target = 'ngc3621_1'
            res_tag = utilsResolutions.get_tag_for_res(this_res)
            #logger.debug('target: '+target+', res_tag: '+str(res_tag))
            #logger.debug('target: '+target+', get_distance: '+str(self.key_handler.get_distance_for_target(target)))
            assert self.key_handler.get_distance_for_target(target) is not None
            distance = self.key_handler.get_distance_for_target(target)
            res_arcsec = utilsResolutions.get_angular_resolution_for_res(this_res, distance = distance)
            #logger.debug('target: '+target+', distance: '+str(distance)+' Mpc, this_res: '+this_res+', res_arcsec: '+str(res_arcsec))
            assert res_arcsec is not None
    
    def test_get_filenames(self):
        product = 'co21'
        if self.key_handler.get_interf_configs() is not None:
            #logger.debug('')
            for config in self.key_handler.get_interf_configs():
                for target in self.key_handler.get_all_targets():
                    for this_target, this_project, this_arraytag, this_obsnum in self.key_handler.loop_over_input_ms():
                        ms_filename = self.key_handler.get_file_for_input_ms(this_target, this_project, this_arraytag, this_obsnum)
                        vis_filename = '%s_%s_%s_%s.ms'%(this_target, this_project, this_arraytag, this_obsnum)
                        #logger.debug('target: %r, product: %r, config: %r, ms_filename: %r, vis_filename: %r'%(target, product, config, ms_filename, vis_filename))
                        assert ms_filename is not None
                        assert ms_filename != ''
                    vis_filename = utilsFilenames.get_vis_filename(target=target, product=product, config=config)
                    assert vis_filename is not None
                    assert vis_filename != ''
                    #logger.debug('target: %r, product: %r, config: %r, vis_filename: %r'%(target, product, config, vis_filename))



if __name__ == '__main__':
    unittest.main(exit=False)
else:
    unittest.main()



