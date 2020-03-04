
# Test under:
#   CASA 5.4.0
# 
# How to run:
#   casapy-5.4.0 --nogui --log2term -c "exec(open('test_scripts/unittest_handlerKeys.py','r').read())"

from __future__ import print_function
import os, sys, unittest
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
        #cls.key_handler = handlerKeys.KeyHandler(master_key = 'test_keys/master_key.txt', dochecks = False)
        cls.key_handler = handlerKeys.KeyHandler(master_key = 'test_keys/master_key.txt', dochecks = False)
        
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
        target = 'ngc3621'
        config = '7m'
        assert self.key_handler.get_res_for_config(config) is not None
        res_list = self._kh.get_res_for_config(config)
        for this_res in res_list:
            res_tag = utilsResolutions.get_tag_for_res(this_res)
            assert self._kh.get_distance_for_target(target) is not None
            distance = self._kh.get_distance_for_target(target)
            res_arcsec = utilsResolutions.get_angular_resolution_for_res(this_res, distance = distance)
            print('this_res', this_res, 'res_arcsec', res_arcsec)
            assert res_arcsec is not None


if __name__ == '__main__':
    unittest.main(exit=False)
else:
    unittest.main()



