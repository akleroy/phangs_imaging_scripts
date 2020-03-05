
# Test under:
#   CASA 5.4.0
# 
# How to run:
#   casapy-5.4.0 --nogui --log2term -c "exec(open('test_scripts/unittest_utilsResolutions.py','r').read())"
#   python2.7 test_scripts/unittest_utilsResolutions.py
# 
# TODO: 
# 

from __future__ import print_function
import os, sys, shutil, logging, unittest
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
sys.path.append('.')
from phangsPipeline import utilsResolutions
if sys.version_info.major <= 2:
    reload(utilsResolutions)



class test_utils_resolutions(unittest.TestCase):
    
    current_dir = None
    
    def setUp(self):
        self.current_dir = os.getcwd()
    
    @classmethod
    def setUpClass(cls):
        # requires Python >= 2.7
        pass
        
    def tearDown(self):
        if self.current_dir is not None:
            os.chdir(self.current_dir)
        pass
    
    def test_1(self):
        assert utilsResolutions.get_tag_for_angular_resolution(3.0) == '3p00'
        assert utilsResolutions.get_tag_for_angular_resolution(25.0) == '25p00'
        assert utilsResolutions.is_angular_resolution(25.0) == True
        assert utilsResolutions.is_angular_resolution('25p0') == True
        assert utilsResolutions.is_angular_resolution('25.0arcsec') == True
        assert utilsResolutions.is_physical_resolution(25.0) == True
        assert utilsResolutions.is_physical_resolution('25.0') == False # physical resolution must have a unit
        assert utilsResolutions.is_physical_resolution('25.0pc') == True
        assert utilsResolutions.is_physical_resolution('25.0kpc') == True
        assert utilsResolutions.is_physical_resolution('25.0 Mpc') == True
    
    def test_2(self):
        #print(str(utilsResolutions.get_angular_resolution_for_res('25.0 Mpc')))
        assert utilsResolutions.get_angular_resolution_for_res('250 pc', distance=30) is not None
        print('250 pc is '+str(utilsResolutions.get_angular_resolution_for_res('250 pc', distance=30))+' arcsec')
    
    


if __name__ == '__main__':
    unittest.main(exit=False)
else:
    unittest.main()





