"""
How to run this test inside CASA:

```
sys.path.append('../casa_analysis_scripts')
sys.path.append('../analysis_scripts')
sys.path.append('.')
import importlib
#importlib.reload = reload
import phangsPipeline
importlib.reload(phangsPipeline)
importlib.reload(phangsPipeline.handlerKeys)
importlib.reload(phangsPipeline.handlerVis)
import phangsPipelineTests
importlib.reload(phangsPipelineTests)
importlib.reload(phangsPipelineTests.test_handlerVis)
phangsPipelineTests.TestingHandlerVisInCasa().run()
```

What will be tested:

```
prep_sd_for_feather

```
"""

import os, sys, shutil
import unittest


class TestingHandlerVis(unittest.TestCase):
    """docstring for TestingHandlerVis"""
    
    def __init__(self, *args, **kwargs):
        super(TestingHandlerVis, self).__init__(*args, **kwargs)
        import phangsPipeline
        self.current_dir = os.getcwd()
        self.module_dir = os.path.dirname(os.path.abspath(phangsPipeline.__path__[0]))
        self.working_dir = os.path.join(self.module_dir, 'phangsPipelineTests')
        self.test_data_dir = os.path.join(self.working_dir, 'test_data')
        self.test_keys_dir = os.path.join(self.working_dir, 'test_keys')
    
    def setup_necessary_directories(self):
        for this_dir in ['uvdata', 'singledish']: 
            assert os.path.isdir(self.test_data_dir+os.sep+this_dir), 'Testing data not found: '+self.test_data_dir+os.sep+this_dir
        for this_dir in ['cleanmasks', 'reduction']: 
            if not os.path.isdir(self.test_data_dir+os.sep+this_dir):
                os.makedirs(self.test_data_dir+os.sep+this_dir)
    
    def test_loop_stage_uvdata(self):
        import phangsPipeline
        from phangsPipeline import handlerKeys as kh
        from phangsPipeline import handlerVis as uvh
        os.chdir(self.working_dir)
        self.setup_necessary_directories()
        this_kh = kh.KeyHandler(master_key=self.test_keys_dir+os.sep+'master_key.txt')
        this_uvh = uvh.VisHandler(key_handler=this_kh)
        this_uvh.loop_stage_uvdata(
            do_all=True,
            )
        os.chdir(self.current_dir)
        assert os.path.isdir(os.path.join(self.test_data_dir, 'reduction', 'imaging', 'ngc4321', 'ngc4321_7m_co10.ms'))
        
    def tearDown(self):
        os.chdir(self.current_dir)
        #for this_dir in ['cleanmasks', 'reduction']: 
        #    if os.path.isdir(self.test_data_dir+os.sep+this_dir):
        #        shutil.rmtree(self.test_data_dir+os.sep+this_dir)



class TestingHandlerVisInCasa():
    """docstring for TestingHandlerVisInCasa"""
    
    def __init__(self):
        pass
    
    def suite(self=None):
        testsuite = unittest.TestSuite()
        testsuite.addTest(unittest.makeSuite(TestingHandlerVis))
        return testsuite
    
    def run(self):
        unittest.main(defaultTest='phangsPipelineTests.TestingHandlerVisInCasa.suite', exit=False)



if __name__ == '__main__':
    unittest.main()


