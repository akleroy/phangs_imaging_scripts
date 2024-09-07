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
importlib.reload(phangsPipeline.handlerImaging)
importlib.reload(phangsPipeline.casaImagingRoutines)
import phangsPipelineTests
importlib.reload(phangsPipelineTests)
importlib.reload(phangsPipelineTests.test_handlerImaging)
phangsPipelineTests.TestingHandlerImagingInCasa().run()
```

What will be tested:

```
prep_sd_for_feather

```
"""

import os, sys, shutil
import unittest


class TestingHandlerImaging(unittest.TestCase):
    """docstring for TestingHandlerImaging"""
    
    def __init__(self, *args, **kwargs):
        super(TestingHandlerImaging, self).__init__(*args, **kwargs)
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
    
    # def test_task_initialize_clean_call(self):
    #     import phangsPipeline
    #     from phangsPipeline import handlerKeys as kh
    #     from phangsPipeline import handlerImaging as imh
    #     os.chdir(self.working_dir)
    #     self.setup_necessary_directories()
    #     this_kh = kh.KeyHandler(master_key=self.test_keys_dir+os.sep+'master_key.txt')
    #     this_imh = imh.ImagingHandler(key_handler=this_kh)
    #     os.chdir(os.path.join(self.test_data_dir, 'reduction', 'imaging', 'ngc4321'))
    #     clean_call = this_imh.task_initialize_clean_call(
    #         target='ngc4321', config='7m', product='co10',
    #         )
    #     #print('clean_call', clean_call)
    #     cell, imsize = this_imh.task_pick_cell_and_imsize(
    #         clean_call=clean_call,
    #         check_overrides=True,
    #         target='ngc4321', config='7m', product='co10',
    #         )
    #     os.chdir(self.current_dir)
    
    def test_loop_imaging(self):
        import phangsPipeline
        from phangsPipeline import handlerKeys as kh
        from phangsPipeline import handlerImaging as imh
        os.chdir(self.working_dir)
        self.setup_necessary_directories()
        this_kh = kh.KeyHandler(master_key=self.test_keys_dir+os.sep+'master_key.txt')
        this_imh = imh.ImagingHandler(key_handler=this_kh)
        this_imh.loop_imaging(
            do_all=True,
            )
        os.chdir(self.current_dir)
    
    def tearDown(self):
        os.chdir(self.current_dir)
        #for this_dir in ['cleanmasks', 'reduction']: 
        #    if os.path.isdir(self.test_data_dir+os.sep+this_dir):
        #        shutil.rmtree(self.test_data_dir+os.sep+this_dir)



class TestingHandlerImagingInCasa():
    """docstring for TestingHandlerImagingInCasa"""
    
    def __init__(self):
        pass
    
    def suite(self=None):
        testsuite = unittest.TestSuite()
        testsuite.addTest(unittest.makeSuite(TestingHandlerImaging))
        return testsuite
    
    def run(self):
        unittest.main(defaultTest='phangsPipelineTests.TestingHandlerImagingInCasa.suite', exit=False)



if __name__ == '__main__':
    unittest.main()


