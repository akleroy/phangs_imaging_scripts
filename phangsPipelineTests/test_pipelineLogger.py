"""
How to run this test inside CASA:

```
sys.path.append('../scripts/analysis_scripts')
sys.path.append('../scripts/phangs_imaging_scripts')
import phangsPipelineTests
importlib.reload(phangsPipelineTests)
phangsPipelineTests.TestingPipelineLoggerInCasa().run()
```
"""

import os, sys
import unittest


class TestingPipelineLogger(unittest.TestCase):
    """docstring for TestingPipelineLogger"""
    
    def set_sys_path(self):
        if '__file__' in globals():
            script_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            if script_dir not in sys.path:
                sys.path.insert(1, script_dir)
        #print('sys.path', sys.path)
    
    def test_logfile(self):
        self.set_sys_path()
        from phangsPipeline.pipelineLogger import PipelineLogger
        with PipelineLogger('TestingPipelineLogger', level='DEBUG', logfile='logfile.txt') as logger:
            logger.info('testing info to logfile')
            logger.warning('testing warning to logfile')
            logger.debug('testing debug to logfile')
            logger.error('testing error to logfile')
    
    def test_all(self):
        self.set_sys_path()
        from phangsPipeline.pipelineLogger import PipelineLogger
        with PipelineLogger('TestingPipelineLogger', level='DEBUG') as logger:
            logger.info('testing info')
            logger.warning('testing warning')
            logger.debug('testing debug')
            logger.error('testing error')


class TestingPipelineLoggerInCasa():
    """docstring for TestingPipelineLoggerInCasa"""
    
    def __init__(self):
        pass
    
    def suite(self=None):
        #modules = (
        #    'TestingPipelineLogger',
        #)
        testsuite = unittest.TestSuite()
        #for module in map(__import__, modules):
        #    testsuite.addTest(unittest.findTestCases(module))
        testsuite.addTest(unittest.makeSuite(TestingPipelineLogger))
        return testsuite
    
    def run(self):
        #del sys.modules['phangsPipelineTests']
        #del sys.modules['phangsPipelineTests.test_pipelineLogger']
        #import phangsPipelineTests; phangsPipelineTests.TestingPipelineLoggerInCasa().run() 
        unittest.main(defaultTest='phangsPipelineTests.TestingPipelineLoggerInCasa.suite', exit=False)


if __name__ == '__main__':
    unittest.main()


