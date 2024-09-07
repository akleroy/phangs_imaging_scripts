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
    
    def __init__(self, *args, **kwargs):
        super(TestingPipelineLogger, self).__init__(*args, **kwargs)
        import phangsPipeline
        self.current_dir = os.getcwd()
        self.module_dir = os.path.dirname(os.path.abspath(phangsPipeline.__path__[0]))
        self.working_dir = os.path.join(self.module_dir, 'phangsPipelineTests')
        self.test_data_dir = os.path.join(self.working_dir, 'test_data')
        self.test_log_file = os.path.join(self.working_dir, 'test_log_file.txt')
    
    def test_logfile(self):
        from phangsPipeline.pipelineLogger import PipelineLogger
        with PipelineLogger('TestingPipelineLogger', level='DEBUG', logfile=self.test_log_file) as logger:
            logger.info('testing info to logfile')
            logger.warning('testing warning to logfile')
            logger.debug('testing debug to logfile')
            logger.error('testing error to logfile')
    
    def test_all(self):
        from phangsPipeline.pipelineLogger import PipelineLogger
        with PipelineLogger('TestingPipelineLogger', level='DEBUG') as logger:
            logger.info('testing info')
            logger.warning('testing warning')
            logger.debug('testing debug')
            logger.error('testing error')
        
    def tearDown(self):
        if os.path.isfile(self.test_log_file):
            os.remove(self.test_log_file)



class TestingPipelineLoggerInCasa():
    """docstring for TestingPipelineLoggerInCasa"""
    
    def __init__(self):
        pass
    
    def suite(self=None):
        testsuite = unittest.TestSuite()
        testsuite.addTest(unittest.makeSuite(TestingPipelineLogger))
        return testsuite
    
    def run(self):
        unittest.main(defaultTest='phangsPipelineTests.TestingPipelineLoggerInCasa.suite', exit=False)


if __name__ == '__main__':
    unittest.main()


