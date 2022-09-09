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
importlib.reload(phangsPipeline.casaCubeRoutines)
importlib.reload(phangsPipeline.casaFeatherRoutines)
import phangsPipelineTests
importlib.reload(phangsPipelineTests)
importlib.reload(phangsPipelineTests.test_casaFeatherRoutines)
phangsPipelineTests.TestingCasaFeatherRoutinesInCasa().run()
```

What will be tested:

```
prep_sd_for_feather

```
"""

import os, sys, shutil
import unittest


class TestingCasaFeatherRoutines(unittest.TestCase):
    """docstring for TestingCasaFeatherRoutines"""
    
    def __init__(self, *args, **kwargs):
        super(TestingCasaFeatherRoutines, self).__init__(*args, **kwargs)
        import phangsPipeline
        self.current_dir = os.getcwd()
        self.module_dir = os.path.dirname(os.path.abspath(phangsPipeline.__path__[0]))
        self.working_dir = os.path.join(self.module_dir, 'phangsPipelineTests')
        self.test_data_dir = os.path.join(self.working_dir, 'test_data')
        self.test_keys_dir = os.path.join(self.working_dir, 'test_keys')
        self.test_sd_file = os.path.join(self.test_data_dir, 'feathering', 'ALMA_TP.NGC1808.CO21_5kmsres.image.fits')
        self.test_interf_file = os.path.join(self.test_data_dir, 'feathering', 'ngc1808_12m+7m_co21_pbcorr_trimmed_k.fits')
    
    def import_test_image(self, filename, outname, overwrite = True):
        from phangsPipeline import casaStuff
        if os.path.exists(outname):
            if overwrite:
                shutil.rmtree(outname)
            else:
                return
        casaStuff.importfits(filename, outname)
        self.assertTrue(os.path.isdir(outname))
    
    # def fetch_test_data(self):
    #     import ssl
    #     import urllib.request
    #     import tarfile
    #     import numpy as np
    #     try:
    #         from astropy.io import fits
    #     except:
    #         import pyfits as fits
        
    #     self.test_sd_file = self.test_data_dir+os.sep+'ALMA_TP.NGC1808.CO21_5kmsres.image.fits'
    #     self.test_interf_file = self.test_data_dir+os.sep+'ngc1808_12m+7m_co21_pbcorr_trimmed_k.fits'
        
    #     # use_M100 = False
    #     # if use_M100:
    #     #     self.test_sd_file = self.test_data_dir+os.sep+'M100_Band3_ACA_ReferenceImages/M100_TP_CO_cube.bl.image.fits'
    #     #     self.test_interf_file = self.test_data_dir+os.sep+'M100_Band3_ACA_ReferenceImages/M100_7m_CO.image.fits'
    #     # 
    #     #     if not os.path.isfile(self.test_sd_file) and \
    #     #        not os.path.isfile(self.test_interf_file):
    #     # 
    #     #         if not os.path.isdir(self.test_data_dir):
    #     #             os.makedirs(self.test_data_dir)
    #     #         if not os.path.isfile(self.test_data_dir+os.sep+'M100_Band3_ACA_ReferenceImages.tgz'):
    #     #             with urllib.request.urlopen(
    #     #                 'https://almascience.eso.org/almadata/sciver/M100Band3ACA/M100_Band3_ACA_ReferenceImages.tgz',
    #     #                 context = ssl._create_unverified_context()
    #     #                 ) as fp_from, open(
    #     #                 self.test_data_dir+os.sep+'M100_Band3_ACA_ReferenceImages.tgz', 'wb'
    #     #                 ) as fp_to:
    #     #                 fp_to.write(fp_from.read())
    #     # 
    #     #         self.assertTrue(os.path.isfile(self.test_data_dir+os.sep+'M100_Band3_ACA_ReferenceImages.tgz'))
    #     # 
    #     #         with tarfile.open(self.test_data_dir+os.sep+'M100_Band3_ACA_ReferenceImages.tgz', 'r:gz') as tgz:
    #     #             #tgz.extract(
    #     #             #    'M100_Band3_ACA_ReferenceImages/M100_TP_CO_cube.bl.image.fits', 
    #     #             #    self.test_data_dir
    #     #             #    )
    #     #             #tgz.extract(
    #     #             #    'M100_Band3_ACA_ReferenceImages/M100_7m_CO.image.fits', 
    #     #             #    self.test_data_dir
    #     #             #    )
    #     #             tgz.extractall(
    #     #                 self.test_data_dir, 
    #     #                 members=[
    #     #                     tgz.getmember('M100_Band3_ACA_ReferenceImages/M100_TP_CO_cube.bl.image.fits'), 
    #     #                     tgz.getmember('M100_Band3_ACA_ReferenceImages/M100_7m_CO.image.fits')
    #     #                     ],
    #     #                 )
        
    #     if not os.path.isfile(self.test_sd_file):
    #         print('Error! Test data file is not found: '+os.path.isfile(self.test_sd_file))
    #     if not os.path.isfile(self.test_interf_file):
    #         print('Error! Test data file is not found: '+os.path.isfile(self.test_interf_file))
    #     self.assertTrue(os.path.isfile(self.test_sd_file))
    #     self.assertTrue(os.path.isfile(self.test_interf_file))
    
    def test_prep_sd_for_feather(self):
        from phangsPipeline import casaFeatherRoutines as cfr
        assert os.path.exists(self.test_sd_file), 'Testing data not found: '+self.test_sd_file
        assert os.path.exists(self.test_interf_file), 'Testing data not found: '+self.test_interf_file
        if self.test_interf_file.endswith('.fits'):
            self.import_test_image(
                self.test_interf_file, 
                self.test_data_dir+os.sep+'interf_data_to_feather.image', 
                )
        else:
            shutil.copy2(
                self.test_interf_file, 
                self.test_data_dir+os.sep+'interf_data_to_feather.image', 
                )
        from casatasks import imhead
        if imhead(self.test_data_dir+os.sep+'interf_data_to_feather.image')['ndim'] == 3:
            do_dropdeg = True
        else:
            do_dropdeg = False
        cfr.prep_sd_for_feather(
                sdfile_in=self.test_sd_file,
                sdfile_out=self.test_data_dir+os.sep+'sd_data_to_feather.image',
                interf_file=self.test_data_dir+os.sep+'interf_data_to_feather.image',
                do_dropdeg=do_dropdeg,
                overwrite=True,
            )
    
    def test_feather_two_cubes(self):
        from phangsPipeline import casaFeatherRoutines as cfr
        self.test_prep_sd_for_feather()
        cfr.feather_two_cubes(
                interf_file=self.test_data_dir+os.sep+'interf_data_to_feather.image',
                sd_file=self.test_data_dir+os.sep+'sd_data_to_feather.image',
                out_file=self.test_data_dir+os.sep+'feathered.image',
                do_blank=True,
                do_apodize=False,
                overwrite=True,
            )
        
    def tearDown(self):
        if os.path.exists(self.test_data_dir+os.sep+'sd_data_to_feather.image'):
           shutil.rmtree(self.test_data_dir+os.sep+'sd_data_to_feather.image')
        if os.path.exists(self.test_data_dir+os.sep+'interf_data_to_feather.image'):
           shutil.rmtree(self.test_data_dir+os.sep+'interf_data_to_feather.image')
        if os.path.exists(self.test_data_dir+os.sep+'feathered.image'):
           shutil.rmtree(self.test_data_dir+os.sep+'feathered.image')



class TestingCasaFeatherRoutinesInCasa():
    """docstring for TestingCasaFeatherRoutinesInCasa"""
    
    def __init__(self):
        pass
    
    def suite(self=None):
        testsuite = unittest.TestSuite()
        testsuite.addTest(unittest.makeSuite(TestingCasaFeatherRoutines))
        return testsuite
    
    def run(self):
        unittest.main(defaultTest='phangsPipelineTests.TestingCasaFeatherRoutinesInCasa.suite', exit=False)



if __name__ == '__main__':
    unittest.main()


