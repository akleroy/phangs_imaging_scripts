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
importlib.reload(phangsPipelineTests.test_casaCubeRoutines)
phangsPipelineTests.TestingCasaCubeRoutinesInCasa().run()
```
"""

import os, sys, shutil
import unittest


class TestingCasaCubeRoutines(unittest.TestCase):
    """docstring for TestingCasaCubeRoutines"""
    
    def __init__(self, *args, **kwargs):
        super(TestingCasaCubeRoutines, self).__init__(*args, **kwargs)
        import phangsPipeline
        self.current_dir = os.getcwd()
        self.module_dir = os.path.dirname(os.path.abspath(phangsPipeline.__path__[0]))
        self.working_dir = os.path.join(self.module_dir, 'phangsPipelineTests')
        self.test_data_dir = os.path.join(self.working_dir, 'test_data')
        self.test_keys_dir = os.path.join(self.working_dir, 'test_keys')
        for filename in ['test.fits', 'test.image', 'test2.image', 'test3.image']:
            filepath = os.path.join(self.test_data_dir, filename)
            if os.path.isfile(filepath):
                os.remove(filepath)
            elif os.path.isdir(filepath):
                shutil.rmtree(filepath)
    
    def make_test_image(self, filename, nx = 3000, ny = 3000, nchan = 393, nstokes = None, radius = 800.0, overwrite = True):
        import numpy as np
        try:
            from astropy.io import fits
        except:
            import pyfits as fits
        
        if os.path.isfile(filename): 
            if overwrite:
                shutil.move(filename, filename+'.backup')
            else:
                return
        
        data = np.zeros([nchan, ny, nx], dtype=np.float32)
        header = fits.Header()
        header['NAXIS'] = 3
        header['NAXIS1'] = nx
        header['NAXIS2'] = ny
        header['NAXIS3'] = nchan
        header['RADESYS'] = 'FK5'
        header['SPECSYS'] = 'LSRK'
        header['EQUINOX'] = 2000.0
        header['CTYPE1'] = 'RA---TAN'
        header['CTYPE2'] = 'DEC--TAN'
        header['CTYPE3'] = 'VRAD'
        header['CUNIT1'] = 'deg'
        header['CUNIT2'] = 'deg'
        header['CUNIT3'] = 'm/s'
        header['CRVAL1'] = 1.0
        header['CRVAL2'] = 1.0
        header['CRVAL3'] = 1000.0
        header['CRPIX1'] = 1.0
        header['CRPIX2'] = 1.0
        header['CRPIX3'] = 1.0
        header['CDELT1'] = -0.2/3600.0
        header['CDELT2'] = 0.2/3600.0
        header['CDELT3'] = 1000.0
        header['RESTFRQ'] = 2.30538e11
        header['TELESCOP'] = 'ALMA'
        header['OBJECT'] = ' '
        header['BUNIT'] = 'Jy/beam'
        header['BTYPE'] = 'Intensity'
        header['BMAJ'] = 0.3/3600.0
        header['BMIN'] = 0.3/3600.0
        header['BPA'] = 0.0
        if nstokes is not None:
            header['NAXIS'] = 4
            header['NAXIS4'] = nstokes
            header['CTYPE4'] = 'STOKES'
            header['CUNIT4'] = ''
            header['CRVAL4'] = 1.0
            header['CRPIX4'] = 1.0
            header['CDELT4'] = 1.0
        gy, gx = np.mgrid[:ny, :nx]
        cy, cx = (ny-1)/2.0, (nx-1)/2.0
        gr2 = ((gy-cy)**2+(gx-cx)**2)
        radius2 = radius**2
        for ichan in range(nchan):
            #sig = 30.0 - np.abs(ichan-(nchan-1)/2.0)/nchan*15.0
            #data[ichan] = np.exp(-0.5*gr2/(sig**2))
            #data[ichan][gr2>radius**2] = np.nan
            data[ichan][gr2<=radius2] = 1.0
        hdu = fits.PrimaryHDU(data=data, header=header)
        hdu.writeto(filename)
        
        self.assertTrue(os.path.isfile(filename))
        
        if os.path.isfile(filename) and os.path.isfile(filename+'.backup'):
            os.remove(filename+'.backup')
    
    def import_test_image(self, filename, outname, overwrite = True):
        from phangsPipeline import casaStuff
        if os.path.exists(outname):
            if overwrite:
                shutil.rmtree(outname)
            else:
                return
        casaStuff.importfits(filename, outname)
        self.assertTrue(os.path.isdir(outname))
    
    def test_get_mask(self):
        from phangsPipeline import casaCubeRoutines as ccr
        self.make_test_image(self.test_data_dir+os.sep+'test.fits', overwrite=False)
        self.import_test_image(self.test_data_dir+os.sep+'test.fits', 
                               self.test_data_dir+os.sep+'test.image', overwrite=False)
        self.assertTrue(os.path.isdir(self.test_data_dir+os.sep+'test.image'))
        ccr.get_mask(self.test_data_dir+os.sep+'test.image')
        shutil.rmtree(self.test_data_dir+os.sep+'test.image')
    
    def test_copy_mask(self):
        from phangsPipeline import casaCubeRoutines as ccr
        self.make_test_image(self.test_data_dir+os.sep+'test.fits', overwrite=False)
        self.import_test_image(self.test_data_dir+os.sep+'test.fits', 
                               self.test_data_dir+os.sep+'test.image', overwrite=False)
        self.import_test_image(self.test_data_dir+os.sep+'test.fits', 
                               self.test_data_dir+os.sep+'test2.image', overwrite=False)
        ccr.copy_mask(self.test_data_dir+os.sep+'test.image', 
                      self.test_data_dir+os.sep+'test2.image')
        shutil.rmtree(self.test_data_dir+os.sep+'test.image')
        shutil.rmtree(self.test_data_dir+os.sep+'test2.image')
    
    def test_multiply_cube_by_value(self):
        from phangsPipeline import casaStuff
        from phangsPipeline import casaCubeRoutines as ccr
        self.make_test_image(self.test_data_dir+os.sep+'test.fits', overwrite=False)
        self.import_test_image(self.test_data_dir+os.sep+'test.fits', 
                               self.test_data_dir+os.sep+'test3.image', overwrite=False)
        ccr.multiply_cube_by_value(self.test_data_dir+os.sep+'test3.image', 
                                   5.0, 'Jy/beam')
        shutil.rmtree(self.test_data_dir+os.sep+'test3.image')
        
    def tearDown(self):
        if os.path.isfile(self.test_data_dir+os.sep+'test.fits'):
            os.remove(self.test_data_dir+os.sep+'test.fits')



class TestingCasaCubeRoutinesInCasa():
    """docstring for TestingCasaCubeRoutinesInCasa"""
    
    def __init__(self):
        pass
    
    def suite(self=None):
        testsuite = unittest.TestSuite()
        testsuite.addTest(unittest.makeSuite(TestingCasaCubeRoutines))
        return testsuite
    
    def run(self):
        unittest.main(defaultTest='phangsPipelineTests.TestingCasaCubeRoutinesInCasa.suite', exit=False)



if __name__ == '__main__':
    unittest.main()


