# This script demonstrates generation of a "PHANGS product suite" for
# a data set that doesn't have a key file suite set up.

from spectral_cube import SpectralCube, BooleanArrayMask
from astropy.io import fits
from astropy.wcs import wcs
import scipy.ndimage as nd
import numpy as np
import scConvolution as scc
import scNoiseRoutines as scn
import scMaskingRoutines as scm
import scMoments as sco
from astropy import units as u

# ####################################################
# Definitions/tuning
# ####################################################

working_dir = '../working_data/'

infile = '12co10-05kms-d-30m-clean.fits'

outfile = 'ic342_noema+30m_co10_5kms.fits'

dist = "3.45Mpc"
target_angular = ['8arcsec','15arcsec','21arcsec']
target_physical = ['150pc','500pc','1000pc','1500pc']

# ####################################################
# Convolution
# ####################################################

this_infile = working_dir+outfile[ii]

for this_ext in target_angular:
    print("Convolution to "+this_ext)
    this_outfile = this_infile.replace('.fits','_'+this_ext+'.fits')
    scc.smooth_cube(
        incube=this_infile,
        outfile=this_outfile,
        angular_resolution=this_ext,
        linear_resolution=None,
        distance=None,
        velocity_resolution=None,
        nan_treatment='interpolate',
        tol=None,
        make_coverage_cube=True,
        collapse_coverage=True,
        coveragefile=this_outfile.replace('.fits','_coverage.fits'),
        coverage2dfile=this_outfile.replace('.fits','_coverage2d.fits'),
        dtype=np.float32,
        overwrite=True
    )

for this_ext in target_physical:
    print("Convolution to "+this_ext)
    this_outfile = this_infile.replace('.fits','_'+this_ext+'.fits')
    scc.smooth_cube(
        incube=this_infile,
        outfile=this_outfile,
        angular_resolution=None,
        linear_resolution=this_ext,
        distance=dist,
        velocity_resolution=None,
        nan_treatment='interpolate',
        tol=None,
        make_coverage_cube=True,
        collapse_coverage=True,
        coveragefile=this_outfile.replace('.fits','_coverage.fits'),
        coverage2dfile=this_outfile.replace('.fits','_coverage2d.fits'),
        dtype=np.float32,
        overwrite=True
    )

# ####################################################
# Noise Estimation
# ####################################################

# native
this_base_infile = working_dir+outfile[ii]
print("Noise native res for ", this_base_infile)
scn.recipe_phangs_noise(
    incube=this_base_infile,
    outfile=this_base_infile.replace('.fits','_noise.fits'),
    overwrite=True)
    
# angular resolutions
for this_ext in target_angular:
    this_infile = this_base_infile.replace('.fits','_'+this_ext+'.fits')
    print("Noise for "+this_ext, this_infile)
    scn.recipe_phangs_noise(
        incube=this_infile,
        outfile=this_infile.replace('.fits','_noise.fits'),
        overwrite=True)

# physical resolutions
for this_ext in target_physical:
    this_infile = this_base_infile.replace('.fits','_'+this_ext+'.fits')
    print("Noise for "+this_ext, this_infile)
    scn.recipe_phangs_noise(
        incube=this_infile,
        outfile=this_infile.replace('.fits','_noise.fits'),
        overwrite=True)

# ####################################################
# Strict Mask Creation
# ####################################################

mask_kwargs = {'hi_thresh':5.0, 'lo_thresh':2.0, 
               'hi_nchan':2, 'lo_nchan':2}

this_base_infile = working_dir+outfile[ii]

# native
this_infile = this_base_infile
print("Strict mask for ", this_infile)
scm.recipe_phangs_strict_mask(
    this_infile, 
    this_infile.replace('.fits','_noise.fits'),
    outfile=this_infile.replace('.fits','_strictmask.fits'),
    mask_kwargs=mask_kwargs,
    overwrite=True)

for this_ext in target_angular:
    this_infile = this_base_infile.replace('.fits','_'+this_ext+'.fits')
    print("Strict mask for ", this_infile)
    scm.recipe_phangs_strict_mask(
        this_infile,
        this_infile.replace('.fits','_noise.fits'),
        outfile=this_infile.replace('.fits','_strictmask.fits'),
        mask_kwargs=mask_kwargs,
        overwrite=True)

for this_ext in target_physical:
    this_infile = this_base_infile.replace('.fits','_'+this_ext+'.fits')
    print("Strict mask for ", this_infile)
    scm.recipe_phangs_strict_mask(
        this_infile, 
        this_infile.replace('.fits','_noise.fits'),
        outfile=this_infile.replace('.fits','_strictmask.fits'),
        mask_kwargs=mask_kwargs,
        overwrite=True)
    
# ####################################################
# Broad Mask Creation
# ####################################################

this_base_infile = working_dir+outfile[ii]
template_mask = this_base_infile.replace('.fits','_strictmask.fits')
broad_mask_file = this_base_infile.replace('.fits','_broadmask.fits')

list_of_masks = [template_mask]

for this_ext in target_angular:
    this_infile = this_base_infile.replace('.fits','_'+this_ext+'.fits')
    this_mask = this_infile.replace('.fits','_strictmask.fits')
    list_of_masks.append(this_mask)
    
for this_ext in target_physical:
    this_infile = this_base_infile.replace('.fits','_'+this_ext+'.fits')
    this_mask = this_infile.replace('.fits','_strictmask.fits')
    list_of_masks.append(this_mask)
    
scm.recipe_phangs_broad_mask(
    template_mask, 
    outfile=broad_mask_file, 
    list_of_masks = list_of_masks,
    grow_xy = 2, 
    grow_v = 2,
    return_spectral_cube=False, overwrite=True,
    recipe='anyscale', fraction_of_scales=0.25)
 
# ####################################################
# Moments
# ####################################################

mom_list = {
    'strictmom0':{
        'algorithm':'mom0',
        'mask':'_strictmask.fits',
        'ext':'_strict_mom0.fits',
        'ext_error':'_strict_emom0.fits',
        'kwargs':{},
        },
    'broadmom0':{
        'algorithm':'mom0',
        'mask':'_broadmask.fits',
        'ext':'_broad_mom0.fits',
        'ext_error':'_broad_emom0.fits',
        'kwargs':{},
        },
    'strictmom1':{
        'algorithm':'mom1',
        'mask':'_strictmask.fits',
        'ext':'_strict_mom1.fits',
        'ext_error':'_strict_emom1.fits',
        'kwargs':{},
        },
    'broadmom1':{
        'algorithm':'mom1',
        'mask':'_broadmask.fits',
        'ext':'_broad_mom1.fits',
        'ext_error':'_broad_emom1.fits',
        'kwargs':{},
        },
    'strictmom2':{
        'algorithm':'mom2',
        'mask':'_strictmask.fits',
        'ext':'_strict_mom2.fits',
        'ext_error':'_strict_emom2.fits',
        'kwargs':{},
        },
    'strictew':{
        'algorithm':'ew',
        'mask':'_strictmask.fits',
        'ext':'_strict_ew.fits',
        'ext_error':'_strict_eew.fits',
        'kwargs':{},
        },
    'tpeak':{
        'algorithm':'tpeak',
        'mask':'_broadmask.fits',
        'ext':'_broad_tpeak.fits',
        'ext_error':'_broad_etpeak.fits',
        'kwargs':{},
        },
    }

# native
this_base_infile = working_dir+outfile[ii]
this_infile = this_base_infile
print("Moments for ", this_infile)

for this_mom_key in mom_list:      
    this_mom = mom_list[this_mom_key]
    if this_mom['mask'] == '_broadmask.fits':
        this_mask_file = this_base_infile.replace('.fits',this_mom['mask'])
    else:
        this_mask_file = this_infile.replace('.fits',this_mom['mask'])
    this_noise_file = this_infile.replace('.fits','_noise.fits')
    this_out_file = this_infile.replace('.fits',this_mom['ext'])
    this_error_file = this_infile.replace('.fits',this_mom['ext_error'])
    sco.moment_generator(
        this_infile,
        mask=this_mask_file, noise=this_noise_file,
        moment=this_mom['algorithm'], momkwargs=this_mom['kwargs'],
        outfile=this_out_file, errorfile=this_error_file,
        channel_correlation=None,
    )

# angular resolutions
for this_ext in target_angular:
    this_infile = this_base_infile.replace('.fits','_'+this_ext+'.fits')
    for this_mom_key in mom_list:      
        this_mom = mom_list[this_mom_key]
        if this_mom['mask'] == '_broadmask.fits':
            this_mask_file = this_base_infile.replace('.fits',this_mom['mask'])
        else:
            this_mask_file = this_infile.replace('.fits',this_mom['mask'])
        this_noise_file = this_infile.replace('.fits','_noise.fits')
        this_out_file = this_infile.replace('.fits',this_mom['ext'])
        this_error_file = this_infile.replace('.fits',this_mom['ext_error'])
        sco.moment_generator(
            this_infile,
            mask=this_mask_file, noise=this_noise_file,
            moment=this_mom['algorithm'], momkwargs=this_mom['kwargs'],
            outfile=this_out_file, errorfile=this_error_file,
            channel_correlation=None,
        )

# physical resolutions
for this_ext in target_physical:
    this_infile = this_base_infile.replace('.fits','_'+this_ext+'.fits')
    for this_mom_key in mom_list:      
        this_mom = mom_list[this_mom_key]
        if this_mom['mask'] == '_broadmask.fits':
            this_mask_file = this_base_infile.replace('.fits',this_mom['mask'])
        else:
            this_mask_file = this_infile.replace('.fits',this_mom['mask'])
        this_noise_file = this_infile.replace('.fits','_noise.fits')
        this_out_file = this_infile.replace('.fits',this_mom['ext'])
        this_error_file = this_infile.replace('.fits',this_mom['ext_error'])
        sco.moment_generator(
            this_infile,
            mask=this_mask_file, noise=this_noise_file,
            moment=this_mom['algorithm'], momkwargs=this_mom['kwargs'],
            outfile=this_out_file, errorfile=this_error_file,
            channel_correlation=None,
        )
        
