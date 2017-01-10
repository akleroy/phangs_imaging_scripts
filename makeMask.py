tested_versions = ['4.6.0','4.7.0']
this_version = (casa['build']['version']).split('-')[0]
if this_version not in tested_versions:
    print "The script hasn't been verified for this version of CASA."
    print "This version of CASA is "+this_version
    print "Tested versions are "+str(tested_versions)

import glob
import numpy as np
import scipy

# Hate this ... but directories
execfile('../scripts/auExtensions.py')

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# TIMER
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

import time
total_start_time = time.time()

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# CONTROL FLOW
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

try:
    do_mask
except NameError:    
    do_mask = True

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# INPUTS
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

abort = False

print ""
print "....................................."
print "makeMask: begins."

# ......................................
# Input/output files
# ......................................

try:
    cube_root
except NameError:
    print "makeMask: Please define a root output name via cube_root=XXX."
    abort = True

# ......................................
# Recipes
# ......................................

try:
    recipe
except NameError:
    recipe = ""

if recipe == "smoothandclip":
    print "Using recipe: "+recipe
    do_convolve = True
    do_scale_beam = True
    scale_factor = 8.0
    hi_thresh = 10.0
    lo_thresh = 10.0
    do_specsmooth = False
    spec_width = 1
    do_smooth_mask = False

if recipe == "clipandsmooth":

    # This recipe tries to make a mask that includes only signal and
    # then smooths it very heavily AFTER thresholding. The result
    # should be appropriate for a very agnostic clean. The basic
    # pipeline 

    print "Using recipe: "+recipe
    do_convolve = False
    do_specsmooth = False
    do_scale_beam = False

    hi_thresh = 5.0
    lo_thresh = 2.0

    do_reject_small_regions = True
    reject_area_in_beams = 1.0

    do_spectral_dilation = True
    spectral_dilation_iters = 1

    do_smooth_mask = True
    smooth_mask_factor = 20.

if recipe == "cprops":
    pass

if recipe == "modelmask":

    print "Using recipe: "+recipe

    use_model = True
    
    do_convolve = True
    do_scale_beam = True
    scale_factor = 2.0
    
    hi_thresh = 1.0
    lo_thresh = 1.0

    do_spectral_dilation = True
    spectral_dilation_iters = 1

    do_smooth_mask = False

# ------------
# Image to use
# ------------

try:
    use_resid
except NameError:
    use_resid = False

try:
    use_model
except NameError:
    use_model = False

# -----------
# Convolution
# -----------

try:
    do_convolve
except NameError:
    do_convolve = False

try:
    do_scale_beam 
except NameError:
    do_scale_beam = True

try:
    scale_factor 
except NameError:
    scale_factor = 2.0

try:
    target_beam 
except NameError:
    target_beam = '3.0arcsec'

try:
    do_specsmooth
except NameError:
    do_specsmooth = False

try:
    spec_width
except NameError:
    spec_width = 3

# -----------
# Threshold
# -----------

try:
    hi_thresh
except NameError:
    hi_thresh = 3

try:
    lo_thresh
except NameError:
    lo_thresh = 3

# -----------
# Pruning
# -----------

try:
    do_reject_small_regions
except NameError:
    do_reject_small_regions = False

try:
    reject_area_in_beams
except NameError:
    reject_area_in_beams = 1.0

# ----------------------
# Dilation and smoothing
# ----------------------

try:
    do_spectral_dilation
except NameError:
    do_spectral_dilation = False

try:
    spectral_dilation_iters
except NameError:
    spectral_dilation_iters = 1

try:
    do_smooth_mask
except NameError:
    do_smooth_mask = False

try:
    smooth_mask_to_res
except NameError:
    smooth_mask_to_res = 20.

# ......................................
# Other options
# ......................................

try:
    working_file
except NameError:
    working_file = "working_mask_file.image"

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# PRE-PROCESSING
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

# ...........................................................................
# Either the cube or the model or the residual serves as the base.
# ...........................................................................

print "makeMask: copying the base data."

if use_resid == True:
    base_image = cube_root+'.residual'
elif use_model == True:
    base_image = cube_root+'.model'
else:
    base_image = cube_root+'.image'

test_list = glob.glob(base_image)
if len(test_list) == 0:
    print "makeMask: cannot find the target file."

# ...........................................................................
# Get the image statistics, prefer the residual cube if available.
# ...........................................................................

# ... check for the existence of a residual file
resid_files = glob.glob(cube_root+'.residual')
if len(resid_files) == 1:
    stats = imstat(cube_root+'.residual')
else:
    stats = imstat(cube_root+'.image')

# ... the RMS comes from the median absolute deviation, rescaled
rms = stats['medabsdevmed'][0]/0.6745
print "makeMask: found a noise of "+str(rms)+" Jy/beam"

# ...........................................................................    
# Get the beam.
# ...........................................................................

header = imhead(cube_root+'.image')

# ... this is complicated by the possible presence of per-plane beams
if (header.keys()).count('restoringbeam') == 1:
    # ... the simple case
    beam = str(header['restoringbeam']['major']['value'])+'arcsec'
elif (header.keys()).count('perplanebeams') == 1:
    # ... per plane beams, pick the largest beam
    ppbdict = header['perplanebeams']['beams']
    beam = 0.0
    for plane in ppbdict.keys():
        this_plane = ppbdict[plane]
        for key in this_plane.keys():
            this_major = this_plane[key]["major"]["value"]
            if this_major > beam:
                beam = this_major
    beam = str(beam)+'arcsec'
else:
    print "makeMask: could not find a beam."
    beam = None

pix_arcsec = abs(header['incr'][0]*180./np.pi*3600.)
beam_arcsec = float((beam.split('arcsec')[0]))
pix_per_beam = beam_arcsec/pix_arcsec
beam_area_pix = (beam_arcsec/pix_arcsec/2.)**2*np.pi/log(2)

print "makeMask: found a beam of "+str(beam_arcsec)+" arcseconds"
print "makeMask: found "+str(pix_per_beam)+" pixels per beam."
print "makeMask: found a beam area of "+str(beam_area_pix)+" pixels."

# ...........................................................................
# Figure out if the data are a cube or an image.
# ...........................................................................

this_shape = header['shape']
non_degenerate_axes = np.sum(this_shape > 1)
if non_degenerate_axes == 3:
    is_cube = True
else:
    is_cube = False

if is_cube == True:
    print "makeMask: Identified the data as a cube."
else:
    print "makeMask: Identified the data as an image."

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# CONVOLUTION
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

# Convolve the image.
if do_convolve == True:

    print "makeMask: convolving the image."
    
    # Before masking, convolve the image or model to lower resolution,
    # useful to recover faint, extended emission.

    if do_scale_beam == True:
        # Scale from the current beam.
        target_beam = str(beam_arcsec*scale_factor)+'arcsec'
        os.system('rm -rf '+working_file)
        imsmooth(imagename=base_image,
                 targetres=True, major=target_beam, minor=target_beam, pa='0deg',
                 outfile=working_file, overwrite=True)
    else:
        # Use a fixed target resolution.
        os.system('rm -rf '+working_file)
        imsmooth(imagename=base_image,
                 targetres=True, major=target_beam, minor=target_beam, pa='0deg',
                 outfile=working_file, overwrite=True)

else:

    # Don't carry out a convolution, but shift the name of the file so
    # that we reach the same point.

    print "makeMask: no convolution."

    os.system('rm -rf '+working_file)
    os.system('cp -r '+base_image+' '+working_file)

# Smooth the cube spectrally
if do_specsmooth == True:

    if is_cube == True:
        print "makeMask: spectrally smoothing the image."
        os.system('rm -rf '+working_file+'.temp')
        os.system('cp -r '+working_file+' '+working_file+'.temp')

        os.system('rm -rf '+working_file)
        specsmooth(imagename=working_file+'.temp',
                   outfile=working_file,
                   function='boxcar',
                   width=spec_width,
                   dmethod="")
    else:
        print "makeMask: can't spectrally smooth an image, only a cube."
        pass

else:
    
    print "makeMask: no spectral smoothing."

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# GET STATISTICS OF THE WORKING IMAGE
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

print "makeMask: found a noise of "+str(rms)+" Jy/beam in the working image"

working_stats = imstat(working_file)
working_rms = stats['medabsdevmed'][0]/0.6745

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# THRESHOLD THE WORKING DATA
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

print "makeMask: making a high threshold image."

if use_model == False:
    rms_for_mask = working_rms
else:
    rms_for_mask = rms

os.system('rm -rf '+working_file+'.hi_mask')
immath(imagename = working_file,
       outfile = working_file+'.hi_mask',
       expr = 'iif(IM0 > '+str(hi_thresh*rms_for_mask) +',1.0,0.0)')

if lo_thresh < hi_thresh:

    print "makeMask: also making a low threshold image."

    os.system('rm -rf '+working_file+'.lo_mask')
    immath(imagename = working_file,
           outfile = working_file+'.lo_mask',
           expr = 'iif(IM0 > '+str(lo_thresh*rms_for_mask) +',1.0,0.0)')

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# READ IN THE THRESHOLDED MASK
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

print "makeMask: reading the data in to an array."

ia.open(working_file+'.hi_mask')
mask = ia.getchunk()
ia.done()

if lo_thresh < hi_thresh:
    ia.open(working_file+'.lo_mask')
    lo_mask = ia.getchunk()
    ia.done()

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# REQUIRE SUCCESSIVE CHANNELS
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

print "makeMask: requiring multiple adjacent channels."

#    new_mask = mask*(np.roll(mask,1,3) + np.roll(mask,-1,3)) >= 1

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# SMALL REGION REJECTION
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

if do_reject_small_regions:

    print "makeMask: removing small regions from the mask."

    regions, n_regions = scipy.ndimage.label(mask)                     
    myhistogram = scipy.ndimage.measurements.histogram(regions,0,n_regions+1,n_regions+1)
    object_slices = scipy.ndimage.find_objects(regions)
    for i in range(n_regions):
        if myhistogram[i+1] < reject_area_in_beams*beam_area_pix:
            mask[object_slices[i]] = 0

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# DILATION
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

if lo_thresh < hi_thresh:
    print "makeMask: dilating the high S/N mask into the low S/N mask."
    new_mask = scipy.ndimage.binary_dilation(mask,iterations=-1,mask=lo_mask)    
    mask = new_mask

#print "makeMask: dilating the mask spatially."
#spatial_dilation = 

if do_spectral_dilation:
    if is_cube == True:
        print "makeMask: dilating the mask spectrally."
        if header['axisnames'][2] == 'Frequency':
            spec_axis = 2
        else:
            spec_axis = 3
        for ii in range(spectral_dilation_iters):
            new_mask = \
                (mask + np.roll(mask,1,axis=spec_axis) + \
                     np.roll(mask,-1,axis=spec_axis)) >= 1
            mask = new_mask

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# PUT THE DATA INTO THE MASK FOR THE IMAGE
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

print "makeMask: removing the old mask."
os.system('rm -rf '+cube_root+'.mask')

print "makeMask: copying the image."
os.system('cp -r '+cube_root+'.image '+cube_root+'.mask')

print "makeMask: putting the data into a copy of the image."
ia.open(cube_root+'.mask')
ia.putchunk(mask*1.0)
ia.done()

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# SMOOTH THE MASK
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

if do_smooth_mask:
    print "makeMask: smoothing the mask then clipping."
    print "makeMask: (this operation assumes Jy/beam math)."

    os.system('rm -rf '+cube_root+'.mask.temp')
    smooth_to_beam = str(beam_arcsec*smooth_mask_factor)+'arcsec'
    imsmooth(imagename=cube_root+'.mask',
             targetres=True,
             major=smooth_to_beam, minor=smooth_to_beam, pa='0deg',
             outfile=cube_root+'.mask.temp',
             overwrite=True)

    cutoff = 1.0
    os.system('rm -rf '+cube_root+'.mask.temp.temp')
    immath(imagename = cube_root+'.mask.temp',
           outfile = cube_root+'.mask.temp.temp',
           expr = 'iif(IM0 > '+str(cutoff) +',1.0,0.0)')
    
    print "makeMask: copying the data from the smoothed mask to the original."

    ia.open(cube_root+'.mask.temp.temp')
    mask = ia.getchunk()
    ia.done()

    ia.open(cube_root+'.mask')
    ia.putchunk(mask*1.0)
    ia.done()

    os.system('rm -rf '+cube_root+'.mask.temp')
    os.system('rm -rf '+cube_root+'.mask.temp.temp')

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# CLEAN UP
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

print "makeMask: removing temporary files."

os.system('rm -rf '+working_file)
os.system('rm -rf '+working_file+'.temp')
os.system('rm -rf '+working_file+'.lo_mask')
os.system('rm -rf '+working_file+'.hi_mask')

print "makeMask: ends."
print "....................................."
