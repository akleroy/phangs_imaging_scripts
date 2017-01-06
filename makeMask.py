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
print "makeMask: Begins."

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
    do_convolve = True
    do_scale_beam = True
    scale_factor = 3.0
    hi_thresh = 5.0
    do_specsmooth = True
    spec_width = 3

if recipe == "cprops":
    pass

if recipe == "modelmask":
    use_model = True

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
else:
    print "makeMask: could not find a beam."
    beam = None

print "makeMask: found a beam of "+str(beam)+" arcseconds"

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
        target_beam = str(beam*scale_factor)+'arcsec'
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

print "makeMask: thresholding the working image."

print "makeMask: making a high threshold image."

rms_for_mask = working_rms
os.system('rm -rf '+working_file+'.hi_mask')

rms_for_mask = working_rms
immath(imagename = working_file,
       outfile = working_file+'.hi_mask',
       expr = 'iif(IM0 > '+str(hi_thresh*rms_for_mask) +',1.0,0.0)')

if lo_thresh < hi_thresh:

    print "makeMask: making a low threshold image."

    os.system('rm -rf '+working_file+'.lo_mask')

    rms_for_mask = working_rms
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

print "makeMask: removing small regions from the mask."

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# DILATION
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

if lo_thresh < hi_thresh:
    print "makeMask: dilating the high S/N mask into the low S/N mask."
    new_mask = scipy.ndimage.binary_dilation(mask,iterations=-1,mask=lo_mask)    
    mask = new_mask

print "makeMask: dilating the mask spatially."

print "makeMask: dilating the mask spectrally."

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# PUT THE DATA INTO THE MASK FOR THE IMAGE
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

print "makeMask: putting the data into existing mask."

mask_file_list = glob.glob(cube_root+'.mask')

if len(mask_file_list) == 0:
    print "makeMask: no mask found. Making a new mask file."
    os.system('cp -r '+cube_root+'.image '+cube_root+'.mask')

ia.open(cube_root+'.mask')
ia.putchunk(mask*1.0)
ia.done()

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# CLEAN UP
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

print "makeMask: removing temporary files."

os.system('rm -rf '+working_file)
os.system('rm -rf '+working_file+'.temp')
os.system('rm -rf '+working_file+'.lo_mask')
os.system('rm -rf '+working_file+'.hi_mask')

print "makeMask: Ends."
print "....................................."
