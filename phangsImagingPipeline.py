tested_versions = ['4.6.0','4.7.0','4.7.1','4.7.2', '5.0.0']
this_version = (casa['build']['version']).split('-')[0]
if this_version not in tested_versions:
    print "The script hasn't been verified for this version of CASA."
    print "This version of CASA is "+this_version
    print "Tested versions are "+str(tested_versions)

execfile('../scripts/auExtensions.py')

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# Inputs
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

print "****************************************"
print "********** phangsImagingPipeline.py ****"
print "****************************************"

abort = False

try:
    input_vis
except NameError:
    print "Specify an input visibility. Aborting."
    abort = True

try:
    cube_root
except NameError:
    print "Specify a cube root. Aborting."
    abort = True

try:
    array
except NameError:
    print "Specify the array being used. Aborting."
    array = None
    abort = True

try:
    clean_mask
except NameError:
    print "No clean mask. Will use a primary beam based mask."
    clean_mask = "None"

try:
    pb_limit
except NameError:
    print "Primary beam coverage defaulting to 0.25."
    pb_limit = 0.25

try:
    smallscalebias
except NameError:
    print "Small scale bias defaulting to 0.6."
    smallscalebias = 0.6

try:
    snr_thresh
except NameError:
    print "S/N threshold defaulting to 4."
    snr_thresh = 4.0

# Look up the scales to use
if array == '7m':
    scales_as_angle = [0, 5, 10]
elif array == '12m':
    scales_as_angle = [0, 1, 2.5, 5]
elif array == '12m+7m':
    scales_as_angle = [0, 1, 2.5, 5, 10]
else:
    print "Array unrecognized. Aborting."
    abort = True

do_image = True
if abort:
    do_image = False

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# Image
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

if do_image:

    print ""
    print "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+"
    print "Imaging the data."
    print "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+"
    print ""

    do_make_dirty_cube = True
    do_revert_to_dirty = False

    if (clean_mask == "None") or (clean_mask == None):
        do_read_in_clean_mask = False
        do_use_pbmask = True
    else:
        do_read_in_clean_mask = True

    do_multiscale_clean = True
    do_singlescale_clean = True
    do_postprocess = True
    
    multiscale_delta_thresh = 0.01
    singlescale_delta_thresh = 0.01

    execfile('../scripts/imageMultiscale.py')
