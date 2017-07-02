tested_versions = ['4.6.0','4.7.0','4.7.1','4.7.2']
this_version = (casa['build']['version']).split('-')[0]
if this_version not in tested_versions:
    print "The script hasn't been verified for this version of CASA."
    print "This version of CASA is "+this_version
    print "Tested versions are "+str(tested_versions)

execfile('../scripts/auExtensions.py')

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# Inputs
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

print ""
print "---------- phangsImagingPipeline.py ----------"
print ""

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# Image the ACA 7m only
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

if do_image_7m:

    print ""
    print "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+"
    print "Imaging the 7m data by itself."
    print "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+"
    print ""

    input_vis = input_vis_7m
    cube_root = cube_root_7m

    do_make_dirty_cube = True
    do_revert_to_dirty = False

    if do_use_pbmask == False:
        do_read_in_clean_mask = True
    else:
        do_read_in_clean_mask = False

    try:
        do_multiscale_clean
    except:
        do_multiscale_clean = True

    do_singlescale_clean = True
    do_postprocess = True
    
    try:
        smallscalebias_7m
    except NameError:
        smallscalebias = 0.6
    else:
        smallscalebias = smallscalebias_7m

    try:
        snr_thresh
    except NameError:
        snr_thresh = 0.0

    try:
        multiscale_delta_thresh
    except NameError:
        multiscale_delta_thresh = 0.01

    try:
        singlescale_delta_thresh
    except NameError:
        singlescale_delta_thresh = 0.01

    scales_as_angle = [0, 5, 10]

    execfile('../scripts/imageMultiscale.py')

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# Image the combined data set
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

if do_image_combo:

    print ""
    print "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+"
    print "Imaging the 7m+12m data together."
    print "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+"
    print ""

    input_vis = input_vis_combo
    cube_root = cube_root_combo
    
    do_make_dirty_cube = True
    do_revert_to_dirty = False

    if do_use_pbmask == False:
        do_read_in_clean_mask = True
    else:
        do_read_in_clean_mask = False

    try:
        do_multiscale_clean
    except:
        do_multiscale_clean = True

    do_singlescale_clean = True
    do_postprocess = True

    try:
        smallscalebias_combo
    except NameError:
        smallscalebias = 0.6
    else:
        smallscalebias = smallscalebias_combo

    try:
        snr_thresh
    except NameError:
        snr_thresh = 0.0

    try:
        multiscale_delta_thresh
    except NameError:
        multiscale_delta_thresh = 0.01

    try:
        singlescale_delta_thresh
    except NameError:
        singlescale_delta_thresh = 0.01

    scales_as_angle = [0, 1, 2.5, 5, 10]

    execfile('../scripts/imageMultiscale.py')   

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# Image the 12-m only data set
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

if do_image_12m:

    print ""
    print "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+"
    print "Imaging the 12m data by itself."
    print "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+"
    print ""

    input_vis = input_vis_12m
    cube_root = cube_root_12m

    do_make_dirty_cube = True
    do_revert_to_dirty = False

    if do_use_pbmask == False:
        do_read_in_clean_mask = True
    else:
        do_read_in_clean_mask = False

    try:
        do_multiscale_clean
    except:
        do_multiscale_clean = True

    do_singlescale_clean = True
    do_postprocess = True

    try:
        smallscalebias_12m
    except NameError:
        smallscalebias = 0.6
    else:
        smallscalebias = smallscalebias_12m

    try:
        snr_thresh
    except NameError:
        snr_thresh = 0.0

    try:
        multiscale_delta_thresh
    except NameError:
        multiscale_delta_thresh = 0.01

    try:
        singlescale_delta_thresh
    except NameError:
        singlescale_delta_thresh = 0.01

    scales_as_angle = [0, 1, 2.5, 5]

    execfile('../scripts/imageMultiscale.py')

