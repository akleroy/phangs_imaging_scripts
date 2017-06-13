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
    do_multiscale_clean = True
    do_revert_to_multiscale = False
    do_singlescale_clean = True
    do_postprocess = True
    
    smallscalebias = 0.6
    outerscale = 15.

    snr_thresh = 3.0
    multiscale_delta_thresh = 0.02
    singlescale_delta_thresh = 0.005
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
    
    imstat_7m = imstat(cube_root_7m+'.residual')
    multiscale_threshold = str(4.0*(imstat_7m['medabsdevmed'][0]/0.6745))+'Jy/beam'
    snr_thresh = 4.0
    singlescale_threshold = None
   
    smallscalebias = 0.6
    outerscale = 12.0
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

    multiscale_threshold = str(4.0*(imstat_7m['medabsdevmed'][0]/0.6745))+'Jy/beam'
    outerscale = 6.2
    smallscalebias = 0.6
    snr_thresh = 4.0
    singlescale_threshold = None

    execfile('../scripts/imageMultiscale.py')

