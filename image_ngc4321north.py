# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# PREPARATION AND CONTROL FLOW
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

# --------------------------------------
# User inputs
# --------------------------------------

out_root = 'ngc4321north'
tag = '956'
phase_center = 'J2000 12h22m54.8s +15d50m05s'
source_vel_kms = 1571
vwidth_kms = 500

calibrated_files = {'12m':'../../2015.1.00956.S/science_goal.uid___A001_X2fb_X2b7/group.uid___A001_X2fb_X2b8/member.uid___A001_X2fb_X2b9/calibrated/calibrated_final.ms',
                    '7m':'../../2015.1.00956.S/science_goal.uid___A001_X2fb_X2b7/group.uid___A001_X2fb_X2b8/member.uid___A001_X2fb_X2bb/calibrated/calibrated_final.ms',
                    }

clean_mask_file = '../clean_masks/ngc4321_co21_clean_mask.fits'

# --------------------------------------
# Overall control flow
# --------------------------------------

execfile('../scripts/line_list.py')

# Extract data
script_copy = False
script_extract_co21 = False
script_extract_c18o21 = False
script_extract_continuum = False

# Image data>
script_image_cube = True

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# EXTRACTION
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

# --------------------------------------
# Copy the data
# --------------------------------------

if script_copy:
    do_copy = True
    do_split = True
    do_extract = False
    do_combine = False
    execfile('../scripts/extractLineData.py')

# --------------------------------------
# Extract line data
# --------------------------------------

# 12CO 2-1
linetag = 'co21'
restfreq_ghz = line_list[linetag]
chan_dv_kms = 2.5

if script_extract_co21:
    do_copy = False
    do_split = False
    do_extract = True
    do_combine = True
    execfile('../scripts/extractLineData.py')

# C18O 2-1
linetag = 'c18o21'
restfreq_ghz = line_list[linetag]
chan_dv_kms = 5.0

if script_extract_c18o21:
    do_copy = False
    do_split = False
    do_extract = True
    do_combine = True
    execfile('../scripts/extractLineData.py')

# --------------------------------------
# Extract continuum data
# --------------------------------------

if script_extract_continuum:
    do_recopy = True
    do_flag = True
    do_average = True
    do_statwt = True
    lines_to_flag = lines_co+lines_13co+lines_c18o
    execfile('../scripts/extractContinuum.py')

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# IMAGING
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

if script_image_cube:

    do_use_pbmask = False
    linetag = 'co21'
    specmode = 'cube'    
    restfreq_ghz = line_list[linetag]
    max_loop = 20
    pb_limit = 0.25
    uvtaper = None    
    
    input_vis_7m = 'ngc4321north_7m_co21.ms'
    cube_root_7m = 'ngc4321north_co21_7m'

    input_vis_combo = 'ngc4321north_956_co21.ms'
    cube_root_combo = 'ngc4321north_co21_12m+7m'

    input_vis_12m = 'ngc4321north_12m_co21.ms'
    cube_root_12m = 'ngc4321north_co21_12m'

    do_image_7m = False
    do_image_combo = True
    do_image_12m = True

    smallscalebias_7m = 0.6
    smallscalebias_combo = 0.8
    smallscalebias_12m = 0.8

    execfile('../scripts/phangsImagingPipeline.py')
