# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# PREPARATION AND CONTROL FLOW
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

# --------------------------------------
# User inputs
# --------------------------------------

out_root = 'ngc3627south'
tag = '956'
# overall center: phase_center = 'J2000 11h20m14.9s +12d59m30'
phase_center = 'J2000 11h20m14.9s +12d57m48'
source_vel_kms = 727
vwidth_kms = 700

calibrated_files = {'12m':'../../2015.1.00956.S/science_goal.uid___A001_X2fb_X28f/group.uid___A001_X2fb_X290/member.uid___A001_X2fb_X291/calibrated/uid___A002_Xb12f3b_X92ce.ms.split.cal/',
                    '7m':'../../2015.1.00956.S/science_goal.uid___A001_X2fb_X28f/group.uid___A001_X2fb_X290/member.uid___A001_X2fb_X293/calibrated/calibrated_final.ms'}

clean_mask_file = '../clean_masks/ngc3627_co21_clean_mask.fits'

# --------------------------------------
# Overall control flow
# --------------------------------------

execfile('../scripts/line_list.py')

# Extract data
script_copy = False
script_extract_co21 = False
script_extract_c18o21 = False
script_extract_continuum = False

# Image data
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
    
    input_vis_7m = 'ngc3627south_7m_co21.ms'
    cube_root_7m = 'ngc3627south_co21_7m'

    input_vis_combo = 'ngc3627south_956_co21.ms'
    cube_root_combo = 'ngc3627south_co21_12m+7m'

    input_vis_12m = 'ngc3627south_12m_co21.ms'
    cube_root_12m = 'ngc3627south_co21_12m'

    do_image_7m = True
    do_image_combo = False
    do_image_12m = False

    execfile('../scripts/phangsImagingPipeline.py')
