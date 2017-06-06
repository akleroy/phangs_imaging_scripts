# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# PREPARATION AND CONTROL FLOW
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

# --------------------------------------
# User inputs
# --------------------------------------

out_root = 'ic5332'
tag = '925'
phase_center = 'J2000 23h34m27.5s -36d06m04s'
source_vel_kms = 701.
vwidth_kms = 250.

calibrated_files = {
    '7m':'../../2015.1.00925.S/science_goal.uid___A001_X2fe_X2ba/group.uid___A001_X2fe_X2bb/member.uid___A001_X2fe_X2be/calibrated/calibrated_final.ms'
    , '12m':'../../2015.1.00925.S/science_goal.uid___A001_X2fe_X2ba/group.uid___A001_X2fe_X2bb/member.uid___A001_X2fe_X2bc/calibrated/calibrated_final.ms'
    }

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
chan_dv_kms = 6.0

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

    do_use_pbmask = True
    linetag = 'co21'
    specmode = 'cube'    
    restfreq_ghz = line_list[linetag]
    max_loop = 10
    pb_limit = 0.25
    uvtaper = None    
    
    input_vis_7m = 'ic5332_7m_co21.ms'
    cube_root_7m = 'ic5332_co21_7m'

    input_vis_combo = 'ic5332_925_co21.ms'
    cube_root_combo = 'ic5332_co21_12m'

    input_vis_12m = 'ic5332_12m_co21.ms'
    cube_root_12m = 'ic5332_co21_12m'

    do_image_7m = True
    do_image_combo = False
    do_image_12m = False

    execfile('../scripts/phangsImagingPipeline.py')
