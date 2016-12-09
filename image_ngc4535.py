# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# PREPARATION AND CONTROL FLOW
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

# --------------------------------------
# User inputs
# --------------------------------------

out_root = 'ngc4535'
tag = '956'
#phase_center = 'J2000 12h34m20.3s +08d11m52'
source_vel_kms = 1970
vwidth_kms = 500

calibrated_files = {'12m':'../../2015.1.00956.S/science_goal.uid___A001_X2fb_X2cb/group.uid___A001_X2fb_X2cc/member.uid___A001_X2fb_X2cd/calibrated/calibrated_final.ms',
              '7m':'../../2015.1.00956.S/science_goal.uid___A001_X2fb_X2cb/group.uid___A001_X2fb_X2cc/member.uid___A001_X2fb_X2cf/calibrated/calibrated_final.ms'}

# --------------------------------------
# Overall control flow
# --------------------------------------

execfile('../scripts/line_list.py')

# Extract data
script_copy = True
script_extract_co21 = True
script_extract_c18o21 = True
script_extract_continuum = True

# Image data
script_image_co21 = False
script_image_c18o21 = False

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
# IMAGING AND EXPORT
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

# --------------------------------------
# Image channel 0
# --------------------------------------

if script_image_chan0:
    do_pickcellsize = True
    do_init = False
    do_mask = False
    do_clean = False
    do_postprocess = False

    input_vis = 'ngc4535_956_co21_chan0.ms'
    cube_root = 'ngc4535_co21_tapercube'
    uvtaper = ["2.5arcsec", "2.5arcsec", "0.0deg"]
    force_cell_size = 0.5
    niter = 1E6
    execfile('../scripts/imageImage.py')    

# --------------------------------------
# Image cubes
# --------------------------------------

if script_image_co21:
    do_init = True
    do_mask = False
    do_clean = False
    do_postprocess = False

    im_size = [600, 600]
    cell_size = '0.5arcsec'
    input_vis = 'ngc4535_co21_regrid.ms'
    cube_root = 'ngc4535_co21_tapercube'
    uvtaper = ["2.5arcsec", "2.5arcsec", "0.0deg"]
    niter = 1E6
    execfile('../scripts/imageLine.py')

if script_image_c18o21:
    pass
