# --------------------------------------
# Overall control flow
# --------------------------------------

script_copy = False
script_extract_co21 = False
script_extract_c18o21 = False
script_extract_continuum = False
script_image_co21 = True
script_image_c18o21 = False

# --------------------------------------
# Copy the data
# --------------------------------------

out_root = 'ngc4535'
tag = '956'

calibrated_files = {'12m':'../../2015.1.00956.S/science_goal.uid___A001_X2fb_X2cb/group.uid___A001_X2fb_X2cc/member.uid___A001_X2fb_X2cd/calibrated/calibrated_final.ms',
              '7m':'../../2015.1.00956.S/science_goal.uid___A001_X2fb_X2cb/group.uid___A001_X2fb_X2cc/member.uid___A001_X2fb_X2cf/calibrated/calibrated_final.ms'}

fields_to_use = {'12m':'',
                 '7m':''}

if script_copy:
    do_copy = True
    do_concat = True
    do_extract = False
    execfile('../scripts/extractLineData.py')

# --------------------------------------
# Extract line data
# --------------------------------------

source_vel_kms = 1970
vwidth_kms = 500

# 12CO 2-1
linetag = 'co21'
restfreq_ghz = 230.53800
chan_dv_kms = 2.5
nchan = 240

if script_extract_co21:
    do_copy = False
    do_concat = False
    do_extract = True
    execfile('../scripts/extractLineData.py')

# C18O 2-1
linetag = 'c18o21'
restfreq_ghz = 219.56035
chan_dv_kms = 5.0

if script_extract_c18o21:
    do_copy = False
    do_concat = False
    do_extract = True    
    execfile('../scripts/extractLineData.py')

# --------------------------------------
# Extract continuum data
# --------------------------------------

if script_extract_continuum:
    pass

# --------------------------------------
# Image the data
# --------------------------------------

phase_center = 'J2000 12h34m20.3s +08d11m52'

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
