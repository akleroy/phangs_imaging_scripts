# --------------------------------------
# Overall controll flow
# --------------------------------------

script_copy = False
script_extract_co21 = False
script_extract_c18o21 = False
script_extract_cont = False
script_image_co21 = False
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

start_vel = '1600km/s'

# 12CO 2-1
linetag = 'co21'
restfreq = '230.53800GHz'
chan_dv = '2.5km/s'
nchan = 240
chanbin = 4
target_freq_hz = 230.538e9*(1.-1900./2.99e5)

if script_extract_co21:
    do_copy = False
    do_concat = False
    do_extract = True
    execfile('../scripts/extractLineData.py')

# C18O 2-1
linetag = 'c18o21'
restfreq = '219.56035GHz'
chan_dv = '5.0km/s'
nchan = 120
chanbin = 1
target_freq_hz = 219.560e9*(1.-1900./2.99e5)

if script_extract_c18o21:
    do_copy = False
    do_concat = False
    do_extract = True    
    execfile('../scripts/extractLineData.py')

# --------------------------------------
# Extract continuum data
# --------------------------------------

if script_extract_continuum:
    # ... flagging covering the line in the continuum image
    flag_co21 = '2:200~550'
    flag_c18o21 = '3:15~50'
    flag_co21_7m = '0:200~550'
    flag_c18o21_7m = '2:15~50'

# --------------------------------------
# Image the data
# --------------------------------------

if script_image_co21:
    pass

if script_image_co21:
    pass

# ... target resolution used in post-processing
target_beam_co21 = '1.4arcsec'
target_beam_c18o21 = '1.5arcsec'
target_beam_cont = '1.5arcsec'

# imaging parameters
phase_center = 'J2000 12h34m20.3s +08d11m52'
imsize = [2048, 2048]
imsize_lowres = [600, 600]

# threshold parameters for cleaning
thresh_factor = 1.5
