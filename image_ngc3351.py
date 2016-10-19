# --------------------------------------
# Control Flow
# --------------------------------------

do_copy = False
do_cleanup = False
split_mosaic = False

do_mask = True
do_co21 = True
do_c18o = True
do_cont = True

do_native = True
do_taper = True

# --------------------------------------
# Inputs
# --------------------------------------

# ... original data
calibrated_file = '../../2015.1.00956.S/science_goal.uid___A001_X2fb_X27b/group.uid___A001_X2fb_X27c/member.uid___A001_X2fb_X27d/calibrated/calibrated_final.ms'
field = ''

# ... parameters for regridding and binning
gal = 'ngc3351'
start_vel = '550km/s'
nchan = 160
flag_co21 = '2:225~525'
flag_c18o21 = '3:15~50'

# ... target resolution used in post-processing
target_beam_co21 = '1.3arcsec'
target_beam_c18o21 = '1.4arcsec'
target_beam_cont = '1.4arcsec'

# imaging parameters
phase_center = 'J2000 10h43m57s +11d42m14s'
imsize = [2048, 2048]
imsize_lowres = [600, 600]

# factor to multiply by RMS beam for clean threshold
thresh_factor = 1.5

# --------------------------------------
# Execute Script
# --------------------------------------

execfile('../scripts/image1mmGalaxy.py')
