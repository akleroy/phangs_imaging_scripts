#############
### SOUTH ###
#############

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
calibrated_file = '../../2015.1.00956.S/science_goal.uid___A001_X2fb_X2df/group.uid___A001_X2fb_X2e0/member.uid___A001_X2fb_X2e1/calibrated/calibrated_final.ms'
field = ''

# ... parameters for regridding and binning
gal = 'ngc6744south'
start_vel = '650km/s'
nchan = 150
flag_co21 = '2:225~525'
flag_c18o21 = '3:15~50'

# ... target resolution used in post-processing
target_beam_co21 = '1.0arcsec'
target_beam_c18o21 = '1.1arcsec'
target_beam_cont = '1.1arcsec'

# imaging parameters
phase_center = 'J2000 19h09m46.1s -63d53m21.8'
imsize = [2000,1500]
imsize_lowres = [600, 600]

# factor to multiply by RMS beam for clean threshold
thresh_factor = 1.5

# --------------------------------------
# Execute Script
# --------------------------------------

execfile('../scripts/image1mmGalaxy.py')
