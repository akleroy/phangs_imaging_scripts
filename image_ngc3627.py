# --------------------------------------
# Control Flow
# --------------------------------------

do_copy = False
do_cleanup = False
split_mosaic = True

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
calibrated_file = ['../../2015.1.00956.S/science_goal.uid___A001_X2fb_X28f/group.uid___A001_X2fb_X290/member.uid___A001_X2fb_X291/calibrated/uid___A002_Xb12f3b_X92ce.ms.split.cal/',
                   '../../2015.1.00956.S/science_goal.uid___A001_X2fb_X285/group.uid___A001_X2fb_X286/member.uid___A001_X2fb_X287/calibrated/uid___A002_Xaf5c32_X6d2.ms.split.cal']
field = 'NGC_3627'

# ... parameters for regridding and binningn_
gal = 'ngc3627'
start_vel = '475km/s'
nchan = 200
flag_co21 = '2:200~550'
flag_c18o21 = '3:15~50'

# ... target resolution used in post-processing
target_beam_co21 = '1.4arcsec'
target_beam_c18o21 = '1.5arcsec'
target_beam_cont = '1.5arcsec'

# imaging parameters
phase_center = 'J2000 11h20m14.9s +12d59m30'
#phase_center = 'J2000 11h20m14.9s +12d58m30'
imsize = [1500, 3000]
imsize_lowres = [600, 1000]

# threshold parameters for cleaning
thresh_factor = 1.5

skip_mask = True

# --------------------------------------
# Execute Script
# --------------------------------------

execfile('../scripts/image1mmGalaxy.py')
