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
calibrated_file = ['../../2015.1.00956.S/science_goal.uid___A001_X2fb_X2b7/group.uid___A001_X2fb_X2b8/member.uid___A001_X2fb_X2b9/calibrated/calibrated_final.ms',
                   '../../2015.1.00956.S/science_goal.uid___A001_X2fb_X2c1/group.uid___A001_X2fb_X2c2/member.uid___A001_X2fb_X2c3/calibrated/calibrated_final.ms']
field = ''

# ... parameters for regridding and binning
gal = 'ngc4321'
start_vel = '1400km/s'
nchan = 150
flag_co21 = '2:200~550'
flag_c18o21 = '3:15~50'

# ... target resolution used in post-processing
target_beam_co21 = '1.0arcsec'
target_beam_c18o21 = '1.1arcsec'
target_beam_cont = '1.1arcsec'

# imaging parameters
phase_center = 'J2000 12h22m54.8s +15d49m19s'
imsize = [2048, 2048]
imsize_lowres = [600, 600]

# threshold parameters for cleaning
thresh_factor = 1.5

# --------------------------------------
# Execute Script
# --------------------------------------

execfile('../scripts/image1mmGalaxy.py')
