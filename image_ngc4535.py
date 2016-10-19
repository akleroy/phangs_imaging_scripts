# --------------------------------------
# Control Flow
# --------------------------------------

do_copy = True
do_cleanup = False
split_mosaic = False
has_7m = True

do_mask = False
do_co21 = False
do_c18o = False
do_cont = False

do_native = False
do_taper = True

# --------------------------------------
# Inputs
# --------------------------------------

# ... original data
calibrated_file = '../../2015.1.00956.S/science_goal.uid___A001_X2fb_X2cb/group.uid___A001_X2fb_X2cc/member.uid___A001_X2fb_X2cd/calibrated/calibrated_final.ms'
calibrated_7m_file = '../../2015.1.00956.S/science_goal.uid___A001_X2fb_X2cb/group.uid___A001_X2fb_X2cc/member.uid___A001_X2fb_X2cf/calibrated/calibrated_final.ms'

# ... SPW and field labeling
field = ''
field_7m = ''
co21_spw_7m = '0,4,8,12'
c18o_spw_7m = '2,6,10,14'

flag_co21 = '2:200~550'
flag_c18o21 = '3:15~50'
flag_co21_7m = '0:200~550'
flag_c18o21_7m = '2:15~50'

# ... parameters for regridding and binning
gal = 'ngc4535'
start_vel = '1800km/s'
nchan = 120

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

# --------------------------------------
# Execute Script
# --------------------------------------

execfile('../scripts/image1mmGalaxy.py')
