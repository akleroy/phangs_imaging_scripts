# Note that the NGC 0628 setup not cover C18O emission.

# --------------------------------------
# Control Flow
# --------------------------------------

do_copy = True
do_cleanup = False
split_mosaic = True
has_c18o = False

do_mask = True
do_co21 = True
do_c18o = True
do_cont = True

do_native = True
do_taper = True

# --------------------------------------
# Inputs
# --------------------------------------
calibrated_dir = '../../../Cycle1_650/2012.1.00650.S/science_goal.uid___A002_X5a9a13_X579/group.uid___A002_X5a9a13_X57a/member.uid___A002_X5a9a13_X57b/calibrated/'

co21_spw = '0'

# ... original data
calibrated_file = [
    calibrated_dir + 'uid___A002_X5b2f01_X3f.ms.split.cal.M74',
    calibrated_dir + 'uid___A002_X8081ba_X4018.ms.split.cal.M74',
    calibrated_dir + 'uid___A002_X8081ba_X4527.ms.split.cal.M74',
    calibrated_dir + 'uid___A002_X80c782_X1911.ms.split.cal.M74',
    calibrated_dir + 'uid___A002_X8204db_X611.ms.split.cal.M74',
    calibrated_dir + 'uid___A002_X95b353_X471.ms.split.cal.M74',
    calibrated_dir + 'uid___A002_X960614_X2d5b.ms.split.cal.M74',
    calibrated_dir + 'uid___A002_X966cea_X96c.ms.split.cal.M74'
    ]
field = ''

# ... parameters for regridding and binning
gal = 'ngc0628'
start_vel = '540km/s'
nchan = 80
flag_co21 = '0:200~550'
flag_c18o21 = '0:200~550'

# ... target resolution used in post-processing
target_beam_co21 = '1.0arcsec'
target_beam_c18o21 = '1.0arcsec'
target_beam_cont = '1.0arcsec'

# imaging parameters
phase_center = 'J2000 01h36m41.7s +15d47m01'
imsize = [2500, 2500]
imsize_lowres = [800, 800]

# threshold parameters for cleaning
thresh_factor = 1.5

# --------------------------------------
# Execute Script
# --------------------------------------

execfile('../scripts/image1mmGalaxy.py')
