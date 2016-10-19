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
calibrated_file = ['../../2015.1.00925.S/science_goal.uid___A001_X2fe_X328/group.uid___A001_X2fe_X329/member.uid___A001_X2fe_X32a/calibrated/uid___A002_Xb148a2_X255e.ms.split.cal',
                   '../../2015.1.00925.S/science_goal.uid___A001_X2fe_X332/group.uid___A001_X2fe_X333/member.uid___A001_X2fe_X334/calibrated/uid___A002_Xb148a2_X2c49.ms.split.cal']
field = 'NGC_5068'

# ... parameters for regridding and binning
gal = 'ngc5068'
co21_spw = '0'
c18o_spw = '2'
c18o_deltav = '6km/s'
start_vel = '545km/s'
nchan = 100
flag_co21 = '0:200~550'
flag_c18o21 = '2:15~50'

# ... target resolution used in post-processing
target_beam_co21 = '1.0arcsec'
target_beam_c18o21 = '1.1arcsec'
target_beam_cont = '1.1arcsec'

# imaging parameters
phase_center = 'J2000 13h18m54.8s -21d02m21s'
imsize = [2500, 2500]
imsize_lowres = [800, 800]

# threshold parameters for cleaning
thresh_factor = 1.5

# --------------------------------------
# Execute Script
# --------------------------------------

execfile('../scripts/image1mmGalaxy.py')
