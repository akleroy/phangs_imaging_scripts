# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# PREPARATION AND CONTROL FLOW
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

# --------------------------------------
# User inputs
# --------------------------------------

out_root = 'ngc0628'
tag = '650'
phase_center = 'J2000 01h36m41.7s +15d47m01'
source_vel_kms = 657
vwidth_kms = 200

dir_12m = '../../Cycle1_650/2012.1.00650.S/science_goal.uid___A002_X5a9a13_X579/group.uid___A002_X5a9a13_X57a/member.uid___A002_X5a9a13_X57b/calibrated/'

dir_7m = '../../Cycle1_650/2012.1.00650.S/science_goal.uid___A002_X5a9a13_X579/group.uid___A002_X5a9a13_X57a/member.uid___A002_X5a9a13_X57d/calibrated/'

calibrated_files = {'12m_1':dir_12m+'uid___A002_X5b2f01_X3f.ms.split.cal.M74/',
                    '12m_2':dir_12m+'uid___A002_X8081ba_X4018.ms.split.cal.M74/',
                    '12m_3':dir_12m+'uid___A002_X8081ba_X4527.ms.split.cal.M74/',
                    '12m_4':dir_12m+'uid___A002_X80c782_X1911.ms.split.cal.M74/',
                    '12m_5':dir_12m+'uid___A002_X8204db_X611.ms.split.cal.M74/',
                    '12m_6':dir_12m+'uid___A002_X95b353_X471.ms.split.cal.M74/',
                    '12m_7':dir_12m+'uid___A002_X960614_X2d5b.ms.split.cal.M74/',
                    '12m_8':dir_12m+'uid___A002_X966cea_X96c.ms.split.cal.M74/',
                    '7m_1':dir_7m+'uid___A002_X6f1341_X80b.ms.split.cal',
                    '7m_2':dir_7m+'uid___A002_X6f2c6e_Xb0a.ms.split.cal',
                    '7m_3':dir_7m+'uid___A002_X7fc9da_X1f0d.ms.split.cal',
                    '7m_4':dir_7m+'uid___A002_X7fc9da_X4b45.ms.split.cal',
                    '7m_5':dir_7m+'uid___A002_X8081ba_X11fb.ms.split.cal',
                    '7m_6':dir_7m+'uid___A002_X8081ba_X3e04.ms.split.cal',
                    '7m_7':dir_7m+'uid___A002_X8081ba_X44b8.ms.split.cal',
                    '7m_8':dir_7m+'uid___A002_X8081ba_X85e.ms.split.cal',
                    '7m_9':dir_7m+'uid___A002_X8081ba_Xce1.ms.split.cal',
                    '7m_10':dir_7m+'uid___A002_X8204db_X4f.ms.split.cal',
                    }

clean_mask_file = '../clean_masks/ngc0628_co21_widemask.fits'

# --------------------------------------
# Overall control flow
# --------------------------------------

execfile('../scripts/line_list.py')

# Extract data
script_copy = False
script_extract_co21 = False
script_extract_c18o21 = False
script_extract_continuum = False

# Image data
script_image_chan0 = False
script_image_cube = True

script_image_co21 = True
script_image_c18o21 = False
script_image_cont = False

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
# IMAGING
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

# --------------------------------------
# Image channel 0
# --------------------------------------

if script_image_chan0:
    
    if script_image_co21:
        do_end_to_end = True
        do_start_with_pbmask = True
        input_vis = out_root+'_'+tag+'_co21_chan0.ms'
        cube_root = out_root+'_co21_chan0'
        uvtaper = None
        linetag = 'co21'
        specmode = 'mfs'
        restfreq_ghz = line_list[linetag]
        execfile('../scripts/imageImage.py')

    if script_image_c18o21:
        do_end_to_end = True
        do_start_with_pbmask = True
        input_vis = out_root+'_'+tag+'_c18o21_chan0.ms'
        cube_root = out_root+'_c18o21_chan0'
        uvtaper = None
        linetag = 'c18o21'
        specmode = 'mfs'
        restfreq_ghz = line_list[linetag]
        execfile('../scripts/imageImage.py')

    if script_image_cont:
        do_end_to_end = True
        do_start_with_pbmask = True
        input_vis = out_root+'_'+tag+'_cont.ms'
        cube_root = out_root+'_cont'
        uvtaper = None
        specmode = 'mfs'
        restfreq_ghz = ''
        execfile('../scripts/imageImage.py')

# --------------------------------------
# Image cubes
# --------------------------------------

if script_image_cube:

    if script_image_co21:
        do_end_to_end = False
        do_revert_to_dirty = True
        do_read_in_clean_mask = True
        do_clean = True
       
        input_vis = out_root+'_'+tag+'_co21.ms'
        cube_root = out_root+'_co21'
        uvtaper = None
        linetag = 'co21'
        specmode = 'cube'
        scales_to_use=[0,2,4,8,16,32,64]
        restfreq_ghz = line_list[linetag]

        execfile('../scripts/imageMultiscale.py')

    if script_image_c18o21:
        do_end_to_end = True
        do_start_with_pbmask = False
        
        input_vis = out_root+'_'+tag+'_c18o21.ms'
        cube_root = out_root+'_c18o21'
        uvtaper = None
        linetag = 'c18o21'
        specmode = 'cube'
        restfreq_ghz = line_list[linetag]

        execfile('../scripts/imageImage.py')
