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
special_concat = False

# Image data
script_image_cube = True

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


# --------------------------------------
# Some special concatenation
# --------------------------------------

if special_concat:
    files_to_concat = []
    for this_tag in calibrated_files.keys():
        if this_tag[0:2] == '7m':
            files_to_concat.append(out_root+'_'+this_tag+'_co21.ms')
    out_file = out_root+'_7m_co21.ms'
    os.system('rm -rf '+out_file)
    os.system('rm -rf '+out_file+'.flagversions')
    concat(vis=files_to_concat,
           concatvis=out_file)

    files_to_concat = []
    for this_tag in calibrated_files.keys():
        if this_tag[0:3] == '12m':
            files_to_concat.append(out_root+'_'+this_tag+'_co21.ms')
    out_file = out_root+'_12m_co21.ms'
    os.system('rm -rf '+out_file)
    os.system('rm -rf '+out_file+'.flagversions')
    concat(vis=files_to_concat,
           concatvis=out_file)

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# IMAGING
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

if script_image_cube:

    do_use_pbmask = True
    linetag = 'co21'
    specmode = 'cube'    
    restfreq_ghz = line_list[linetag]
    max_loop = 10
    pb_limit = 0.25
    uvtaper = None    
    
    input_vis_7m = 'ngc0628_7m_co21.ms'
    cube_root_7m = 'ngc0628_co21_7m'

    input_vis_combo = 'ngc0628_650_co21.ms'
    cube_root_combo = 'ngc0628_co21_12m'

    input_vis_12m = 'ngc0628_12m_co21.ms'
    cube_root_12m = 'ngc0628_co21_12m'

    do_image_7m = True
    do_image_combo = False
    do_image_12m = False

    execfile('../scripts/phangsImagingPipeline.py')
