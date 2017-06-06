# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# PREPARATION AND CONTROL FLOW
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

# --------------------------------------
# User inputs
# --------------------------------------

out_root = 'ngc1566'
tag = '925'
phase_center = 'J2000 04h20m00.4s -54d56m16s'
source_vel_kms = 1504.
vwidth_kms = 500.

calibrated_files = {
    '7m_1':'../../2015.1.00925.S/science_goal.uid___A001_X2fe_X31e/group.uid___A001_X2fe_X31f/member.uid___A001_X2fe_X322/calibrated/calibrated_final.ms',
    '7m_2':'../../2015.1.00925.S/science_goal.uid___A001_X2fe_X314/group.uid___A001_X2fe_X315/member.uid___A001_X2fe_X318/calibrated/calibrated_final.ms'
    }

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
chan_dv_kms = 6.0

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

    #files_to_concat = []
    #for this_tag in calibrated_files.keys():
    #    if this_tag[0:3] == '12m':
    #        files_to_concat.append(out_root+'_'+this_tag+'_co21.ms')
    #out_file = out_root+'_12m_co21.ms'
    #os.system('rm -rf '+out_file)
    #os.system('rm -rf '+out_file+'.flagversions')
    #concat(vis=files_to_concat,
    #       concatvis=out_file)

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
    
    input_vis_7m = 'ngc1566_7m_co21.ms'
    cube_root_7m = 'ngc1566_co21_7m'

    input_vis_combo = 'ngc1566_925_co21.ms'
    cube_root_combo = 'ngc1566_co21_12m'

    input_vis_12m = 'ngc1566_12m_co21.ms'
    cube_root_12m = 'ngc1566_co21_12m'

    do_image_7m = True
    do_image_combo = False
    do_image_12m = False

    execfile('../scripts/phangsImagingPipeline.py')
