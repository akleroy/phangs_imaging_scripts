execfile('../scripts/header_ngc1566.py')
execfile('../scripts/line_list.py')

script_copy = True
script_extract_co21 = True
script_extract_c18o21 = True
script_extract_continuum = True
special_concat = True

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
