# --------------------------------------
# Check inputs
# --------------------------------------

try:
    out_root
except NameError:
    print "Please specify out_root= a known file name."
    out_root = 'None'

# --------------------------------------
# Read galaxy header file
# --------------------------------------

infile = open('../scripts/ms_file_key.txt', 'r')

tag = None
while True:
    line  = infile.readline()    
    if len(line) == 0:
        break
    if line[0] == '#':
        continue
    words = line.split()
    if len(words) < 4:
        continue

    this_out_root = words[0]
    this_tag = words[1]
    this_key = words[2]
    this_ms = words[3]
    
    if this_out_root == out_root:
        if tag == None:
            tag = this_tag
            calibrated_files = {}
        calibrated_files[this_key] = this_ms

infile.close()

# --------------------------------------
# Load the line list
# --------------------------------------

execfile('../scripts/line_list.py')

# --------------------------------------
# Look up which steps to perform
# --------------------------------------

infile = open('../scripts/extraction_key.txt', 'r')

script_copy = False
script_contsub = False
script_extract_co21 = False
script_extract_c18o21 = False
script_extract_continuum = False
deltav_co21 = 2.5
deltav_c18o21 = 6.0

while True:
    line  = infile.readline()    
    if len(line) == 0:
        break
    if line[0] == '#':
        continue
    words = line.split()
    if len(words) < 8:
        continue

    this_out_root = words[0]
    this_copy = words[1]
    this_contsub = words[2]
    this_extract_co21 = words[3]
    this_extract_c18o21 = words[4]
    this_extract_continuum = words[5]
    this_deltav_co21 = float(words[6])
    this_deltav_c18o21 = float(words[7])
    
    if this_out_root == out_root:
        script_copy = (this_copy == 'True')
        script_contsub = (this_contsub == 'True')
        script_extract_co21 = (this_extract_co21 == 'True')
        script_extract_c18o21 = (this_extract_c18o21 == 'True')
        script_extract_continuum = (this_extract_continuum == 'True')
        deltav_co21 = this_deltav_co21
        deltav_c18o21 = this_deltav_c18o21

infile.close()

infile = open('../scripts/mosaic_definitions.txt', 'r')

tag = None
while True:
    line  = infile.readline()    
    if len(line) == 0:
        break
    if line[0] == '#':
        continue
    words = line.split()
    if len(words) < 5:
        continue

    this_out_root = words[0]
    this_ra = words[1]
    this_dec = words[2]
    this_vsys = words[3]
    this_deltav = words[4]
    
    if this_out_root == out_root:
        source_vel_kms = float(this_vsys)
        vwidth_kms = float(this_deltav)

infile.close()

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
# Option for conitnuum subtraction
# --------------------------------------

if script_contsub:
    pass

# --------------------------------------
# Extract line data
# --------------------------------------

if script_extract_co21:
    linetag = 'co21'
    restfreq_ghz = line_list[linetag]
    chan_dv_kms = deltav_co21

    do_copy = False
    do_split = False
    do_extract = True
    do_combine = True
    execfile('../scripts/extractLineData.py')

if script_extract_c18o21:
    linetag = 'c18o21'
    restfreq_ghz = line_list[linetag]
    chan_dv_kms = deltav_c18o21

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
