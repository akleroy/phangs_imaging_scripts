tested_versions = ['4.6.0','4.7.0','4.7.1','4.7.2', '5.0.0']
this_version = (casa['build']['version']).split('-')[0]
if this_version not in tested_versions:
    print "The script hasn't been verified for this version of CASA."
    print "This version of CASA is "+this_version
    print "Tested versions are "+str(tested_versions)

execfile('../scripts/auExtensions.py')

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# Inputs
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

print "****************************************"
print "********** phangsImagingPipeline.py ****"
print "****************************************"

abort = False

try:
    input_vis
except NameError:
    print "Specify an input visibility. Aborting."
    abort = True

try:
    cube_root
except NameError:
    print "Specify a cube root. Aborting."
    abort = True

try:
    array
except NameError:
    print "Specify the array being used. Aborting."
    array = None
    abort = True

try:
    clean_mask
except NameError:
    print "No clean mask. Will use a primary beam based mask."
    clean_mask = "None"

try:
    pb_limit
except NameError:
    print "Primary beam coverage defaulting to 0.25."
    pb_limit = 0.25

try:
    smallscalebias
except NameError:
    print "Small scale bias defaulting to 0.6."
    smallscalebias = 0.6

try:
    multiscale_snr_thresh
except NameError:
    print "Multiscale S/N threshold defaulting to 4."
    multiscale_snr_thresh = 4.0

try:
    singlescale_snr_thresh
except NameError:
    print "Single scale S/N threshold defaulting to 4."
    singlescale_snr_thresh = 4.0

# Look up the scales to use
if array == '7m':
    scales_as_angle = [0, 5, 10]
elif array == '12m':
    scales_as_angle = [0, 1, 2.5, 5]
elif array == '12m+7m':
    scales_as_angle = [0, 1, 2.5, 5, 10]
else:
    print "Array unrecognized. Aborting."
    abort = True

do_image = True
if abort:
    do_image = False

specmode = 'cube'
max_loop = 20
uvtaper = None

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# Look up the source
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

infile = open('../scripts/mosaic_definitions.txt', 'r')

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
        phase_center = 'J2000 '+this_ra+' '+this_dec
        source_vel_kms = float(this_vsys)
        vwidth_kms = float(this_deltav)

infile.close()

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# Image
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

if do_image:

    print ""
    print "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+"
    print "Imaging the data."
    print "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+"
    print ""

    do_make_dirty_cube = True
    do_revert_to_dirty = False

    if (clean_mask == "None") or (clean_mask == None):
        do_read_in_clean_mask = False
        do_use_pbmask = True
    else:
        do_read_in_clean_mask = True
        clean_mask_file = '../clean_masks/'+clean_mask

    do_multiscale_clean = True
    do_singlescale_clean = True
    do_postprocess = True
    
    multiscale_delta_thresh = 0.01
    singlescale_delta_thresh = 0.01

    execfile('../scripts/imageMultiscale2.py')
