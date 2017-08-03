import os

only_target = ['ngc1365']
skip_target = []

only_array = ['7m']
skip_array = []

infile = open('../scripts/imaging_key.txt', 'r')

input_vis_list = []
cube_root_list = []
pb_limit_list = []
multiscale_snr_thresh_list = []
smallscalebias_list = []
clean_mask_list = []
singlescale_snr_thresh_list = []

dir_list = []
out_root_list = []
linetag_list = []
array_list = []

while True:
    line  = infile.readline()    
    if len(line) == 0:
        break
    if line[0] == '#':
        continue
    words = line.split()
    if len(words) < 7:
        continue

    this_input_vis = words[0]
    this_cube_root = words[1]
    this_pb_limit = words[2]
    this_multiscale_snr_thresh = words[3]
    this_smallscalebias = words[4]
    this_clean_mask = words[5]
    this_singlescale_snr_thresh = words[6]

    split_cube_root = this_cube_root.split('_')
    this_out_root = split_cube_root[0]
    this_linetag = split_cube_root[1]
    this_array = split_cube_root[2]
    
    target_ok = True
    if len(only_target) > 0:
        target_ok = False
        for this_only_target in only_target:
            if this_out_root == this_only_target:
                target_ok = True

    if len(skip_target) > 0:
        for this_skip_target in skip_target:
            if this_out_root == this_skip_target:
                target_ok = False

    array_ok = True
    if len(only_array) > 0:
        array_ok = False
        for this_only_array in only_array:
            if this_array == this_only_array:
                array_ok = True

    if len(skip_array) > 0:
        for this_skip_array in skip_array:
            if this_array == this_skip_array:
                array_ok = False

    add_this = (target_ok == True) and (array_ok == True)
    
    if add_this:
        out_root_list.append(this_out_root)
        if this_out_root.find('north') != -1:
            dir_list.append(this_out_root[0:this_out_root.find('north')])
        elif this_out_root.find('south') != -1:
            dir_list.append(this_out_root[0:this_out_root.find('south')])
        else:
            dir_list.append(this_out_root)
        linetag_list.append(this_linetag)
        array_list.append(this_array)
        
        input_vis_list.append(this_input_vis)
        cube_root_list.append(this_cube_root)

        pb_limit_list.append(float(this_pb_limit))
        multiscale_snr_thresh_list.append(float(this_multiscale_snr_thresh))
        smallscalebias_list.append(float(this_smallscalebias))
        clean_mask_list.append(this_clean_mask)
        singlescale_snr_thresh_list.append(float(this_singlescale_snr_thresh))

infile.close()

for ii in range(len(dir_list)):
    os.chdir('../'+dir_list[ii])
    out_root = out_root_list[ii]
    
    input_vis = input_vis_list[ii]
    cube_root = cube_root_list[ii]
    linetag = linetag_list[ii]
    array = array_list[ii]
    pb_limit = pb_limit_list[ii]
    smallscalebias = smallscalebias_list[ii]
    multiscale_snr_thresh = multiscale_snr_thresh_list[ii]
    clean_mask = clean_mask_list[ii]
    singlescale_snr_thresh = singlescale_snr_thresh_list[ii]

    execfile('../scripts/phangsImagingPipeline2.py')

os.chdir('../scripts/')
