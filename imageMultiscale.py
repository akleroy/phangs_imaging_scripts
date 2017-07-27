# Note that "IMPROVEMENT TBD" items are scattered throughout the
# program. Search on this to find areas where some clean up or a next
# logical step could be worked out.

tested_versions = ['4.6.0','4.7.0','4.7.1','4.7.2']
this_version = (casa['build']['version']).split('-')[0]
if this_version not in tested_versions:
    print "The script hasn't been verified for this version of CASA."
    print "This version of CASA is "+this_version
    print "Tested versions are "+str(tested_versions)

execfile('../scripts/auExtensions.py')

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# TIMER
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

import time
total_start_time = time.time()

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# DEFINITIONS AND INPUTS
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

# ----------------------------------------------------
# CONTROL FLOW
# ----------------------------------------------------

try:
    do_end_to_end
except NameError:
    do_end_to_end = False

try:
    do_use_pbmask
except NameError:
    do_use_pbmask = False
    do_read_in_cleam_mask = False

if do_end_to_end:
    do_make_dirty_cube = True
    do_revert_to_dirty = False
    if do_use_pbmask == False:
        do_read_in_clean_mask = True
    do_multiscale_clean = True
    do_revert_to_multiscale = False
    do_singlescale_clean = True
    do_postprocess = True

try:
    do_pickcellsize
except NameError:
    do_pickcellsize = True

try:
    do_make_dirty_cube
except NameError:
    do_make_dirty_cube = True

try:
    do_revert_to_dirty
except NameError:
    do_revert_to_dirty = False

try:
    do_read_in_clean_mask 
except NameError:
    do_read_in_clean_mask = False

try:
    do_multiscale_clean
except NameError:
    do_multiscale_clean = True

try:
    do_revert_to_multiscale
except NameError:
    do_revert_to_multiscale = False

try:
    do_singlescale
except NameError:
    do_singlescale = False

try:
    do_postprocess
except NameError:
    do_postprocess = False

# ----------------------------------------------------
# TUNING PARAMETERS
# ----------------------------------------------------

try:
    pb_limit
except NameError:
    pb_limit = 0.25

try:
    scales_to_use = [0,2,4,8,16,32,64]
except NameError:
    scales_to_use = [0,2,4,8,16,32,64]

try:
    snr_thresh 
except NameError:
    snr_thresh = 4.0

try:
    max_loop
except NameError:
    max_loop = 10

try:
    single_scale_niter_per_chan
except NameError:
    single_scale_niter_per_chan = 1000

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# Inputs
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

abort = False

print ""
print "---------- imageMultiscale.py ----------"
print ""

# ......................................
# Input/output files
# ......................................

try:
    input_vis
except NameError:
    print "Please define an input visibility via input_vis=XXX."
    abort = True

try:
    cube_root
except NameError:
    print "Please define a root output name for the cube via cube_root=XXX."
    abort = True

# ......................................
# Other option
# ......................................

if abort:
    print "(Turning off other parts of the script)."
    do_pickcellsize = False
    do_make_dirty_cube = False
    do_revert_to_dirty = False
    do_read_in_clean_mask = False
    do_use_pbmask = False
    do_multiscale_clean = False
    do_revert_to_multiscale = False
    do_singlescale_clean = False
    do_postprocess = False
 
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# LOOK UP THE INPUT PARAMETERS FOR CLEAN
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

if do_pickcellsize:
    
    infile = open('../scripts/image_params.txt', 'r')

    vis_list = []
    cell_size_list = []
    x_size_list = []
    y_size_list = []

    while True:
        line  = infile.readline()    
        if len(line) == 0:
            break
        if line[0] == '#':
            continue
        words = line.split()
        if len(words) < 3:
            continue
        vis_list.append(words[0])
        cell_size_list.append(words[1])
        x_size_list.append(words[2])
        y_size_list.append(words[3])

    infile.close()

    for ii in range(len(vis_list)):
        if vis_list[ii] == input_vis:
            cell_size_string = cell_size_list[ii]
            cell_size = float(cell_size_string[0:len(cell_size_string)-len('arcsec')])
            image_size = [int(x_size_list[ii]), int(y_size_list[ii])]

    print "I propose the following:"
    print "... cell size = "+cell_size_string+" or "+str(cell_size)+" arcsec"
    print "... image size = ", image_size

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# MAKE A DIRTY CUBE
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

# Now we make a dirty cube, with 0 iterations supplied. This gives us
# a target astrometry and WCS to work with for external clean
# masks. It also provides a useful point of comparison for checking
# residual rescaling, etc. Save the results separately so that one can
# revert to the dirty map.

if do_make_dirty_cube:

    print ""
    print "+-+-+-+-+-+-+-+-+-+-+-+-+"
    print "Making a dirty cube."
    print "+-+-+-+-+-+-+-+-+-+-+-+-+"
    print ""

    try:
        uvtaper
    except NameError:
        uvtaper = None

    if uvtaper == None:
        uv_taper_string = ''
    else:
        uv_taper_string = [str(uvtaper)+'arcsec',str(uvtaper)+'arcsec','0deg']

    niter = 0
    do_reset = True
    do_callclean = True
    do_savecopy = True

    bkup_ext = "dirty"
    logfile = cube_root+"_dirty.log"

    usemask = 'pb'
    mask = ''
    pbmask=pb_limit
    calcres = True
    calcpsf = True
    execfile('../scripts/callClean.py')

if do_revert_to_dirty:
    
    print ""
    print "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+"
    print "Resetting to the dirty cube or image."
    print "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+"
    print ""
        
    print ""
    print "Copying the previous dirty image to be the main cube."
    print ""
    
    os.system('rm -rf '+cube_root+'.image')
    os.system('rm -rf '+cube_root+'.model')
    os.system('rm -rf '+cube_root+'.mask')
    os.system('rm -rf '+cube_root+'.pb')
    os.system('rm -rf '+cube_root+'.psf')
    os.system('rm -rf '+cube_root+'.residual')
        
    bkup_ext = "dirty"
    os.system('cp -r '+cube_root+'_'+bkup_ext+'.image '+cube_root+'.image')
    os.system('cp -r '+cube_root+'_'+bkup_ext+'.model '+cube_root+'.model')
    os.system('cp -r '+cube_root+'_'+bkup_ext+'.mask '+cube_root+'.mask')
    os.system('cp -r '+cube_root+'_'+bkup_ext+'.pb '+cube_root+'.pb')
    os.system('cp -r '+cube_root+'_'+bkup_ext+'.psf '+cube_root+'.psf')
    os.system('cp -r '+cube_root+'_'+bkup_ext+'.residual '+cube_root+'.residual ')

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# GENERATE OR READ THE CLEAN MASK
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

# Generate the mask that will be used for cleaning. This is either a
# pure primary-beam based mask or an externally defined clean mask
# that we align to the astrometry of the input cube.

if do_read_in_clean_mask:

    print ""
    print "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+"
    print "Aligning the external clean mask to the dirty cube."
    print "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+"
    print ""

    mask_in = clean_mask_file
    mask_root = cube_root
    execfile('../scripts/alignMask.py')

    ia.open(cube_root+".pb")
    pb = ia.getchunk()
    ia.close()
    
    ia.open(cube_root+".mask")
    mask = ia.getchunk()
    mask *= (pb > pb_limit)
    ia.putchunk(mask)
    ia.close()

    mask = ''
    pb = 0.0
    
if do_use_pbmask and do_read_in_clean_mask == False:

    print ""
    print "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+"
    print "Using a primary beam based mask only."
    print "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+"
    print ""
 
    ia.open(cube_root+".pb")
    pb = ia.getchunk()
    ia.close()

    os.system('rm -rf '+cube_root+'.mask')
    os.system('cp -r '+cube_root+'.image '+cube_root+'.mask')
    
    ia.open(cube_root+".mask")
    mask = (pb > pb_limit)*1.0
    ia.putchunk(mask)
    ia.close()

    mask = ''
    pb = 0.0

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# RUN THE FIRST MULTISCALE CLEAN
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

# Run the main clean loop using the multiscale cleaning approach. This
# picks up from the dirty cube and executes a series of clean calls
# with the "multiscale" algorithm. Successive calls to clean allocate
# more and more iterations.

if do_multiscale_clean:

    print ""
    print "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+"
    print "Multiscale cleaning."
    print "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+"
    print ""

    print ""
    print "... Clean using MULTISCALE."
    print ""

    # Calculate the scales to use

    #try:
    #    outerscale
    #except NameError:
    #    outerscale = oversamp*cell_size*4.5
    #    print "Defaulting to an outer scale of "+str(outerscale)+" for multiscale."

    #if outerscale != None:
    #    outerscale_in_pix = ceil(outerscale / cell_size)
    #    scales_as_pix = [0]
    #    scales_as_angle = [0.0]
    #    for factor in range(10):
    #        scale = oversamp*2.0**(factor)
    #        if scale < outerscale_in_pix:
    #            scales_as_pix.append(int(round(scale)))
    #            scales_as_angle.append(scale*cell_size)
                
    scales_as_pix = []
    for scale in scales_as_angle:
        scales_as_pix.append(int(scale/cell_size))

    print "I will use the following scales: "
    print "... as pixels: ", str(scales_as_pix)
    print "... as arcseconds: ", str(scales_as_angle)

    # Calculate the threshold for cleaning

    try:
        multiscale_threshold
    except NameError:
        multiscale_threshold = None

    if multiscale_threshold == None:
        execfile('../scripts/statCleanCube.py')    
        this_threshold = snr_thresh* \
            imstat_residual['medabsdevmed'][0]/0.6745
        multiscale_threshold = str(this_threshold)+'Jy/beam'
    print "I will use threshold: ", multiscale_threshold
    threshold = multiscale_threshold

    # Note the number of channels

    vm = au.ValueMapping(input_vis)
    nchan = vm.spwInfo[0]['numChannels']

    # Figure out the number of iterations we will use    
    base_niter = 10*nchan
    base_cycle_niter = 100
    loop = 1

    # Initialize our output file and tracking of the flux in the model
    this_flux = 0.0
    try:
        multiscale_delta_thresh
    except NameError:
        multiscale_delta_thresh = 0.02
    delta_thresh = multiscale_delta_thresh
    proceed = True

    loop_record_file = cube_root+"_multiscale_record.txt"
    f = open(loop_record_file,'w')    
    f.write("# column 1: loop type\n")
    f.write("# column 2: loop number\n")
    f.write("# column 3: supplied threshold\n")
    f.write("# column 4: model flux at end of this clean\n")
    f.write("# column 5: fractional change in flux (current-previous)/current\n")
    f.write("# column 6: number of iterations allocated (not necessarily used)\n")
    f.close()

    while proceed == True and loop <= max_loop:
        
        # Steadily increase the iterations between statistical checks.        
        do_reset = False
        do_callclean = True
        
        if loop > 5:
            factor = 5
        else:
            factor = (loop-1)
        this_niter = base_niter*(2**factor)
        cycle_niter = base_cycle_niter*factor
        logfile = cube_root+"_loop_"+str(loop)+"_multiscale.log"        
            
        # Clean and associated tuning parameters
        calcres = False
        calcpsf = False
        minpsffraction = 0.5
        usemask = 'user'
        mask = ''
        deconvolver = 'multiscale'
        restoringbeam = 'common'
        niter = this_niter

        do_savecopy = False
        bkup_ext = "multiscaleloop"+str(loop)        

        print ""
        print "CALLING MULTISCALE CLEAN."
        print ""

        scales = scales_as_pix

        execfile('../scripts/callClean.py')

        # Run stats after the clean and write to the log file.
        
        execfile('../scripts/statCleanCube.py')    
        
        prev_flux = this_flux
        this_flux = imstat_model['sum'][0]
        delta_flux = (this_flux-prev_flux)/this_flux
        proceed = \
            (delta_flux > delta_thresh) and \
            (this_flux > 0.0)
        
        print ""
        print "******************************"
        print "MULTISCALE CLEAN LOOP "+str(loop)
        print "... threshold "+threshold
        print "... old flux "+str(prev_flux)
        print "... new flux "+str(this_flux)
        print "... fractional change "+str(delta_flux)+ \
            " compare to stopping criterion of "+str(delta_thresh)
        print "... proceed? "+str(proceed)
        print "******************************"
        print ""

        line = 'MULTISCALE '+str(loop)+' '+threshold+' '+str(this_flux)+ \
            ' '+str(delta_flux) + ' ' + str(this_niter)+ '\n' 
        f = open(loop_record_file,'a')
        f.write(line)
        f.close()

        if proceed == False:
            break
        loop += 1

    # Save a final copy
    do_savecopy = True
    do_callclean = False
    bkup_ext = "multiscale"
    execfile('../scripts/callClean.py')

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# REVERT TO THE MULTISCALE CASE IF DESIRED
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

if do_revert_to_multiscale:
    
    print ""
    print "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+"
    print "Resetting to the multiscale cube or image."
    print "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+"
    print ""
        
    print ""
    print "Copying the previous multiscale image to be the main cube."
    print ""
    
    os.system('rm -rf '+cube_root+'.image')
    os.system('rm -rf '+cube_root+'.model')
    os.system('rm -rf '+cube_root+'.mask')
    os.system('rm -rf '+cube_root+'.pb')
    os.system('rm -rf '+cube_root+'.psf')
    os.system('rm -rf '+cube_root+'.residual')
        
    bkup_ext = "multiscale"
    os.system('cp -r '+cube_root+'_'+bkup_ext+'.image '+cube_root+'.image')
    os.system('cp -r '+cube_root+'_'+bkup_ext+'.model '+cube_root+'.model')
    os.system('cp -r '+cube_root+'_'+bkup_ext+'.mask '+cube_root+'.mask')
    os.system('cp -r '+cube_root+'_'+bkup_ext+'.pb '+cube_root+'.pb')
    os.system('cp -r '+cube_root+'_'+bkup_ext+'.psf '+cube_root+'.psf')
    os.system('cp -r '+cube_root+'_'+bkup_ext+'.residual '+cube_root+'.residual ')
        
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# RUN THE SINGLE SCALE CLEAN
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

# Run a second clean loop using a single scale cleaning approach. This
# picks up from the multiscale cube and executes a series of clean
# calls with the "hogbom" algorithm. Successive calls to clean
# allocate more and more iterations.

if do_singlescale_clean:

    print ""
    print "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+"
    print "Single scale cleaning."
    print "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+"
    print ""

    print ""
    print "... Clean using one scale."
    print ""

    # Calculate the threshold for cleaning

    try:
        singlescale_threshold
    except NameError:
        singlescale_threshold = None

    if singlescale_threshold == None:
        execfile('../scripts/statCleanCube.py')    
        this_threshold = snr_thresh* \
            imstat_residual['medabsdevmed'][0]/0.6745
        singlescale_threshold = str(this_threshold)+'Jy/beam'
    print "I will use threshold: ", singlescale_threshold
    threshold = singlescale_threshold

    # Note the number of channels

    vm = au.ValueMapping(input_vis)
    nchan = vm.spwInfo[0]['numChannels']

    # Figure out the number of iterations we will use    
    base_niter = 10*nchan
    base_cycle_niter = 100
    loop = 1

    # Initialize our output file and tracking of the flux in the model
    this_flux = 0.0
    try:
        singlescale_delta_thresh
    except NameError:
        singlescale_delta_thresh = 0.02
    delta_thresh = singlescale_delta_thresh
    proceed = True

    loop_record_file = cube_root+"_singlescale_record.txt"
    f = open(loop_record_file,'w')    
    f.write("# column 1: loop type\n")
    f.write("# column 2: loop number\n")
    f.write("# column 3: supplied threshold\n")
    f.write("# column 4: model flux at end of this clean\n")
    f.write("# column 5: fractional change in flux (current-previous)/current\n")
    f.write("# column 6: number of iterations allocated (not necessarily used)\n")
    f.close()

    while proceed == True and loop <= max_loop:
        
        # Steadily increase the iterations between statistical checks.        
        do_reset = False
        do_callclean = True
        
        if loop > 5:
            factor = 5
        else:
            factor = (loop-1)
        this_niter = base_niter*(2**factor)
        cycle_niter = base_cycle_niter*factor
        logfile = cube_root+"_loop_"+str(loop)+"_singlescale.log"
            
        # Clean and associated tuning parameters
        calcres = False
        calcpsf = False
        minpsffraction = 0.5
        usemask = 'user'
        mask = ''
        deconvolver = 'hogbom'
        scales = [0]
        restoringbeam = 'common'
        niter = this_niter

        do_savecopy = False
        bkup_ext = "singlescaleloop"+str(loop)        

        print ""
        print "CALLING SINGLE SCALE CLEAN."
        print ""

        execfile('../scripts/callClean.py')

        # Run stats after the clean and write to the log file.
        
        execfile('../scripts/statCleanCube.py')    
        
        prev_flux = this_flux
        this_flux = imstat_model['sum'][0]
        delta_flux = (this_flux-prev_flux)/this_flux
        proceed = \
            (delta_flux > delta_thresh) and \
            (this_flux > 0.0)
        
        print ""
        print "******************************"
        print "SINGLESCALE CLEAN LOOP "+str(loop)
        print "... threshold "+threshold
        print "... old flux "+str(prev_flux)
        print "... new flux "+str(this_flux)
        print "... fractional change "+str(delta_flux)+ \
            " compare to stopping criterion of "+str(delta_thresh)
        print "... proceed? "+str(proceed)
        print "******************************"
        print ""

        line = 'SINGLESCALE '+str(loop)+' '+threshold+' '+str(this_flux)+ \
            ' '+str(delta_flux) + ' ' + str(this_niter)+ '\n' 
        f = open(loop_record_file,'a')
        f.write(line)
        f.close()

        if proceed == False:
            break
        loop += 1

    # Save a final copy
    do_savecopy = True
    do_callclean = False
    bkup_ext = "singlescale"
    execfile('../scripts/callClean.py')

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# POST PROCESS
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

if do_postprocess:

    print ""
    print "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+"
    print "Exporting the data."
    print "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+"
    print ""

    execfile('../scripts/exportToFITS.py')

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# TIMER
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

total_stop_time = time.time()

total_elapsed_time = (total_stop_time - total_start_time)/60.
print "This run took "+str(total_elapsed_time)+" minutes"
