# Note that "IMPROVEMENT TBD" items are scattered throughout the
# program. Search on this to find areas where some clean up or a next
# logical step could be worked out.

tested_versions = ['4.6.0','4.7.0','4.7.1']
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
    do_clean = True
    do_singlescale = True
    do_revert_to_multiscale = False
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
    do_clean
except NameError:
    do_clean = True

try:
    do_revert_to_multiscale
except NameError:
    revert_to_multiscale = False

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
    scales_to_use
except NameError:
    scales_to_use = [0,2,4,8,16,32,64]

try:
    bright_snr_thresh 
except NameError:
    bright_snr_thresh = 4.0

try:
    pb_limit
except NameError:
    pb_limit = 0.25

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
    do_clean = False
    do_revert_to_multiscale = False
    do_singlescale = False
    do_postprocess = False
 
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# PICK THE INPUT PARAMETERS FOR CLEAN
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

# Here, determine the choices for the CLEAN image - the cell size,
# image size, etc. Use the analysis utilities and the input UV data
# set to do this.

# IMPROVEMENT TBD: Pick up these numbers from an existing image. Or
# perhaps CLEAN doesn't need them if it resumes from an existing image
# on disk? In any case, figure out a way to avoid the case where this
# step misfires and down some processes down the line come tumbling
# down.

# IMPROVEMENT TBD: Right now we specify the field center by hand. This
# is actually probably okay. But doing this automatically might be
# useful down the line. This is not high priority.

if do_pickcellsize:
    # Factor by which to oversample the beam
    oversamp = 5

    # Cell size implied by baseline distribution
    au_cellsize, au_imsize, au_centralField = \
        au.pickCellSize(input_vis, imsize=True, npix=oversamp)
    xextent = au_cellsize*au_imsize[0]*1.2
    yextent = au_cellsize*au_imsize[1]*1.2

    # If a u-v taper is present, then take that as an alternate beam
    # size.
    try:
        uvtaper
    except NameError:
        uvtaper = None

    if uvtaper != None:
        print "I account for a uv taper of "+str(uvtaper)+" arcsec"
        cell_size_taper = uvtaper / oversamp
        au_cellsize = sqrt(au_cellsize**2 + cell_size_taper**2)

    # Make the cell size a nice round number
    if au_cellsize < 0.1:
        cell_size = au_cellsize
    if au_cellsize >= 0.1 and au_cellsize < 0.5:
        cell_size = floor(au_cellsize/0.05)*0.05
    if au_cellsize >= 0.5 and au_cellsize < 1.0:
        cell_size = floor(au_cellsize/0.1)*0.1
    if au_cellsize >= 1.0 and au_cellsize < 2.0:
        cell_size = floor(au_cellsize/0.25)*0.25
    if au_cellsize >= 2.0 and au_cellsize < 5.0:
        cell_size = floor(au_cellsize/0.5)*0.5
    if au_cellsize >= 5.0:
        cell_size = floor(au_cellsize/1.0)*0.5

    # If we taper, then the cell size may be forced
    print "I calculated cell size: ", cell_size
    try:
        force_cell_size
    except NameError:
        force_cell_size = None

    if force_cell_size != None:
        cell_size = force_cell_size
        print "I will force cell size: ", cell_size

    # Now make the image size a nice round number
    need_cells_x = xextent / cell_size
    need_cells_y = yextent / cell_size

    # Valid image sizes are even and multiples of 3, 5, 7
    valid_sizes = []
    for ii in range(10):
        for kk in range(3):
            for jj in range(3):
                valid_sizes.append(2**(ii+1)*5**(jj)*3**(kk))
    valid_sizes.sort()
    valid_sizes = np.array(valid_sizes)

    cells_x = np.min(valid_sizes[valid_sizes > need_cells_x])
    cells_y = np.min(valid_sizes[valid_sizes > need_cells_y])

    image_size = [int(cells_x), int(cells_y)]
    cell_size_string = str(cell_size)+'arcsec'

    print "I propose the following:"
    print "... cell size = "+cell_size_string
    print "... image size = ", image_size

    # This is not currently working ... we specify the field center by hand.

    #print "... central field = ", au_centralField
    #phase_center = au_centralField

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
# RUN THE MULTISCALE CLEAN
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

# Run the main clean loop using the multiscale cleaning approach. This
# picks up from the dirty cube and executes a series of clean calls
# with the "multiscale" algorithm. Successive calls to clean allocate
# more and more iterations. Between executions, curate the model
# slightly. We don't expect negative components, so delete negative
# values from the model. 

if do_clean:

    print ""
    print "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+"
    print "Multiscale cleaning."
    print "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+"
    print ""

    print ""
    print "... Cleaning using MULTISCALE. Looking for 2% convergence between iterations."
    print ""

    vm = au.ValueMapping(input_vis)
    nchan = vm.spwInfo[0]['numChannels']

    # Figure out the number of iterations we will use
    
    base_niter = 10*nchan
    base_cycle_niter = 100
    loop = 1

    # Initialize our output file and tracking of the flux in the model

    this_flux = 0.0
    delta_thresh = 0.02
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
        
        # Figure out a current threshold

        execfile('../scripts/statCleanCube.py')    
        this_threshold = bright_snr_thresh* \
            imstat_residual['medabsdevmed'][0]/0.6745
        threshold = str(this_threshold)+'Jy/beam'
            
        # Clean.

        calcres = False
        minpsffraction = 0.5
        usemask = 'user'
        mask = ''
        do_savecopy = False
        deconvolver = 'multiscale'

        scales = scales_to_use
        niter = this_niter

        print ""
        print "CALLING MULTISCALE CLEAN."
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

    # Save a final copy of the multiscale
        
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
