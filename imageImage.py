tested_versions = ['4.6.0','4.7.0']
this_version = (casa['build']['version']).split('-')[0]
if this_version not in tested_versions:
    print "The script hasn't been verified for this version of CASA."
    print "This version of CASA is "+this_version
    print "Tested versions are "+str(tested_versions)

# Hate this ... but directories
execfile('../scripts/auExtensions.py')

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# TIMER
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

import time
total_start_time = time.time()

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# CONTROL FLOW
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

try:
    do_pickcellsize
except NameError:
    do_pickcellsize = True

try:
    do_init
except NameError:
    do_init = True

try:
    do_mask
except NameError:    
    do_mask = True

try:
    resume_from_singlescale
except NameError:
    resume_from_singlescale = False

try:
    do_clean
except NameError:
    do_clean = True

try:
    do_postprocess
except NameError:
    do_postprocess = False

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# Inputs
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

abort = False

print ""
print "---------- imageImage.py ----------"
print ""

bright_snr_thresh = 3.0
thresh_step = 0.5

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
    do_init = False
    do_mask = False
    do_clean = False
    do_postprocess = False
 
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# MAKE THE DIRTY CUBE
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

if do_pickcellsize:

    # Factor by which to oversample the beam
    oversamp = 5

    # Cell size implied by baseline distribution
    au_cellsize, au_imsize, au_centralField = \
        au.pickCellSize(input_vis, imsize=True, npix=oversamp)
    xextent = au_cellsize*au_imsize[0]
    yextent = au_cellsize*au_imsize[1]

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

    # This is not currently working ...

    #print "... central field = ", au_centralField
    #phase_center = au_centralField

if do_init:

    print ""
    print "+-+-+-+-+-+-+-+-+-+-+-+-+"
    print "Making a dirty cube."
    print "+-+-+-+-+-+-+-+-+-+-+-+-+"
    print ""

    # Work out the uv taper
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

    execfile('../scripts/callClean.py')

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# MAKE A MASK
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

if do_makemask:

    # smooth to large (~10-15") scales and make mask. Expand it even a
    # bit more. Also look for very bright things. 

    # Alternatively, regrid an external mask on to the parent cube.

    pass

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# CLEAN
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

if do_cleancube:

    print ""
    print "+-+-+-+-+-+-+-+-+-+"
    print "Resuming cleaning."
    print "+-+-+-+-+-+-+-+-+-+"
    print ""

    if resume_from_singlescale == False:

        print ""
        print "... first, I will clean using a single scale down to a 3sigma threshold"
        print "... and minimal masking. I will stop when the flux in the model converges"
        print "... at the 2% level."
        print ""

        nchan = 1
        base_niter_per_channel = 100
        base_cycle_niter = 100
        loop = 1
        max_loop = 10
        deconvolver = "hogbom"    

        this_flux = 0.0
        delta_thresh = 0.02
        proceed = True

        while proceed == True and loop <= max_loop:

            # Steadily increase the iterations between statistical checks.

            do_reset = False
            do_callclean = True
            do_savecopy = False

            niter = base_niter_per_channel*(2**loop)*nchan
            cycle_niter = base_cycle_niter*(2**loop)/2
            logfile = cube_root+"_loop_"+str(loop)+"_singlescale.log"
            
            # Figure out a current threshold in a very crude way.
            usemask = 'pb'
            mask = ''
            execfile('../scripts/statCleanCube.py')    
            threshold = str(bright_snr_thresh*imstat_residual['medabsdevmed'][0]/0.6745)+'Jy/beam'
            
            # Clean.
            
            calcres = False
            execfile('../scripts/callClean.py')

            # Run stats after the clean.
            
            execfile('../scripts/statCleanCube.py')    
            
            prev_flux = this_flux
            this_flux = imstat_model['sum'][0]
            delta_flux = (this_flux-prev_flux)/this_flux
            proceed = abs(delta_flux) > delta_thresh
            
            print ""
            print "***************"
            print "SINGLE SCALE LOOP "+str(loop)
            print "... threshold "+threshold
            print "... old flux "+str(prev_flux)
            print "... new flux "+str(this_flux)
            print "... fractional change "+str(delta_flux)+ \
                " compare to stopping criterion of "+str(delta_thresh)
            print "... proceed? "+str(proceed)
            print "***************"        
            print ""

            if proceed == False:
                break
            loop += 1

        # Make a copy of the cube cleaned down to S/N 3 using only point sources.
        
        bkup_ext = "singlescale"
        
        do_reset = False
        do_callclean = False
        do_savecopy = True

        execfile('../scripts/callClean.py')
        
    else:
        
        # If the "resume_from_singlescale" option is set, copy the
        # backed up singlescale image to be the new cube. This flow
        # seems likely to be deprecated once we lock the algorithm.
        
        print ""
        print "Copying the previous single scale image to be the main cube."
        print ""
        
        os.system('rm -rf '+cube_root+'.image')
        os.system('rm -rf '+cube_root+'.model')
        os.system('rm -rf '+cube_root+'.mask')
        os.system('rm -rf '+cube_root+'.pb')
        os.system('rm -rf '+cube_root+'.psf')
        os.system('rm -rf '+cube_root+'.residual')
        
        bkup_ext = "singlescale"
        os.system('cp -r '+cube_root+'_'+bkup_ext+'.image '+cube_root+'.image')
        os.system('cp -r '+cube_root+'_'+bkup_ext+'.model '+cube_root+'.model')
        os.system('cp -r '+cube_root+'_'+bkup_ext+'.mask '+cube_root+'.mask')
        os.system('cp -r '+cube_root+'_'+bkup_ext+'.pb '+cube_root+'.pb')
        os.system('cp -r '+cube_root+'_'+bkup_ext+'.psf '+cube_root+'.psf')
        os.system('cp -r '+cube_root+'_'+bkup_ext+'.residual '+cube_root+'.residual ')

    print ""
    print "... now I will clean around the peaks pushing into the noise."
    print ""

    nchan = 1
    base_niter_per_channel = 10000
    base_cycle_niter = 10000
    loop = 1
    max_loop = 10
    deconvolver = "hogbom"    

    print ""
    print "... I am making a mask around the model I have built so far."
    print ""
    
    execfile('../scripts/statCleanCube.py')                

    current_thresh = bright_snr_thresh
    current_rms = imstat_residual['medabsdevmed'][0]/0.6745

    target_beam = str((imhead(cube_root+'.image'))['restoringbeam']['major']['value']*2.0)+'arcsec'

    os.system('rm -rf '+cube_root+'_smoothed_model.image')    
    imsmooth(imagename=cube_root+'.model',
             targetres=True, major=target_beam, minor=target_beam, pa='0deg',
             outfile=cube_root+'_smoothed_model.image', overwrite=True)

    os.system('rm -rf '+cube_root+'_mask.image')
    thresh_for_mask = 0.5*current_rms
    immath(imagename = cube_root+'_smoothed_model.image',
           outfile = cube_root+'_mask.image',
           expr = 'iif(IM0 > '+str(thresh_for_mask) +',1.0,0.0)')

    ia.open(cube_root+'_mask.image')
    mask = ia.getchunk()
    ia.done()

    ia.open(cube_root+'.mask')
    ia.putchunk(mask*1.0)
    ia.done()
    
    print ""
    print "... Proceeding with the clean."
    print ""

    this_flux = imstat_model['sum'][0]
    delta_thresh = 0.02
    proceed = True
    
    while proceed == True and loop <= max_loop:

            # Steadily increase the iterations between statistical checks.

            do_reset = False
            do_callclean = True
            do_savecopy = False

            niter = base_niter_per_channel*(2**loop)*nchan
            cycle_niter = base_cycle_niter*(2**loop)/2
            logfile = cube_root+"_loop_"+str(loop)+"_deepclean.log"
            
            # Figure out a current threshold in a very crude way.
            usemask = 'user'
            mask = ''            
            execfile('../scripts/statCleanCube.py')    
            threshold = '0Jy/beam'
            
            # Clean.
            
            calcres = False            
            mask_file = cube_root+'.mask'
            execfile('../scripts/callClean.py')

            # Run stats after the clean.
            
            execfile('../scripts/statCleanCube.py')    
            
            prev_flux = this_flux
            this_flux = imstat_model['sum'][0]
            delta_flux = (this_flux-prev_flux)/this_flux
            proceed = abs(delta_flux) > delta_thresh
            
            print ""
            print "***************"
            print "DEEP CLEAN SCALE LOOP "+str(loop)
            print "... threshold "+threshold
            print "... old flux "+str(prev_flux)
            print "... new flux "+str(this_flux)
            print "... fractional change "+str(delta_flux)+ \
                " compare to stopping criterion of "+str(delta_thresh)
            print "... proceed? "+str(proceed)
            print "***************"        
            print ""

            if proceed == False:
                break
            loop += 1
    
    # Make a copy of the fully cleaned cube.
    bkup_ext = "deepclean"
        
    do_reset = False
    do_callclean = False
    do_savecopy = True
    
    execfile('../scripts/callClean.py')

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# POST PROCESS
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

if do_postprocess:
    do_process = True
    do_shrink = True
    do_fits = True
    execfile('.../scripts/postProcessCubes.py')

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# TIMER
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

total_stop_time = time.time()

total_elapsed_time = (total_stop_time - total_start_time)/60.
print "This run took "+str(total_elapsed_time)+" minutes"
