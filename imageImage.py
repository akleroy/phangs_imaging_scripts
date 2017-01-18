tested_versions = ['4.6.0','4.7.0']
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
# CONTROL FLOW
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

# Here are the steps in easy form, for copy-and-paste for debugging.

#do_pickcellsize = True
#do_init = False
#do_make_dirty_mask = False
#do_revert_to_dirty = False
#do_start_with_pbmask = False
#do_clean_bright = False
#do_revert_to_bright = True
#do_make_model_mask = True
#do_clean_deep = True
#do_revert_to_deep = False
#do_postprocess = False

try:
    do_end_to_end
except NameError:
    do_end_to_end = False

if do_end_to_end:
    do_pickcellsize = True
    do_init = True
    do_make_dirty_mask = True
    do_revert_to_dirty = False
    do_clean_bright = True
    do_revert_to_bright = False
    do_make_model_mask = True
    do_clean_deep = True
    do_revert_to_deep = False
    do_postprocess = True

try:
    do_pickcellsize
except NameError:
    do_pickcellsize = True

try:
    do_init
except NameError:
    do_init = True

try:
    do_revert_to_dirty
except NameError:
    do_revert_to_dirty = False

try:
    do_make_dirty_mask
except NameError:
    do_make_dirty_mask = True

try:
    do_clean_bright
except NameError:
    do_clean_bright = True

try:
    do_start_with_pbmask
except NameError:
    do_start_with_pbmask = True

try:
    do_revert_to_bright
except NameError:
    do_revert_to_bright = False

try:
    do_make_model_mask
except NameError:
    do_make_model_mask = False

try:
    do_clean_deep
except NameError:
    do_clean_deep = True

try:
    do_revert_to_deep
except NameError:
    do_revert_to_deep = False

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

    # This is not currently working ...

    #print "... central field = ", au_centralField
    #phase_center = au_centralField

if do_init:

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
    pbmask=0.2
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

if do_make_dirty_mask:

    print ""
    print "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+"
    print "Making a mask based on the dirty cube or image."
    print "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+"
    print ""

    recipe = 'clipandsmooth'
    execfile('../scripts/makeMask.py')

    os.system('rm -rf '+cube_root+'_widearea.mask')
    os.system('cp -r '+cube_root+'.mask '+cube_root+'_widearea.mask')

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# CLEAN DOWN TO A BRIGHT THRESHOLD
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

if do_clean_bright and (do_revert_to_bright == False):

    print ""
    print "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+"
    print "Resuming cleaning, targeting bright emission."
    print "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+"
    print ""

    print ""
    print "... first, I will clean using a wide mask and a single scale"
    print "... down to a 3sigma threshold. I will stop when the flux in"
    print "... the model converges for successive cleans at the ~2% level."
    print "... This data product is the 'bright' cube."
    print ""

    vm = au.ValueMapping(input_vis)
    nchan = vm.spwInfo[0]['numChannels']
    
    base_niter = 10*nchan
    base_cycle_niter = 200
    loop = 1
    max_loop = 20
    deconvolver = "hogbom"    

    this_flux = 0.0
    delta_thresh = 0.02
    proceed = True

    loop_record_file = cube_root+"_bright_record.txt"
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
        
        if loop > 6:
            factor = 6
        else:
            factor = (loop-1)
        niter = base_niter*(2**factor)
        cycle_niter = base_cycle_niter*factor
        logfile = cube_root+"_loop_"+str(loop)+"_bright.log"
        
        # Figure out a current threshold in a very crude way.
        execfile('../scripts/statCleanCube.py')    
        this_threshold = bright_snr_thresh* \
            imstat_residual['medabsdevmed'][0]/0.6745
        threshold = str(this_threshold)+'Jy/beam'
            
        # Clean.

        calcres = False
        minpsffraction = 0.5

        if do_start_with_pbmask == False:
            usemask = 'user'
            mask = ''
            mask_file = cube_root+'.mask'
        else:
            usemask = 'pb'
            mask = ''
            pbmask = 0.2

        # Keep a running backup of the previous iteration.

        do_savecopy = True
        bkup_ext = "prev"

        execfile('../scripts/callClean.py')
        
        # Run stats after the clean.
        
        execfile('../scripts/statCleanCube.py')    
        
        prev_flux = this_flux
        this_flux = imstat_model['sum'][0]
        delta_flux = (this_flux-prev_flux)/this_flux
        proceed = \
            (delta_flux > delta_thresh) and \
            (this_flux > 0.0)
        
        print ""
        print "***************"
        print "BRIGHT CLEAN LOOP "+str(loop)
        print "... threshold "+threshold
        print "... old flux "+str(prev_flux)
        print "... new flux "+str(this_flux)
        print "... fractional change "+str(delta_flux)+ \
            " compare to stopping criterion of "+str(delta_thresh)
        print "... proceed? "+str(proceed)
        print "***************"        
        print ""

        line = 'BRIGHT '+str(loop)+' '+threshold+' '+str(this_flux)+ \
            ' '+str(delta_flux) + ' ' + str(niter)+ '\n' 
        f = open(loop_record_file,'a')
        f.write(line)
        f.close()

        if proceed == False:
            break
        loop += 1
    
    # We are done. Now make a copy of the cube cleaned down to S/N 3
    # using only point sources. The next step can be used to revert to
    # this bright emission only.
    
    bkup_ext = "bright"
    
    do_reset = False
    do_callclean = False
    do_savecopy = True
    
    execfile('../scripts/callClean.py')
        
if do_revert_to_bright:
    
    print ""
    print "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+"
    print "Resetting to the bright signal cube or image."
    print "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+"
    print ""
    
    os.system('rm -rf '+cube_root+'.image')
    os.system('rm -rf '+cube_root+'.model')
    os.system('rm -rf '+cube_root+'.mask')
    os.system('rm -rf '+cube_root+'.pb')
    os.system('rm -rf '+cube_root+'.psf')
    os.system('rm -rf '+cube_root+'.residual')
        
    bkup_ext = "bright"
    os.system('cp -r '+cube_root+'_'+bkup_ext+'.image '+cube_root+'.image')
    os.system('cp -r '+cube_root+'_'+bkup_ext+'.model '+cube_root+'.model')
    os.system('cp -r '+cube_root+'_'+bkup_ext+'.mask '+cube_root+'.mask')
    os.system('cp -r '+cube_root+'_'+bkup_ext+'.pb '+cube_root+'.pb')
    os.system('cp -r '+cube_root+'_'+bkup_ext+'.psf '+cube_root+'.psf')
    os.system('cp -r '+cube_root+'_'+bkup_ext+'.residual '+cube_root+'.residual ')

if do_make_model_mask:

    print ""
    print "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+"
    print "Making a mask based on the model."
    print "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+"
    print ""

    recipe = 'modelmask'
    execfile('../scripts/makeMask.py')

    os.system('rm -rf '+cube_root+'_modelbased.mask')
    os.system('cp -r '+cube_root+'.mask '+cube_root+'_modelbased.mask')

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# CLEAN DEEPLY AROUND THE CURRENT MODEL
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

if do_clean_deep and (do_revert_to_deep == False):

    print ""
    print "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+"
    print "Resuming cleaning, pushing in to the noise."
    print "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+"
    print ""

    # First figure out how much flux, approximately, lies inside the
    # mask and has not yet been cleaned. This tells us about what we
    # aim for in this stage (though only within about a factor of
    # two).

    print ""
    print "..............................................."
    print "Estimating the flux that I expect to clean."
    print "..............................................."
    print ""

    stat_image = imstat(cube_root+'.image'
                        , mask=cube_root+'.mask')
    stat_resid = imstat(cube_root+'.residual'
                        , mask=cube_root+'.mask')
    stat_model = imstat(cube_root+'.model')

    frac_to_clean = stat_resid['sum'][0]/stat_image['sum'][0]
    frac_cleaned = 1.0 - frac_to_clean
    this_flux = stat_model['sum'][0]
    flux_at_start = this_flux
    flux_to_clean = this_flux / frac_cleaned

    vm = au.ValueMapping(input_vis)
    nchan = vm.spwInfo[0]['numChannels']

    print "Based on the current residuals I expect to clean about "+ \
        str(flux_to_clean)+" Jy"
    
    print ""
    print "......................................................"
    print "Now I will clean in to deeper levels in a narrow mask."
    print "......................................................"
    print ""    

    proceed = True

    delta_thresh = 0.01
    base_niter = 200*nchan
    base_cycle_niter = 100
    loop = 1
    max_loop = 20
    deconvolver = "hogbom"    

    loop_record_file = cube_root+"_deep_record.txt"
    f = open(loop_record_file,'w')    
    f.write("# column 1: loop type\n")
    f.write("# column 2: loop number\n")
    f.write("# column 3: supplied threshold\n")
    f.write("# column 4: model flux at end of this clean\n")
    f.write("# column 5: fractional change in flux (current-previous)/current\n")
    f.write("# column 6: number of iterations allocated (not necessarily used)\n")
    f.close()
    
    first = True

    while proceed == True and loop <= max_loop:

            # Steadily increase the iterations between statistical checks.

            do_reset = False
            do_callclean = True
            do_savecopy = False

            if loop > 8:
                factor = 8
            else:
                factor = loop
            niter = base_niter
            cycle_niter = base_cycle_niter
            logfile = cube_root+"_loop_"+str(loop)+"_deepclean.log"
            
            # Figure out a current threshold in a very crude way.
            execfile('../scripts/statCleanCube.py')    
            threshold = '0Jy/beam'
            
            # Clean.
            
            calcres = False            
            minpsffraction = 0.8

            usemask = 'user'
            mask = ''
            mask_file = cube_root+'.mask'

            # Keep a running backup of the previous iteration.
            
            do_savecopy = True
            bkup_ext = "prev"

            execfile('../scripts/callClean.py')

            # Run stats after the clean.
            
            execfile('../scripts/statCleanCube.py')    
            
            prev_flux = this_flux
            this_flux = imstat_model['sum'][0]
            delta_flux = (this_flux-prev_flux)/this_flux

            proceed = True

            if this_flux > flux_to_clean:
                proceed = False

            if this_flux <= 0.0:
                proceed = False

            if first:
                first = False
                if delta_flux < delta_thresh:
                    proceed = False
            
            print ""
            print "***************"
            print "DEEP CLEAN LOOP "+str(loop)
            print "... threshold "+threshold
            print "... old flux "+str(prev_flux)
            print "... new flux "+str(this_flux)
            print "... target flux "+str(flux_to_clean)
            print "... fractional change "+str(delta_flux)+ \
                " compare to stopping criterion of "+str(delta_thresh)
            print "... proceed? "+str(proceed)
            print "***************"        
            print ""

            line = 'DEEP '+str(loop)+' '+threshold+' '+str(this_flux)+ \
                ' '+str(delta_flux) + ' ' + str(niter) + '\n'
            f = open(loop_record_file,'a')    
            f.write(line)
            f.close()

            if proceed == False:
                break
            loop += 1
    
    # Make a copy of the fully cleaned cube.
    bkup_ext = "deep"
        
    do_reset = False
    do_callclean = False
    do_savecopy = True
    
    execfile('../scripts/callClean.py')

if do_revert_to_deep:
    
    print ""
    print "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+"
    print "Resetting to the deep cube or image."
    print "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+"
    print ""
    
    os.system('rm -rf '+cube_root+'.image')
    os.system('rm -rf '+cube_root+'.model')
    os.system('rm -rf '+cube_root+'.mask')
    os.system('rm -rf '+cube_root+'.pb')
    os.system('rm -rf '+cube_root+'.psf')
    os.system('rm -rf '+cube_root+'.residual')
        
    bkup_ext = "deep"
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
    print "Post-processing and exporting the data."
    print "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+"
    print ""

    do_calc_beam = True
    do_process = True
    do_shrink = False
    do_fits = True

    execfile('../scripts/postProcessCubes.py')

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# TIMER
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

total_stop_time = time.time()

total_elapsed_time = (total_stop_time - total_start_time)/60.
print "This run took "+str(total_elapsed_time)+" minutes"
