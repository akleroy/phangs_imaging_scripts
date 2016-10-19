# Script to image one of our 1mm galaxies. Creates native resolution
# and tapered imaging of 12CO 2-1, 12C18O 2-1, mm continuum, and
# "channel 0" maps of 12CO 2-1 and C18O 2-1.

tested_versions = ['4.6.0']
this_version = casa['build']['version']
if this_version not in tested_versions:
    print "The script hasn't been verified for this version of CASA."
    print "This version of CASA is "+this_version
    print "Tested versions are "+str(tested_versions)
else:
    print "The script has been verified for this version of CASA."

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# TIMER
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

import time
start_time = time.time()

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# TEST INPUTS
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

# --------------------------------------
# Control flow
# --------------------------------------

try:
    do_copy
except NameError:    
    do_copy = True

try:
    do_mask
except NameError:    
    do_mask = True

try:
    do_co21
except NameError:    
    do_co21 = True

try:
    do_cont
except NameError:    
    do_cont = True

try:
    do_c18o
except NameError:    
    do_c18o = True

try:
    has_c18o
except NameError:
    has_c18o = True

try:
    do_cleanup
except NameError:    
    do_cleanup = True

# --------------------------------------
# Inputs
# --------------------------------------

try:
    skip_mask
except NameError:
    skip_mask = False

try:
    calibrated_file
except NameError:
    print "Please define a calibrated data file (calibrated_file string)."

try:
    split_mosaic
except NameError:
    print "Please indicate if the file is a split mosaic (split_mosaic=False)."

try:
    co21_spw
except NameError:
    print "Defaulting CO21 SPW to 2 (co21_spw='2')."
    co21_spw = '2'

try:
    co21_spw_7m
except NameError:
    print "Defaulting CO21 SPW for the 7m array to 0 (co21_spw_7m='0')."
    co21_spw_7m = '0'

try:
    c18o_spw
except NameError:
    print "Defaulting C18O SPW to 3 (c18o_spw='3')."
    c18o_spw = '3'

try:
    c18o_spw_7m
except NameError:
    print "Defaulting C18O21 SPW for the 7m array to 0 (c18o21_spw_7m='2')."
    c18o21_spw_7m = '2'

try:
    gal
except NameError:
    print "Please define a galaxy name (gal string)."

try:
    start_vel
except NameError:
    print "Please specify a starting velocity (start_vel string)."

try:
    nchan
except NameError:
    print "Please specify a number of channels for CO 2-1 (nchan integer)."

try:
    c18o_deltav
except NameError:
    print "Defaulting C18O channel width to 5km/s (c18o_deltav='5km/s')."
    c18o_deltav = '5km/s'

try:
    flag_co21
except NameError:
    print "Please specify flagging to remove 12CO from the continuum (flag_co21 string)."

try:
    flag_c18o21
except NameError:
    print "Please specify flagging to remove C18O from the continuum (flag_c18o21 string)."

try:
    phase_center
except NameError:
    print "Please specify a phase center for the map (phase_center string)."

try:
    imsize
except NameError:
    print "Please specify an image size (imsize [integer, integer]))."

try:
    imsize_lowres
except NameError:
    print "Please specify an image size low resolution images (imsize_lowres [integer, integer]))."

try:
    thresh_factor
except NameError:
    print "Please specify threshold factor (thresh_factor string)."

try:
    target_beam_co21
except NameError:
    print "Please specify a target round beam for CO21 (target_beam_co21 string)."

try:
    target_beam_c18o
except NameError:
    print "Please specify a target round beam for C18O21 (target_beam_c18o string)."

try:
    target_beam_cont
except NameError:
    print "Please specify a target round beam for continuum (target_beam_cont string)."
 
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# COPY DATA FROM ITS ORIGINAL LOCATION AND CARVE IT UP
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

if do_copy:

    print "--------------------------------------"
    print "(1) Prepping data for imaging."
    print "--------------------------------------"

    # ------------------------------------------------
    # Copy from the pipeline directory
    # ------------------------------------------------

    print "... Copying from original calibrated data."

    if split_mosaic:
        print "... ... treating this as a split mosaic."

        counter = 1
        to_concat = []

        for this_file in calibrated_file:
            infile = calibrated_file[counter-1]
            outfile = gal+'_'+str(counter)+'.ms'
            os.system('rm -rf '+outfile)
            os.system('rm -rf '+outfile+'.flagversions')
            os.system('scp -r '+infile+' '+outfile)
            
            counter += 1
    else:
        print "... ... treating this as a single mosaic."
        os.system('rm -rf '+gal+'.ms')
        os.system('rm -rf '+gal+'.ms.flagversions')
        os.system('scp -r '+calibrated_file+' '+gal+'.ms')

    if has_7m:
        print "... ... copying the 7m data."

        if split_mosaic:
            print "... ... treating this as a split mosaic."

            counter = 1
            to_concat = []

            for this_file in calibrated_file:
                infile = calibrated_file[counter-1]
                outfile = gal+'_'+str(counter)+'_7m.ms'
                os.system('rm -rf '+outfile)
                os.system('rm -rf '+outfile+'.flagversions')
                os.system('scp -r '+infile+' '+outfile)
                
                counter += 1
        else:
            print "... ... treating this as a single mosaic."
            os.system('rm -rf '+gal+'_7m.ms')
            os.system('rm -rf '+gal+'_7m.ms.flagversions')
            os.system('scp -r '+calibrated_7m_file+' '+gal+'_7m.ms')
        
    # ------------------------------------------------
    # CO 2-1
    # ------------------------------------------------

    print "... Processing the CO 2-1 into an MS appropriate for imaging."

    # Split out the CO 2-1, averaging in time and frequency to make
    # the data volume smaller.

    if split_mosaic:

        counter = 1
        to_concat = []

        for this_file in calibrated_file:
                
            this_infile =  gal+'_'+str(counter)+'.ms'
            this_outfile = gal+'_'+str(counter)+'_co21_bin.ms'
            
            os.system('rm -rf '+this_outfile)
            split(vis=this_infile
                  , field=field
                  , datacolumn='DATA'
                  , spw=co21_spw, timebin='120s'
                  , outputvis=this_outfile)
            
            to_concat.append(this_outfile)

            counter += 1
        
        if has_7m:
            
            counter = 1
            for this_file in calibrated_7m_file:

                this_infile =  gal+'_'+str(counter)+'_7m.ms'
                this_outfile = gal+'_'+str(counter)+'_7m_co21_bin.ms'
            
                os.system('rm -rf '+this_outfile)
                split(vis=this_infile
                      , field=field_7m
                      , datacolumn='DATA'
                      , spw=co21_spw_7m, timebin='120s'
                      , outputvis=this_outfile)
            
                to_concat.append(this_outfile)

                counter += 1

        os.system('rm -rf '+gal+'_co21_bin.ms')
        concat(vis=to_concat,
               concatvis=gal+'_co21_bin.ms')

    else:

        print "... ... treating this as a single mosaic."

        if has_7m:
            outputvis = gal+'_12m_co21_bin.ms'
        else:
            outputvis = gal+'_co21_bin.ms'

        os.system('rm -rf '+outputvis)        
        split(vis=gal+'.ms' 
              , datacolumn='DATA'
              , spw=co21_spw, width=5, timebin='120s'
              , outputvis=outputivs
              , field=field)

        if has_7m:
            outputvis = gal+'_7m_co21_bin.ms'
            os.system('rm -rf '+outputvis)        
            split(vis=gal+'_7m.ms'
                  , datacolumn='DATA'
                  , spw=co21_spw_7m, width=5, timebin='120s'
                  , outputvis=outputivs
                  , field=field_7m)
            
            os.system('rm -rf '+gal+'_co21_bin.ms')
            concat(vis=[gal+'_12m_co21_bin.ms',
                        gal+'_7m_co21_bin.ms']
                   concatvis=gal+'_co21_bin.ms')

    # Regrid CO 2-1 to a regular, 2.5km/s velocity grid

    os.system('rm -rf '+gal+'_co21.ms')
    mstransform(vis=gal+'_co21_bin.ms',
                datacolumn='DATA',
                outputvis=gal+'_co21.ms',
                combinespws=True,
                regridms=True,
                mode='velocity',
                width='2.5km/s',
                start=start_vel,
                nchan=nchan,
                restfreq='230.53800GHz' ,
                outframe='lsrk',
                veltype='radio',            
            )

    # ------------------------------------------------
    # C18O 2-1
    # ------------------------------------------------

    if has_c18o:

        # Split out the C18O 2-1

        print "... Processing the C18O 2-1 into an MS appropriate for imaging."

        if split_mosaic:
            print "... ... treating this as a split mosaic."

            counter = 1
            to_concat = []

            for this_file in calibrated_file:
                
                this_infile =  gal+'_'+str(counter)+'.ms'
                this_outfile = gal+'_'+str(counter)+'_c18o_bin.ms'

                os.system('rm -rf '+this_outfile)
                split(vis=this_infile
                      , field=field
                      , datacolumn='DATA'
                      , spw=c18o_spw, timebin='120s'
                      , outputvis=this_outfile)

                to_concat.append(this_outfile)

                counter += 1

            if has_7m:
                
                counter = 1
                for this_file in calibrated_7m_file:

                    this_infile =  gal+'_'+str(counter)+'_7m.ms'
                    this_outfile = gal+'_'+str(counter)+'_7m_c18o_bin.ms'
            
                    os.system('rm -rf '+this_outfile)
                    split(vis=this_infile
                          , field=field_7m
                          , datacolumn='DATA'
                          , spw=c18o_spw_7m, timebin='120s'
                          , outputvis=this_outfile)
            
                    to_concat.append(this_outfile)
                    counter += 1
        
            os.system('rm -rf '+gal+'_c18o_bin.ms')
            concat(vis=to_concat,
                   concatvis=gal+'_c18o_bin.ms')
        else:
            print "... ... treating this as a single mosaic."

            if has_7m:
                outputvis = gal+'_12m_c18o_bin.ms'
            else:
                outputvis = gal+'_c18o_bin.ms'

            os.system('rm -rf '+outputvis)
            split(vis=gal+'.ms' 
                  , field=field
                  , datacolumn='DATA'
                  , spw=c18o_spw, timebin='120s'
                  , outputvis=gal+'_c18o_bin.ms')
            
            if has_7m:
                outputvis = gal+'_7m_c18o_bin.ms'
                os.system('rm -rf '+outputvis)        
                split(vis=gal+'_7m.ms'
                      , datacolumn='DATA'
                      , spw=c18o_spw_7m, timebin='120s'
                      , outputvis=outputivs
                      , field=field_7m)
            
                os.system('rm -rf '+gal+'_c18o_bin.ms')
                concat(vis=[gal+'_12m_c18o_bin.ms',
                            gal+'_7m_c18o_bin.ms']
                       concatvis=gal+'_c18o_bin.ms')
                
        # Regrid C18O to a regular, 5km/s grid

        os.system('rm -rf '+gal+'_c18o21.ms')
        mstransform(vis=gal+'_c18o_bin.ms',
                    datacolumn='DATA',
                    outputvis=gal+'_c18o21.ms',
                    combinespws=True,
                    regridms=True,
                    mode='velocity',
                    width=c18o_deltav,
                    start=start_vel,
                    nchan=nchan/2,
                    restfreq='219.56035GHz',
                    outframe='lsrk',
                    veltype='radio',            
                    )

    # ------------------------------------------------
    # Continuum
    # ------------------------------------------------

    print "... Creating continuum and channel 0 measurement sets."

    if split_mosaic:

        counter = 1
        to_concat = []

        for this_file in calibrated_file:
                
            this_infile =  gal+'_'+str(counter)+'.ms'
            this_outfile = gal+'_'+str(counter)+'_cont_temp.ms'
            
            os.system('rm -rf '+this_outfile)
            os.system('rm -rf '+this_outfile+'.flagversions')
            split(vis=this_infile
                  , field=field
                  , datacolumn='DATA'
                  , width=5
                  , timebin='120s'
                  , outputvis=this_outfile)
            
            flagdata(vis=this_outfile,
                     spw=flag_co21)

            if has_c18o:
                flagdata(vis=this_outfile,
                         spw=flag_c18o21)

            to_concat.append(this_outfile)

            counter += 1

        if has_7m:
                
            counter = 1
            for this_file in calibrated_7m_file:
                
                this_infile =  gal+'_'+str(counter)+'_7m.ms'
                this_outfile = gal+'_'+str(counter)+'_7m_cont_temp.ms'
                
                os.system('rm -rf '+this_outfile)
                split(vis=this_infile
                      , field=field_7m
                      , datacolumn='DATA'
                      , width=5
                      , timebin='120s'
                      , outputvis=this_outfile)
            
                flagdata(vis=this_outfile,
                         spw=flag_co21_7m)

                if has_c18o:
                    flagdata(vis=this_outfile,
                             spw=flag_c18o21_7m)
            
                to_concat.append(this_outfile)

                counter += 1
        
        os.system('rm -rf '+gal+'_cont_temp.ms')
        concat(vis=to_concat,
               concatvis=gal+'_cont_temp.ms')

    else:

        print "... ... treating this as a single mosaic."

        os.system('rm -rf '+gal+'_cont_temp.ms*')
        split(vis=gal+'.ms'
              , field=field 
              , datacolumn='DATA'
              , width=5, timebin='120s'
              , outputvis=gal+'_cont_temp.ms')

        # Flag the CO lines in the continuum
        flagdata(vis=gal+'_cont_temp.ms',
                 spw=flag_co21)
        if has_c18o:
            flagdata(vis=gal+'_cont_temp.ms',
                     spw=flag_c18o21)

    # Average to only one channel per spectral window.
    os.system('rm -rf '+gal+'_cont.ms')
    split(vis=gal+'_cont_temp.ms' 
          , datacolumn='DATA'
          , width=4000 
          , outputvis=gal+'_cont.ms')

    # These two commands allow you to qickly determine by eye what
    # channel range should be flagged in order to blank out 12CO 2-1
    # and C180 1-0 so that they don't contaminate the continuum
    # image. Do them interactively.

    if False:
        plotms(vis=gal+'_cont_temp.ms',
               transform=True,
               freqframe='LSRK',
               restfreq='230.53800GHz',
               xaxis='velocity',
               yaxis='channel',
               spw=co21_spw)

    if False:
        plotms(vis=gal+'_cont_temp.ms',
               transform=True,
               freqframe='LSRK',
               restfreq='219.56035GHz',
               xaxis='velocity',
               yaxis='channel',
               spw=c18o_spw)

    # ......................................
    # Make "Channel 0" CO measurement sets
    # ......................................

    os.system('rm -rf '+gal+'_chan0_co21.ms')
    split(vis=gal+'_co21.ms' 
          , datacolumn='DATA'
          , width=4000, timebin='120s'
          , outputvis=gal+'_chan0_co21.ms')

    if has_c18o:
        os.system('rm -rf '+gal+'_chan0_c18o21.ms')
        split(vis=gal+'_c18o21.ms' 
              , datacolumn='DATA'
              , width=4000, timebin='120s'
              , outputvis=gal+'_chan0_c18o21.ms')

    # ------------------------------------------------
    # Delete intermediate steps
    # ------------------------------------------------

    if do_cleanup:

        print "Cleaning up intermediate files."
        
        os.system('rm -rf '+gal+'_co21_bin.ms')
        os.system('rm -rf '+gal+'_c18o_bin.ms')
        os.system('rm -rf '+gal+'_cont_temp.ms')

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# MAKE CLEAN MASKS FOR LATER IMAGING
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

if do_mask:

    print "--------------------------------------"
    print "(2) Making clean masks."
    print "--------------------------------------"
        
    print "... imaging CO 2-1 with a u-v taper and no mask."
    
    print "... ... making a dirty map."

    os.system('rm -rf '+gal+'_co21_temp*')
    tclean(vis=gal+'_co21.ms',
           imagename=gal+'_co21_temp',
           phasecenter=phase_center,
           gridder='mosaic',
           deconvolver='clark',
           cell='0.5arcsec',
           imsize=imsize_lowres,
           weighting='briggs',
           robust=0.5,
           specmode='cube',
           restfreq='230.53800GHz',
           outframe='lsrk',
           veltype='radio',
           niter=0,
           threshold='0Jy/beam',
           interactive=False,
           usemask='pb', 
           pbmask=0.2,
           uvtaper=['3.0arcsec','3.0arcsec','0deg']
           )

#   Figure out the RMS in the residual image

    print "... ... calculating a threshold. Will clean to S/N 3 at this stage."

    print "... ... finding a threshold."
    
    os.system('rm -rf '+gal+'_co21_pbmask.image')
    immath(imagename = gal+'_co21_temp.pb',
           outfile = gal+'_co21_pbmask.image',
           expr = 'iif(IM0 > 0.2,1.0,0.0)')
        
    cube_stat = imstat(imagename=gal+'_co21_temp.residual', 
                       mask=gal+'_co21_pbmask.image')
    cube_rms = cube_stat['medabsdevmed'][0]/0.6745
    mask_thresh = str(3.0*cube_rms)+"Jy/beam"
    print "... ... I calculated a threshold of "+mask_thresh

#   Resume CLEANing targeting a S/N = 3 threshold

    print "... ... proceeding with clean."

    tclean(vis=gal+'_co21.ms',
           imagename=gal+'_co21_temp',
           phasecenter=phase_center,
           gridder='mosaic',
           deconvolver='clark',
           cell='0.5arcsec',
           imsize=imsize_lowres,
           weighting='briggs',
           robust=0.5,
           specmode='cube',
           restfreq='230.53800GHz',
           outframe='lsrk',
           veltype='radio',
           niter=500000,
           threshold=mask_thresh,
           interactive=False,
           usemask='pb', 
           pbmask=0.2,
           uvtaper=['3.0arcsec','3.0arcsec','0deg']
           )

    # Smooth
    
    print "... smoothing to a 4x4 arcsecond beam."

    os.system('rm -rf '+gal+'_co21_formask.image')
    imsmooth(imagename=gal+'_co21_temp.image',
             outfile=gal+'_co21_formask.image',
             targetres=True,
             major='4.0arcsec', minor='4.0arcsec', pa='0deg',
             overwrite=True)    
    
    # Derive statistics

    print "... calculating image statistics."

    cube_stat = imstat(imagename=gal+'_co21_formask.image', 
                       mask=gal+'_co21_temp.mask')
    cube_rms = cube_stat['medabsdevmed'][0]/0.6745
    high_thresh = cube_rms*5.0
    low_thresh = cube_rms*2.0

    # Make a simple threshold mask

    print "... thresholding the image."

    os.system('rm -rf '+gal+'_co21_mask.image')
    immath(imagename = gal+'_co21_formask.image',
           outfile = gal+'_co21_mask.image',
           expr = 'iif(IM0 > '+str(high_thresh) +',1.0,0.0)')

    os.system('rm -rf '+gal+'_co21_lowmask.image')
    immath(imagename = gal+'_co21_formask.image',
           outfile = gal+'_co21_lowmask.image',
           expr = 'iif(IM0 > '+str(low_thresh) +',1.0,0.0)')
    
    # Manipulate the mask in more detail

    print "... processing the mask."

    import scipy
    import numpy as np
    
    print "... ... require two adjacent channels."
    ia.open(gal+'_co21_mask.image')
    mask = ia.getchunk()    
    new_mask = mask*(np.roll(mask,1,3) + np.roll(mask,-1,3)) >= 1
    ia.putchunk(new_mask*1.0)
    ia.done()

    ia.open(gal+'_co21_lowmask.image')
    mask = ia.getchunk()    
    new_mask = mask*(np.roll(mask,1,3) + np.roll(mask,-1,3)) >= 1
    ia.putchunk(new_mask*1.0)
    ia.done()

    print "... ... require a high significance seed."

    ia.open(gal+'_co21_lowmask.image')
    low_mask = ia.getchunk()    
    ia.done()

    ia.open(gal+'_co21_mask.image')
    mask = ia.getchunk()    
    regions, n_regions = scipy.ndimage.label(mask)                     
    myhistogram = scipy.ndimage.measurements.histogram(regions,0,n_regions+1,n_regions+1)
    object_slices = scipy.ndimage.find_objects(regions)
    threshold=25
    for i in range(n_regions):
        if myhistogram[i+1] < threshold:
            mask[object_slices[i]] = 0

    new_mask = scipy.ndimage.binary_dilation(mask,iterations=-1,mask=low_mask)
    
    structure = scipy.ndimage.generate_binary_structure(4,1)
    structure[:,:,:,0] = False    
    structure[:,:,:,2] = False
    new_mask = scipy.ndimage.binary_dilation(mask,structure,iterations=2,mask=low_mask)
    ia.putchunk(new_mask*1.0)
    ia.done()

    # Collapse the mask into a two dimensional version

    os.system('rm -rf '+gal+'_mask_twod.image')
    immoments(imagename=gal+'_co21_mask.image',
              moments=8,
              outfile=gal+'_mask_twod.image')
    
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# IMAGE THE CO 2-1
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

if do_co21:

    print "--------------------------------------"
    print "(3) Imaging CO 2-1 emission."
    print "--------------------------------------"

    if do_taper:

        # ------------------------------------------------
        # Image with a U-V taper
        # ------------------------------------------------
        
        # Non interactive clean
        
        print "... imaging CO 2-1 with a u-v taper."

        # Make a dirty map to get the astrometry that we use for the
        # mask and figure out the RMS in the cube so that we can
        # figure out the target threshold.

        print "... ... making a dirty map."

        os.system('rm -rf '+gal+'_co21_taper*')

        tclean(vis=gal+'_co21.ms',
               imagename=gal+'_co21_taper',
               phasecenter=phase_center,
               gridder='mosaic',
               deconvolver='clark',
               cell='0.5arcsec',
               imsize=imsize_lowres,
               weighting='briggs',
               robust=0.5,
               specmode='cube',
               restfreq='230.53800GHz',
               outframe='lsrk',
               veltype='radio',
               niter=0,
               threshold="0mJy/beam",
               interactive=False,
               usemask='pb', 
               pbmask=0.2,
               uvtaper=['2.5arcsec','2.5arcsec','0deg']
               )

        # Align the mask that we made above and figure out the
        # statistics of the cube.

        print "... ... aligning the mask to the map."

        os.system('rm -rf '+gal+'_mask_for_taper.image')
        imregrid(imagename=gal+'_co21_mask.image'
                 , output=gal+'_mask_for_taper.image'
                 , template=gal+'_co21_taper.residual'
                 , interpolation='nearest')
        
        print "... ... finding a threshold."

        os.system('rm -rf '+gal+'_co21_pbmask.image')
        immath(imagename = gal+'_co21_taper.pb',
               outfile = gal+'_co21_pbmask.image',
               expr = 'iif(IM0 > 0.2,1.0,0.0)')
        
        cube_stat = imstat(imagename=gal+'_co21_taper.residual', 
                           mask=gal+'_co21_pbmask.image')
        cube_rms = cube_stat['medabsdevmed'][0]/0.6745
        thresh_string = str(thresh_factor*cube_rms)+"Jy/beam"
        print "... ... I calculated a threshold of "+thresh_string

        # Proceed with the clean using the mask and threshold 
        
        print "... ... proceeding with clean."

        os.system('rm -rf '+gal+'_co21_taper.mask')

        tclean(vis=gal+'_co21.ms',
               imagename=gal+'_co21_taper',
               phasecenter=phase_center,
               gridder='mosaic',
               deconvolver='clark',
               cell='0.5arcsec',
               imsize=imsize_lowres,
               weighting='briggs',
               robust=0.5,
               specmode='cube',
               restfreq='230.53800GHz',
               outframe='lsrk',
               veltype='radio',
               niter=1000000,
               threshold=thresh_string,
               interactive=False,
               usemask='user', 
               mask = gal+'_mask_for_taper.image',
               uvtaper=['2.5arcsec','2.5arcsec','0deg']
               )

        print "... post processing and exporting the tapered CO 2-1."

        # Smooth to a round beam
        
        imsmooth(imagename=gal+'_co21_taper.image',
                 outfile=gal+'_co21_taper_round.image',
                 targetres=True,
                 major='3.0arcsec', minor='3.0arcsec', pa='0deg',
                 overwrite=True)
        
        # Correct for the primary beam

        os.system('rm -rf '+gal+'_co21_taper_pbcor.image')
        impbcor(imagename=gal+'_co21_taper.image',
                pbimage=gal+'_co21_taper.pb',
                outfile=gal+'_co21_taper_pbcor.image')
        
        # Smooth the primary beam corrected image

        imsmooth(imagename=gal+'_co21_taper_pbcor.image',
                 outfile=gal+'_co21_taper_round_pbcor.image',
                 targetres=True,
                 major='3.0arcsec', minor='3.0arcsec', pa='0deg',
                 overwrite=True)

        # Export to FITS

        exportfits(imagename=gal+'_co21_taper_round_pbcor.image',
                   fitsimage=gal+'_co21_taper_round_pbcor.fits',
                   velocity=True, overwrite=True, dropstokes=True, dropdeg=True, bitpix=16)

        exportfits(imagename=gal+'_co21_taper.residual',
                   fitsimage=gal+'_co21_taper_residual.fits',
                   velocity=True, overwrite=True, dropstokes=True, dropdeg=True, bitpix=16)
        
        exportfits(imagename=gal+'_co21_taper_round.image',
                   fitsimage=gal+'_co21_taper_round.fits',
                   velocity=True, overwrite=True, dropstokes=True, dropdeg=True, bitpix=16)
        
        exportfits(imagename=gal+'_co21_taper.pb',
                   fitsimage=gal+'_co21_taper_pb.fits',
                   velocity=True, overwrite=True, dropstokes=True, dropdeg=True, bitpix=16)

    if do_native:

        # ------------------------------------------------
        # Image at the native resolution
        # ------------------------------------------------
        
        # Run a non-interactive clean
        
        print "... imaging CO 2-1 at the native resolution."

        # Make a dirty map to get the astrometry that we use for the
        # mask and figure out the RMS in the cube so that we can
        # figure out the target threshold.

        print "... ... making a dirty map."
        
        os.system('rm -rf '+gal+'_co21_cube*')

        if skip_mask:
            print " ... ... I will clean without a mask down to about 3 sigma."
            tclean(vis=gal+'_co21.ms',
                   imagename=gal+'_co21_cube',
                   phasecenter=phase_center,
                   gridder='mosaic',
                   deconvolver='clark',
                   cell='0.15arcsec',
                   imsize=imsize,
                   weighting='briggs',
                   robust=0.5,
                   specmode='cube',
                   restfreq='230.53800GHz' ,
                   outframe='lsrk',
                   veltype='radio',
                   niter=1000000,
                   threshold='0.018Jy/beam',
                   interactive=False,
                   usemask='pb', 
                   pbmask=0.2,
                   )
        else:            
            tclean(vis=gal+'_co21.ms',
                   imagename=gal+'_co21_cube',
                   phasecenter=phase_center,
                   gridder='mosaic',
                   deconvolver='clark',
                   cell='0.15arcsec',
                   imsize=imsize,
                   weighting='briggs',
                   robust=0.5,
                   specmode='cube',
                   restfreq='230.53800GHz' ,
                   outframe='lsrk',
                   veltype='radio',
                   niter=0,
                   threshold='0mJy/beam',
                   interactive=False,
                   usemask='pb', 
                   pbmask=0.2,
                   )

            # Align the mask that we made above and figure out the
            # statistics of the cube.
            
            print "... ... aligning the mask to the map."
        
            os.system('rm -rf '+gal+'_mask_for_native.image')
            imregrid(imagename=gal+'_co21_mask.image'
                     , output=gal+'_mask_for_native.image'
                     , template=gal+'_co21_cube.residual'
                     , interpolation='nearest')
            
            print "... ... finding a threshold."

            os.system('rm -rf '+gal+'_co21_pbmask.image')
            immath(imagename = gal+'_co21_cube.pb',
                   outfile = gal+'_co21_pbmask.image',
                   expr = 'iif(IM0 > 0.2,1.0,0.0)')
            
            cube_stat = imstat(imagename=gal+'_co21_cube.residual', 
                               mask=gal+'_co21_pbmask.image')
            cube_rms = cube_stat['medabsdevmed'][0]/0.6745
            thresh_string = str(thresh_factor*cube_rms)+"Jy/beam"
            
            print "I found a threshold of "+thresh_string
            
            # Proceed with the clean using the mask and threshold 
            
            print "... ... proceeding with clean."

            os.system('rm -rf '+gal+'_co21_cube.mask')

            tclean(vis=gal+'_co21.ms',
                   imagename=gal+'_co21_cube',
                   phasecenter=phase_center,
                   gridder='mosaic',
                   deconvolver='clark',
                   cell='0.15arcsec',
                   imsize=imsize,
                   weighting='briggs',
                   robust=0.5,
                   specmode='cube',
                   restfreq='230.53800GHz',
                   outframe='lsrk',
                   veltype='radio',
                   niter=10000000,
                   threshold=thresh_string,
                   interactive=False,
                   usemask='user', 
                   mask=gal+'_mask_for_native.image',
                   )

        print "... post processing and exporting the native resolution CO 2-1."

        # Smooth to have a round beam
        
        imsmooth(imagename=gal+'_co21_cube.image',
                 outfile=gal+'_co21_cube_round.image',
                 targetres=True,
                 major=target_beam_co21, minor=target_beam_co21, pa='0deg',
                 overwrite=True)
        
        # Primary beam correct
        
        os.system('rm -rf '+gal+'_co21_cube_pbcor.image')
        impbcor(imagename=gal+'_co21_cube.image',
                pbimage=gal+'_co21_cube.pb',
                outfile=gal+'_co21_cube_pbcor.image')
        
        # Smooth the corrected image to have a round beam
        
        imsmooth(imagename=gal+'_co21_cube_pbcor.image',
                 outfile=gal+'_co21_cube_round_pbcor.image',
                 targetres=True,
                 major=target_beam_co21, minor=target_beam_co21, pa='0deg',
                 overwrite=True)
        
        # Shrink the cube to save space
    
        os.system('rm -rf '+gal+'_co21_cube_round_pbcor_rebin.image')
        imrebin(imagename=gal+'_co21_cube_round_pbcor.image',
                outfile=gal+'_co21_cube_round_pbcor_rebin.image',
                factor=[2,2,1,1])

        os.system('rm -rf '+gal+'_co21_cube_rebin.residual')
        imrebin(imagename=gal+'_co21_cube.residual',
                outfile=gal+'_co21_cube_rebin.residual',
                factor=[2,2,1,1])
    
        os.system('rm -rf '+gal+'_co21_cube_round_rebin.image')
        imrebin(imagename=gal+'_co21_cube_round.image',
                outfile=gal+'_co21_cube_round_rebin.image',
                factor=[2,2,1,1])
    
        os.system('rm -rf '+gal+'_co21_cube_rebin.pb')
        imrebin(imagename=gal+'_co21_cube.pb',
                outfile=gal+'_co21_cube_rebin.pb',
                factor=[2,2,1,1])

        # Export to FITS
        
        exportfits(imagename=gal+'_co21_cube_round_pbcor_rebin.image',
                   fitsimage=gal+'_co21_cube_round_pbcor.fits',
                   velocity=True, overwrite=True, dropstokes=True, dropdeg=True, bitpix=16)

        exportfits(imagename=gal+'_co21_cube_rebin.residual',
                   fitsimage=gal+'_co21_cube_residual.fits',
                   velocity=True, overwrite=True, dropstokes=True, dropdeg=True, bitpix=16)
        
        exportfits(imagename=gal+'_co21_cube_round_rebin.image',
                   fitsimage=gal+'_co21_cube_round.fits',
                   velocity=True, overwrite=True, dropstokes=True, dropdeg=True, bitpix=16)
        
        exportfits(imagename=gal+'_co21_cube_rebin.pb',
                   fitsimage=gal+'_co21_cube_pb.fits',
                   velocity=True, overwrite=True, dropstokes=True, dropdeg=True, bitpix=16)
    
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# IMAGE THE C18O 2-1
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

if do_c18o:

    print "--------------------------------------"
    print "(4) Imaging C18O 2-1 emission."
    print "--------------------------------------"

    if do_taper:

        # ------------------------------------------------
        # Image with a U-V taper
        # ------------------------------------------------
        
        # Non interactive clean
        
        print "... imaging C18O 2-1 with a u-v taper."

        # Make a dirty map to get the astrometry that we use for the
        # mask and figure out the RMS in the cube so that we can
        # figure out the target threshold.

        print "... ... making a dirty map."

        os.system('rm -rf '+gal+'_c18o21_taper*')

        tclean(vis=gal+'_c18o21.ms',
               imagename=gal+'_c18o21_taper',
               phasecenter=phase_center,
               gridder='mosaic',
               deconvolver='clark',
               cell='0.5arcsec',
               imsize=imsize_lowres,
               weighting='briggs',
               robust=0.5,
               specmode='cube',
               restfreq='219.56035GHz',
               outframe='lsrk',
               veltype='radio',
               niter=0,
               threshold="0mJy/beam",
               interactive=False,
               usemask='pb', 
               pbmask=0.2,
               uvtaper=['2.5arcsec','2.5arcsec','0deg']
               )

        # Align the mask that we made above and figure out the
        # statistics of the cube.

        print "... ... aligning the mask to the map."

        os.system('rm -rf '+gal+'_mask_for_taper_temp.image')
        imregrid(imagename=gal+'_co21_mask.image'
                 , output=gal+'_mask_for_taper_temp.image'
                 , template=gal+'_c18o21_taper.residual'
                 , interpolation='nearest'
                 , asvelocity=True)
        
        os.system('rm -rf '+gal+'_mask_for_taper.image')
        immath(imagename = gal+'_c18o21_taper.pb',
               outfile = gal+'_mask_for_taper.image',
               expr = 'iif(IM0 > 0.2,1.0,0.0)')
        
        ia.open(gal+'_mask_for_taper_temp.image')
        mask = ia.getchunk()
        ia.close()

        ia.open(gal+'_mask_for_taper.image')
        ia.putchunk(mask)
        ia.close()
        
        print "... ... finding a threshold."

        os.system('rm -rf '+gal+'_c18o21_pbmask.image')
        immath(imagename = gal+'_c18o21_taper.pb',
               outfile = gal+'_c18o21_pbmask.image',
               expr = 'iif(IM0 > 0.2,1.0,0.0)')
        
        cube_stat = imstat(imagename=gal+'_c18o21_taper.residual', 
                           mask=gal+'_c18o21_pbmask.image')
        cube_rms = cube_stat['medabsdevmed'][0]/0.6745
        thresh_string = str(thresh_factor*cube_rms)+"Jy/beam"
        print "... ... I calculated a threshold of "+thresh_string

        # Proceed with the clean using the mask and threshold 
        
        print "... ... proceeding with clean."

        os.system('rm -rf '+gal+'_c18o21_taper.mask')

        tclean(vis=gal+'_c18o21.ms',
               imagename=gal+'_c18o21_taper',
               phasecenter=phase_center,
               gridder='mosaic',
               deconvolver='clark',
               cell='0.5arcsec',
               imsize=imsize_lowres,
               weighting='briggs',
               robust=0.5,
               specmode='cube',
               restfreq='219.56035GHz',
               outframe='lsrk',
               veltype='radio',
               niter=1000000,
               threshold=thresh_string,
               interactive=False,
               usemask='user', 
               mask = gal+'_mask_for_taper.image',
               uvtaper=['2.5arcsec','2.5arcsec','0deg']
               )

        print "... post processing and exporting the tapered CO 2-1."

        # Smooth to a round beam
        
        imsmooth(imagename=gal+'_c18o21_taper.image',
                 outfile=gal+'_c18o21_taper_round.image',
                 targetres=True,
                 major='3.0arcsec', minor='3.0arcsec', pa='0deg',
                 overwrite=True)
        
        # Correct for the primary beam

        os.system('rm -rf '+gal+'_c18o21_taper_pbcor.image')
        impbcor(imagename=gal+'_c18o21_taper.image',
                pbimage=gal+'_c18o21_taper.pb',
                outfile=gal+'_c18o21_taper_pbcor.image')
        
        # Smooth the primary beam corrected image

        imsmooth(imagename=gal+'_c18o21_taper_pbcor.image',
                 outfile=gal+'_c18o21_taper_round_pbcor.image',
                 targetres=True,
                 major='3.0arcsec', minor='3.0arcsec', pa='0deg',
                 overwrite=True)

        # Export to FITS

        exportfits(imagename=gal+'_c18o21_taper_round_pbcor.image',
                   fitsimage=gal+'_c18o21_taper_round_pbcor.fits',
                   velocity=True, overwrite=True, dropstokes=True, dropdeg=True, bitpix=16)

        exportfits(imagename=gal+'_c18o21_taper.residual',
                   fitsimage=gal+'_c18o21_taper_residual.fits',
                   velocity=True, overwrite=True, dropstokes=True, dropdeg=True, bitpix=16)
        
        exportfits(imagename=gal+'_c18o21_taper_round.image',
                   fitsimage=gal+'_c18o21_taper_round.fits',
                   velocity=True, overwrite=True, dropstokes=True, dropdeg=True, bitpix=16)
        
        exportfits(imagename=gal+'_c18o21_taper.pb',
                   fitsimage=gal+'_c18o21_taper_pb.fits',
                   velocity=True, overwrite=True, dropstokes=True, dropdeg=True, bitpix=16)


    # ------------------------------------------------
    # Image at the native resolution
    # ------------------------------------------------

    if do_native:
        
        # Run a non-interactive clean
        
        print "... imaging C18O 2-1 at the native resolution."

        # Make a dirty map to get the astrometry that we use for the
        # mask and figure out the RMS in the cube so that we can
        # figure out the target threshold.

        print "... ... making a dirty map."
        
        os.system('rm -rf '+gal+'_c18o21_cube*')

        if skip_mask:
            print " ... ... I will clean without a mask down to about 3 sigma."
            tclean(vis=gal+'_c18o21.ms',
                   imagename=gal+'_c18o21_cube',
                   phasecenter=phase_center,
                   gridder='mosaic',
                   deconvolver='clark',
                   cell='0.15arcsec',
                   imsize=imsize,
                   weighting='briggs',
                   robust=0.5,
                   specmode='cube',
                   restfreq='219.56035GHz',
                   outframe='lsrk',
                   veltype='radio',
                   niter=1000000,
                   threshold='0.018Jy/beam',
                   interactive=False,
                   usemask='pb', 
                   pbmask=0.2,
                   )
        else:            
            tclean(vis=gal+'_c18o21.ms',
                   imagename=gal+'_c18o21_cube',
                   phasecenter=phase_center,
                   gridder='mosaic',
                   deconvolver='clark',
                   cell='0.15arcsec',
                   imsize=imsize,
                   weighting='briggs',
                   robust=0.5,
                   specmode='cube',
                   restfreq='219.56035GHz',
                   outframe='lsrk',
                   veltype='radio',
                   niter=0,
                   threshold='0mJy/beam',
                   interactive=False,
                   usemask='pb', 
                   pbmask=0.2,
                   )

            # Align the mask that we made above and figure out the
            # statistics of the cube.
            
            print "... ... aligning the mask to the map."
                    
            os.system('rm -rf '+gal+'_mask_for_native_temp.image')
            imregrid(imagename=gal+'_co21_mask.image'
                     , output=gal+'_mask_for_native_temp.image'
                     , template=gal+'_c18o21_cube.residual'
                     , interpolation='nearest'
                     , asvelocity=True)
            
            os.system('rm -rf '+gal+'_mask_for_native.image')
            immath(imagename = gal+'_c18o21_cube.pb',
                   outfile = gal+'_mask_for_native.image',
                   expr = 'iif(IM0 > 0.2,1.0,0.0)')
        
            ia.open(gal+'_mask_for_native_temp.image')
            mask = ia.getchunk()
            ia.close()

            ia.open(gal+'_mask_for_native.image')
            ia.putchunk(mask)
            ia.close()
 
            print "... ... finding a threshold."

            os.system('rm -rf '+gal+'_c18o21_pbmask.image')
            immath(imagename = gal+'_c18o21_cube.pb',
                   outfile = gal+'_c18o21_pbmask.image',
                   expr = 'iif(IM0 > 0.2,1.0,0.0)')
            
            cube_stat = imstat(imagename=gal+'_c18o21_cube.residual', 
                               mask=gal+'_c18o21_pbmask.image')
            cube_rms = cube_stat['medabsdevmed'][0]/0.6745
            thresh_string = str(thresh_factor*cube_rms)+"Jy/beam"
            
            print "I found a threshold of "+thresh_string
            
            # Proceed with the clean using the mask and threshold 
            
            print "... ... proceeding with clean."

            os.system('rm -rf '+gal+'_c18o21_cube.mask')

            tclean(vis=gal+'_c18o21.ms',
                   imagename=gal+'_c18o21_cube',
                   phasecenter=phase_center,
                   gridder='mosaic',
                   deconvolver='clark',
                   cell='0.15arcsec',
                   imsize=imsize,
                   weighting='briggs',
                   robust=0.5,
                   specmode='cube',
                   restfreq='219.56035GHz',
                   outframe='lsrk',
                   veltype='radio',
                   niter=10000000,
                   threshold=thresh_string,
                   interactive=False,
                   usemask='user', 
                   mask=gal+'_mask_for_native.image',
                   )

        print "... post processing and exporting the native resolution C18O 2-1."

        # Smooth to have a round beam
        
        imsmooth(imagename=gal+'_c18o21_cube.image',
                 outfile=gal+'_c18o21_cube_round.image',
                 targetres=True,
                 major=target_beam_c18o21, minor=target_beam_c18o21, pa='0deg',
                 overwrite=True)
        
        # Primary beam correct
        
        os.system('rm -rf '+gal+'_c18o21_cube_pbcor.image')
        impbcor(imagename=gal+'_c18o21_cube.image',
                pbimage=gal+'_c18o21_cube.pb',
                outfile=gal+'_c18o21_cube_pbcor.image')
        
        # Smooth the corrected image to have a round beam
        
        imsmooth(imagename=gal+'_c18o21_cube_pbcor.image',
                 outfile=gal+'_c18o21_cube_round_pbcor.image',
                 targetres=True,
                 major=target_beam_c18o21, minor=target_beam_c18o21, pa='0deg',
                 overwrite=True)
        
        # Shrink the cube to save space
    
        os.system('rm -rf '+gal+'_c18o21_cube_round_pbcor_rebin.image')
        imrebin(imagename=gal+'_c18o21_cube_round_pbcor.image',
                outfile=gal+'_c18o21_cube_round_pbcor_rebin.image',
                factor=[2,2,1,1])

        os.system('rm -rf '+gal+'_c18o21_cube_rebin.residual')
        imrebin(imagename=gal+'_c18o21_cube.residual',
                outfile=gal+'_c18o21_cube_rebin.residual',
                factor=[2,2,1,1])
    
        os.system('rm -rf '+gal+'_c18o21_cube_round_rebin.image')
        imrebin(imagename=gal+'_c18o21_cube_round.image',
                outfile=gal+'_c18o21_cube_round_rebin.image',
                factor=[2,2,1,1])
    
        os.system('rm -rf '+gal+'_c18o21_cube_rebin.pb')
        imrebin(imagename=gal+'_c18o21_cube.pb',
                outfile=gal+'_c18o21_cube_rebin.pb',
                factor=[2,2,1,1])

        # Export to FITS
        
        exportfits(imagename=gal+'_c18o21_cube_round_pbcor_rebin.image',
                   fitsimage=gal+'_c18o21_cube_round_pbcor.fits',
                   velocity=True, overwrite=True, dropstokes=True, dropdeg=True, bitpix=16)

        exportfits(imagename=gal+'_c18o21_cube_rebin.residual',
                   fitsimage=gal+'_c18o21_cube_residual.fits',
                   velocity=True, overwrite=True, dropstokes=True, dropdeg=True, bitpix=16)
        
        exportfits(imagename=gal+'_c18o21_cube_round_rebin.image',
                   fitsimage=gal+'_c18o21_cube_round.fits',
                   velocity=True, overwrite=True, dropstokes=True, dropdeg=True, bitpix=16)
        
        exportfits(imagename=gal+'_c18o21_cube_rebin.pb',
                   fitsimage=gal+'_c18o21_cube_pb.fits',
                   velocity=True, overwrite=True, dropstokes=True, dropdeg=True, bitpix=16)
        
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# CONTINUUM AND "CHANNEL0" 12CO AND C18O
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

if do_cont:

    print "--------------------------------------"
    print "(5) Imaging continuum and channel zero."
    print "--------------------------------------"

    if do_taper:

        # ------------------------------------------------
        # Image using a U-V taper
        # ------------------------------------------------

        print "... imaging continuum and channel 0 with a u-v taper."

        # ................................................
        # Image the continuum
        # ................................................

        print "... ... making a dirty map."        

        os.system('rm -rf '+gal+'_cont_taper*')
        tclean(vis=gal+'_cont.ms',
               imagename=gal+'_cont_taper',
               phasecenter=phase_center,
               gridder='mosaic',
               deconvolver='clark',
               cell='0.5arcsec',
               imsize=imsize_lowres,
               weighting='briggs',
               robust=0.5,
               specmode='mfs',
               niter=0,
               threshold='0.0Jy/beam',
               interactive=False,
               usemask='pb', 
               pbmask=0.2,
               uvtaper=['2.5arcsec','2.5arcsec','0deg']
               )

        # Align the mask that we made above and figure out the
        # statistics of the cube.

        print "... ... aligning the mask to the map."

        os.system('rm -rf '+gal+'_mask_for_taper_temp.image')
        imregrid(imagename=gal+'_mask_twod.image'
                 , output=gal+'_mask_for_taper_temp.image'
                 , template=gal+'_cont_taper.residual'
                 , interpolation='nearest'
                 , asvelocity=True)
        
        os.system('rm -rf '+gal+'_mask_for_taper.image')
        immath(imagename = gal+'_cont_taper.pb',
               outfile = gal+'_mask_for_taper.image',
               expr = 'iif(IM0 > 0.2,1.0,0.0)')
        
        ia.open(gal+'_mask_for_taper_temp.image')
        mask = ia.getchunk()
        ia.close()

        ia.open(gal+'_mask_for_taper.image')
        ia.putchunk(mask)
        ia.close()
        
        print "... ... finding a threshold."

        os.system('rm -rf '+gal+'_cont_pbmask.image')
        immath(imagename = gal+'_cont_taper.pb',
               outfile = gal+'_cont_pbmask.image',
               expr = 'iif(IM0 > 0.2,1.0,0.0)')
        
        cube_stat = imstat(imagename=gal+'_cont_taper.residual', 
                           mask=gal+'_cont_pbmask.image')
        cube_rms = cube_stat['medabsdevmed'][0]/0.6745
        thresh_string = str(thresh_factor*cube_rms)+"Jy/beam"
        print "... ... I calculated a threshold of "+thresh_string

        # Proceed with the clean using the mask and threshold 
        
        print "... ... proceeding with clean."

        os.system('rm -rf '+gal+'_cont_taper.mask')

        tclean(vis=gal+'_cont.ms',
               imagename=gal+'_cont_taper',
               phasecenter=phase_center,
               gridder='mosaic',
               deconvolver='clark',
               cell='0.5arcsec',
               imsize=imsize_lowres,
               weighting='briggs',
               robust=0.5,
               specmode='mfs',
               niter=100000,
               threshold=thresh_string,
               interactive=False,
               usemask='user', 
               mask=gal+'_mask_for_taper.image',
               uvtaper=['2.5arcsec','2.5arcsec','0deg']
               )

        # ................................................
        # Image the Channel 0 of the CO 2-1
        # ................................................

        print "... ... making a 12CO Channel 0 dirty map."        

        os.system('rm -rf '+gal+'_chan0_co21_taper*')
        tclean(vis=gal+'_chan0_co21.ms',
               imagename=gal+'_chan0_co21_taper',
               phasecenter=phase_center,
               gridder='mosaic',
               deconvolver='clark',
               cell='0.5arcsec',
               imsize=imsize_lowres,
               weighting='briggs',
               robust=0.5,
               specmode='mfs',
               niter=0,
               threshold='0.0Jy/beam',
               interactive=False,
               usemask='pb', 
               pbmask=0.2,
               uvtaper=['2.5arcsec','2.5arcsec','0deg']
               )

        # Align the mask that we made above and figure out the
        # statistics of the cube.

        print "... ... aligning the mask to the map."

        os.system('rm -rf '+gal+'_mask_for_taper_temp.image')
        imregrid(imagename=gal+'_mask_twod.image'
                 , output=gal+'_mask_for_taper_temp.image'
                 , template=gal+'_cont_taper.residual'
                 , interpolation='nearest'
                 , asvelocity=True)
        
        os.system('rm -rf '+gal+'_mask_for_taper.image')
        immath(imagename = gal+'_chan0_co21_taper.pb',
               outfile = gal+'_mask_for_taper.image',
               expr = 'iif(IM0 > 0.2,1.0,0.0)')
        
        ia.open(gal+'_mask_for_taper_temp.image')
        mask = ia.getchunk()
        ia.close()

        ia.open(gal+'_mask_for_taper.image')
        ia.putchunk(mask)
        ia.close()
        
        print "... ... finding a threshold."

        os.system('rm -rf '+gal+'_chan0_co21_pbmask.image')
        immath(imagename = gal+'_chan0_co21_taper.pb',
               outfile = gal+'_chan0_co21_pbmask.image',
               expr = 'iif(IM0 > 0.2,1.0,0.0)')
        
        cube_stat = imstat(imagename=gal+'_chan0_co21_taper.residual', 
                           mask=gal+'_chan0_co21_pbmask.image')
        cube_rms = cube_stat['medabsdevmed'][0]/0.6745
        thresh_string = str(thresh_factor*cube_rms)+"Jy/beam"
        print "... ... I calculated a threshold of "+thresh_string

        # Proceed with the clean using the mask and threshold 
        
        print "... ... proceeding with clean."

        os.system('rm -rf '+gal+'_chan0_co21_taper.mask')

        tclean(vis=gal+'_chan0_co21.ms',
               imagename=gal+'_chan0_co21_taper',
               phasecenter=phase_center,
               gridder='mosaic',
               deconvolver='clark',
               cell='0.5arcsec',
               imsize=imsize_lowres,
               weighting='briggs',
               robust=0.5,
               specmode='mfs',
               niter=100000,
               threshold=thresh_string,
               interactive=False,
               usemask='user', 
               mask=gal+'_mask_for_taper.image',
               uvtaper=['2.5arcsec','2.5arcsec','0deg']
               )

        if has_c18o:

            # ...........................................
            # Image the Channel 0 of the C18O 2-1
            # ...........................................

            print "... ... making a C18O Channel 0 dirty map."        

            os.system('rm -rf '+gal+'_chan0_c18o21_taper*')
            tclean(vis=gal+'_chan0_c18o21.ms',
                   imagename=gal+'_chan0_c18o21_taper',
                   phasecenter=phase_center,
                   gridder='mosaic',
                   deconvolver='clark',
                   cell='0.5arcsec',
                   imsize=imsize_lowres,
                   weighting='briggs',
                   robust=0.5,
                   specmode='mfs',
                   niter=0,
                   threshold='0.0Jy/beam',
                   interactive=False,
                   usemask='pb', 
                   pbmask=0.2,
                   uvtaper=['2.5arcsec','2.5arcsec','0deg']
                   )

            # Align the mask that we made above and figure out the
            # statistics of the cube.

            print "... ... aligning the mask to the map."
        
            os.system('rm -rf '+gal+'_mask_for_taper_temp.image')
            imregrid(imagename=gal+'_mask_twod.image'
                     , output=gal+'_mask_for_taper_temp.image'
                     , template=gal+'_cont_taper.residual'
                     , interpolation='nearest'
                     , asvelocity=True)
        
            os.system('rm -rf '+gal+'_mask_for_taper.image')
            immath(imagename = gal+'_chan0_c18o21_taper.pb',
                   outfile = gal+'_mask_for_taper.image',
                   expr = 'iif(IM0 > 0.2,1.0,0.0)')
            
            ia.open(gal+'_mask_for_taper_temp.image')
            mask = ia.getchunk()
            ia.close()
            
            ia.open(gal+'_mask_for_taper.image')
            ia.putchunk(mask)
            ia.close()
        
            print "... ... finding a threshold."

            os.system('rm -rf '+gal+'_chan0_c18o21_pbmask.image')
            immath(imagename = gal+'_chan0_c18o21_taper.pb',
                   outfile = gal+'_chan0_c18o21_pbmask.image',
                   expr = 'iif(IM0 > 0.2,1.0,0.0)')
        
            cube_stat = imstat(imagename=gal+'_chan0_c18o21_taper.residual', 
                               mask=gal+'_chan0_c18o21_pbmask.image')
            cube_rms = cube_stat['medabsdevmed'][0]/0.6745
            thresh_string = str(thresh_factor*cube_rms)+"Jy/beam"
            print "... ... I calculated a threshold of "+thresh_string

            # Proceed with the clean using the mask and threshold 
        
            print "... ... proceeding with clean."

            os.system('rm -rf '+gal+'_chan0_c18o21_taper.mask')
            
            tclean(vis=gal+'_chan0_c18o21.ms',
                   imagename=gal+'_chan0_c18o21_taper',
                   phasecenter=phase_center,
                   gridder='mosaic',
                   deconvolver='clark',
                   cell='0.5arcsec',
                   imsize=imsize_lowres,
                   weighting='briggs',
                   robust=0.5,
                   specmode='mfs',
                   niter=100000,
                   threshold=thresh_string,
                   interactive=False,
                   usemask='user', 
                   mask=gal+'_mask_for_taper.image',
                   uvtaper=['2.5arcsec','2.5arcsec','0deg']
                   )

        # Smooth to a round beam

        print "... ... processing and exporting continuum data."        

        imsmooth(imagename=gal+'_cont_taper.image',
                 outfile=gal+'_cont_taper_round.image',
                 targetres=True,
                 major='3.0arcsec', minor='3.0arcsec', pa='0deg',
                 overwrite=True)
        
        imsmooth(imagename=gal+'_chan0_co21_taper.image',
                 outfile=gal+'_chan0_co21_taper_round.image',
                 targetres=True,
                 major='3.0arcsec', minor='3.0arcsec', pa='0deg',
                 overwrite=True)

        if has_c18o:

            imsmooth(imagename=gal+'_chan0_c18o21_taper.image',
                     outfile=gal+'_chan0_c18o21_taper_round.image',
                     targetres=True,
                     major='3.0arcsec', minor='3.0arcsec', pa='0deg',
                     overwrite=True)

        # ... primary beam correct

        os.system('rm -rf '+gal+'_cont_taper_pbcor.image')
        impbcor(imagename=gal+'_cont_taper.image',
                pbimage=gal+'_cont_taper.pb',
                outfile=gal+'_cont_taper_pbcor.image')
        
        os.system('rm -rf '+gal+'_chan0_co21_taper_pbcor.image')
        impbcor(imagename=gal+'_chan0_co21_taper.image',
                pbimage=gal+'_chan0_co21_taper.pb',
                outfile=gal+'_chan0_co21_taper_pbcor.image')
        
        if has_c18o:

            os.system('rm -rf '+gal+'_chan0_c18o21_taper_pbcor.image')
            impbcor(imagename=gal+'_chan0_c18o21_taper.image',
                    pbimage=gal+'_chan0_c18o21_taper.pb',
                    outfile=gal+'_chan0_c18o21_taper_pbcor.image')
        
        # ... smooth the primary beam corrected image

        imsmooth(imagename=gal+'_cont_taper_pbcor.image',
                 outfile=gal+'_cont_taper_round_pbcor.image',
                 targetres=True,
                 major='3.0arcsec', minor='3.0arcsec', pa='0deg',
                 overwrite=True)

        imsmooth(imagename=gal+'_chan0_co21_taper_pbcor.image',
                 outfile=gal+'_chan0_co21_taper_round_pbcor.image',
                 targetres=True,
                 major='3.0arcsec', minor='3.0arcsec', pa='0deg',
                 overwrite=True)

        if has_c18o:

            imsmooth(imagename=gal+'_chan0_c18o21_taper_pbcor.image',
                     outfile=gal+'_chan0_c18o21_taper_round_pbcor.image',
                     targetres=True,
                     major='3.0arcsec', minor='3.0arcsec', pa='0deg',
                     overwrite=True)
        
        # ... export to FITS

        exportfits(imagename=gal+'_cont_taper_round_pbcor.image',
                   fitsimage=gal+'_cont_taper_round_pbcor.fits',
                   velocity=True, overwrite=True, dropstokes=True, dropdeg=True)
        
        exportfits(imagename=gal+'_cont_taper_round.image',
                   fitsimage=gal+'_cont_taper_round.fits',
                   velocity=True, overwrite=True, dropstokes=True, dropdeg=True)
        
        exportfits(imagename=gal+'_cont_taper.pb',
                   fitsimage=gal+'_cont_taper_pb.fits',
                   velocity=True, overwrite=True, dropstokes=True, dropdeg=True)
        
        exportfits(imagename=gal+'_chan0_co21_taper_round_pbcor.image',
                   fitsimage=gal+'_chan0_co21_taper_round_pbcor.fits',
                   velocity=True, overwrite=True, dropstokes=True, dropdeg=True)
        
        exportfits(imagename=gal+'_chan0_co21_taper_round.image',
                   fitsimage=gal+'_chan0_co21_taper_round.fits',
                   velocity=True, overwrite=True, dropstokes=True, dropdeg=True)
    
        exportfits(imagename=gal+'_chan0_co21_taper.pb',
                   fitsimage=gal+'_chan0_co21_taper_pb.fits',
                   velocity=True, overwrite=True, dropstokes=True, dropdeg=True)
        
        if has_c18o:

            exportfits(imagename=gal+'_chan0_c18o21_taper_round_pbcor.image',
                       fitsimage=gal+'_chan0_c18o21_taper_round_pbcor.fits',
                       velocity=True, overwrite=True, dropstokes=True, dropdeg=True)
            
            exportfits(imagename=gal+'_chan0_c18o21_taper_round.image',
                       fitsimage=gal+'_chan0_c18o21_taper_round.fits',
                       velocity=True, overwrite=True, dropstokes=True, dropdeg=True)
            
            exportfits(imagename=gal+'_chan0_c18o21_taper.pb',
                       fitsimage=gal+'_chan0_c18o21_taper_pb.fits',
                       velocity=True, overwrite=True, dropstokes=True, dropdeg=True)
            
    if do_native:

        # ------------------------------------------------
        # Image at the native resolution
        # ------------------------------------------------

        print "... imaging continuum and channel 0 at the native resolution."

        # ................................................
        # Image the continuum
        # ................................................

        print "... ... making a dirty map."        

        os.system('rm -rf '+gal+'_cont_map*')
        tclean(vis=gal+'_cont.ms',
               imagename=gal+'_cont_map',
               phasecenter=phase_center,
               gridder='mosaic',
               deconvolver='clark',
               cell='0.15arcsec',
               imsize=imsize,
               weighting='briggs',
               robust=0.5,
               specmode='mfs',
               niter=0,
               threshold='0.0Jy/beam',
               interactive=False,
               usemask='pb', 
               pbmask=0.2,
               )

        # Align the mask that we made above and figure out the
        # statistics of the cube.

        print "... ... aligning the mask to the map."

        os.system('rm -rf '+gal+'_mask_for_native_temp.image')
        imregrid(imagename=gal+'_mask_twod.image'
                 , output=gal+'_mask_for_native_temp.image'
                 , template=gal+'_cont_map.residual'
                 , interpolation='nearest'
                 , asvelocity=True)
        
        os.system('rm -rf '+gal+'_mask_for_native.image')
        immath(imagename = gal+'_cont_map.pb',
               outfile = gal+'_mask_for_native.image',
               expr = 'iif(IM0 > 0.2,1.0,0.0)')
        
        ia.open(gal+'_mask_for_native_temp.image')
        mask = ia.getchunk()
        ia.close()

        ia.open(gal+'_mask_for_native.image')
        ia.putchunk(mask)
        ia.close()
        
        print "... ... finding a threshold."

        os.system('rm -rf '+gal+'_cont_pbmask.image')
        immath(imagename = gal+'_cont_map.pb',
               outfile = gal+'_cont_pbmask.image',
               expr = 'iif(IM0 > 0.2,1.0,0.0)')
        
        cube_stat = imstat(imagename=gal+'_cont_map.residual', 
                           mask=gal+'_cont_pbmask.image')
        cube_rms = cube_stat['medabsdevmed'][0]/0.6745
        thresh_string = str(thresh_factor*cube_rms)+"Jy/beam"
        print "... ... I calculated a threshold of "+thresh_string

        # Proceed with the clean using the mask and threshold 
        
        print "... ... proceeding with clean."

        os.system('rm -rf '+gal+'_cont_map.mask')

        tclean(vis=gal+'_cont.ms',
               imagename=gal+'_cont_map',
               phasecenter=phase_center,
               gridder='mosaic',
               deconvolver='clark',
               cell='0.15arcsec',
               imsize=imsize,
               weighting='briggs',
               robust=0.5,
               specmode='mfs',
               niter=100000,
               threshold=thresh_string,
               interactive=False,
               usemask='user', 
               mask=gal+'_mask_for_native.image',
               )

        # ................................................
        # Image the Channel 0 of the CO 2-1
        # ................................................

        print "... ... making a 12CO Channel 0 dirty map."        

        os.system('rm -rf '+gal+'_chan0_co21_map*')
        tclean(vis=gal+'_chan0_co21.ms',
               imagename=gal+'_chan0_co21_map',
               phasecenter=phase_center,
               gridder='mosaic',
               deconvolver='clark',
               cell='0.15arcsec',
               imsize=imsize,
               weighting='briggs',
               robust=0.5,
               specmode='mfs',
               niter=0,
               threshold='0.0Jy/beam',
               interactive=False,
               usemask='pb', 
               pbmask=0.2,
               uvtaper=['2.5arcsec','2.5arcsec','0deg']
               )

        # Align the mask that we made above and figure out the
        # statistics of the cube.

        print "... ... aligning the mask to the map."

        os.system('rm -rf '+gal+'_mask_for_native_temp.image')
        imregrid(imagename=gal+'_mask_twod.image'
                 , output=gal+'_mask_for_native_temp.image'
                 , template=gal+'_cont_map.residual'
                 , interpolation='nearest'
                 , asvelocity=True)
        
        os.system('rm -rf '+gal+'_mask_for_native.image')
        immath(imagename = gal+'_chan0_co21_map.pb',
               outfile = gal+'_mask_for_native.image',
               expr = 'iif(IM0 > 0.2,1.0,0.0)')
        
        ia.open(gal+'_mask_for_native_temp.image')
        mask = ia.getchunk()
        ia.close()

        ia.open(gal+'_mask_for_native.image')
        ia.putchunk(mask)
        ia.close()
        
        print "... ... finding a threshold."

        os.system('rm -rf '+gal+'_chan0_co21_pbmask.image')
        immath(imagename = gal+'_chan0_co21_map.pb',
               outfile = gal+'_chan0_co21_pbmask.image',
               expr = 'iif(IM0 > 0.2,1.0,0.0)')
        
        cube_stat = imstat(imagename=gal+'_chan0_co21_map.residual', 
                           mask=gal+'_chan0_co21_pbmask.image')
        cube_rms = cube_stat['medabsdevmed'][0]/0.6745
        thresh_string = str(thresh_factor*cube_rms)+"Jy/beam"
        print "... ... I calculated a threshold of "+thresh_string

        # Proceed with the clean using the mask and threshold 
        
        print "... ... proceeding with clean."

        os.system('rm -rf '+gal+'_chan0_co21_map.mask')

        tclean(vis=gal+'_chan0_co21.ms',
               imagename=gal+'_chan0_co21_map',
               phasecenter=phase_center,
               gridder='mosaic',
               deconvolver='clark',
               cell='0.15arcsec',
               imsize=imsize,
               weighting='briggs',
               robust=0.5,
               specmode='mfs',
               niter=100000,
               threshold=thresh_string,
               interactive=False,
               usemask='user', 
               mask=gal+'_mask_for_native.image',
               )

        if has_c18o:

            # ...........................................
            # Image the Channel 0 of the C18O 2-1
            # ...........................................

            print "... ... making a C18O Channel 0 dirty map."        

            os.system('rm -rf '+gal+'_chan0_c18o21_map*')
            tclean(vis=gal+'_chan0_c18o21.ms',
                   imagename=gal+'_chan0_c18o21_map',
                   phasecenter=phase_center,
                   gridder='mosaic',
                   deconvolver='clark',
                   cell='0.15arcsec',
                   imsize=imsize,
                   weighting='briggs',
                   robust=0.5,
                   specmode='mfs',
                   niter=0,
                   threshold='0.0Jy/beam',
                   interactive=False,
                   usemask='pb', 
                   pbmask=0.2,
                   )

            # Align the mask that we made above and figure out the
            # statistics of the cube.

            print "... ... aligning the mask to the map."
        
            os.system('rm -rf '+gal+'_mask_for_native_temp.image')
            imregrid(imagename=gal+'_mask_twod.image'
                     , output=gal+'_mask_for_native_temp.image'
                     , template=gal+'_cont_map.residual'
                     , interpolation='nearest'
                     , asvelocity=True)
        
            os.system('rm -rf '+gal+'_mask_for_native.image')
            immath(imagename = gal+'_chan0_c18o21_map.pb',
                   outfile = gal+'_mask_for_native.image',
                   expr = 'iif(IM0 > 0.2,1.0,0.0)')
            
            ia.open(gal+'_mask_for_native_temp.image')
            mask = ia.getchunk()
            ia.close()
            
            ia.open(gal+'_mask_for_native.image')
            ia.putchunk(mask)
            ia.close()
        
            print "... ... finding a threshold."

            os.system('rm -rf '+gal+'_chan0_c18o21_pbmask.image')
            immath(imagename = gal+'_chan0_c18o21_map.pb',
                   outfile = gal+'_chan0_c18o21_pbmask.image',
                   expr = 'iif(IM0 > 0.2,1.0,0.0)')
        
            cube_stat = imstat(imagename=gal+'_chan0_c18o21_map.residual', 
                               mask=gal+'_chan0_c18o21_pbmask.image')
            cube_rms = cube_stat['medabsdevmed'][0]/0.6745
            thresh_string = str(thresh_factor*cube_rms)+"Jy/beam"
            print "... ... I calculated a threshold of "+thresh_string

            # Proceed with the clean using the mask and threshold 
        
            print "... ... proceeding with clean."
            
            os.system('rm -rf '+gal+'_chan0_c18o21_map.mask')

            tclean(vis=gal+'_chan0_c18o21.ms',
                   imagename=gal+'_chan0_c18o21_map',
                   phasecenter=phase_center,
                   gridder='mosaic',
                   deconvolver='clark',
                   cell='0.15arcsec',
                   imsize=imsize,
                   weighting='briggs',
                   robust=0.5,
                   specmode='mfs',
                   niter=100000,
                   threshold=thresh_string,
                   interactive=False,
                   usemask='user', 
                   mask=gal+'_mask_for_native.image',
                   )

        # Smooth to a round beam

        print "... ... processing and exporting continuum data."        

        # Smooth to a round beam

        imsmooth(imagename=gal+'_cont_map.image',
                 outfile=gal+'_cont_round.image',
                 targetres=True,
                 major=target_beam_cont, minor=target_beam_cont, pa='0deg',
                 overwrite=True)
        
        imsmooth(imagename=gal+'_chan0_co21_map.image',
                 outfile=gal+'_chan0_co21_round.image',
                 targetres=True,
                 major=target_beam_cont, minor=target_beam_cont, pa='0deg',
                 overwrite=True)

        if has_c18o:

            imsmooth(imagename=gal+'_chan0_c18o21_map.image',
                     outfile=gal+'_chan0_c18o21_round.image',
                     targetres=True,
                     major=target_beam_cont, minor=target_beam_cont, pa='0deg',
                     overwrite=True)

        # ... primary beam correct

        os.system('rm -rf '+gal+'_cont_map_pbcor.image')
        impbcor(imagename=gal+'_cont_map.image',
                pbimage=gal+'_cont_map.pb',
                outfile=gal+'_cont_map_pbcor.image')
        
        os.system('rm -rf '+gal+'_chan0_co21_map_pbcor.image')
        impbcor(imagename=gal+'_chan0_co21_map.image',
                pbimage=gal+'_chan0_co21_map.pb',
                outfile=gal+'_chan0_co21_map_pbcor.image')
        
        if has_c18o:

            os.system('rm -rf '+gal+'_chan0_c18o21_map_pbcor.image')
            impbcor(imagename=gal+'_chan0_c18o21_map.image',
                    pbimage=gal+'_chan0_c18o21_map.pb',
                    outfile=gal+'_chan0_c18o21_map_pbcor.image')

        # ... smooth the primary beam corrected image
        
        imsmooth(imagename=gal+'_cont_map_pbcor.image',
                 outfile=gal+'_cont_round_pbcor.image',
                 targetres=True,
                 major=target_beam_cont, minor=target_beam_cont, pa='0deg',
                 overwrite=True)
        
        imsmooth(imagename=gal+'_chan0_co21_map_pbcor.image',
                 outfile=gal+'_chan0_co21_round_pbcor.image',
                 targetres=True,
                 major=target_beam_cont, minor=target_beam_cont, pa='0deg',
                 overwrite=True)

        if has_c18o:

            imsmooth(imagename=gal+'_chan0_c18o21_map_pbcor.image',
                     outfile=gal+'_chan0_c18o21_round_pbcor.image',
                     targetres=True,
                     major=target_beam_cont, minor=target_beam_cont, pa='0deg',
                     overwrite=True)
        
        # ... shrink the images to save space

        os.system('rm -rf '+gal+'_cont_round_pbcor_rebin.image')
        imrebin(imagename=gal+'_cont_round_pbcor.image',
                outfile=gal+'_cont_round_pbcor_rebin.image',
                factor=[2,2,1,1])
        
        os.system('rm -rf '+gal+'_cont_round_rebin.image')
        imrebin(imagename=gal+'_cont_round.image',
                outfile=gal+'_cont_round_rebin.image',
                factor=[2,2,1,1])
        
        os.system('rm -rf '+gal+'_cont_rebin.pb')
        imrebin(imagename=gal+'_cont_map.pb',
                outfile=gal+'_cont_rebin.pb',
                factor=[2,2,1,1])
        
        os.system('rm -rf '+gal+'_chan0_co21_round_pbcor_rebin.image')
        imrebin(imagename=gal+'_chan0_co21_round_pbcor.image',
                outfile=gal+'_chan0_co21_round_pbcor_rebin.image',
                factor=[2,2,1,1])
        
        os.system('rm -rf '+gal+'_chan0_co21_round_rebin.image')
        imrebin(imagename=gal+'_chan0_co21_round.image',
                outfile=gal+'_chan0_co21_round_rebin.image',
                factor=[2,2,1,1])
    
        os.system('rm -rf '+gal+'_chan0_co21_rebin.pb')
        imrebin(imagename=gal+'_chan0_co21_map.pb',
                outfile=gal+'_chan0_co21_rebin.pb',
                factor=[2,2,1,1])
       
        if has_c18o:
 
            os.system('rm -rf '+gal+'_chan0_c18o21_round_pbcor_rebin.image')
            imrebin(imagename=gal+'_chan0_c18o21_round_pbcor.image',
                    outfile=gal+'_chan0_c18o21_round_pbcor_rebin.image',
                    factor=[2,2,1,1])
            
            os.system('rm -rf '+gal+'_chan0_c18o21_round_rebin.image')
            imrebin(imagename=gal+'_chan0_c18o21_round.image',
                    outfile=gal+'_chan0_c18o21_round_rebin.image',
                    factor=[2,2,1,1])
        
            os.system('rm -rf '+gal+'_chan0_c18o21_rebin.pb')
            imrebin(imagename=gal+'_chan0_c18o21_map.pb',
                    outfile=gal+'_chan0_c18o21_rebin.pb',
                    factor=[2,2,1,1])

        # ... export to FITS

        exportfits(imagename=gal+'_cont_round_pbcor_rebin.image',
                   fitsimage=gal+'_cont_round_pbcor.fits',
                   velocity=True, overwrite=True, dropstokes=True, dropdeg=True)
        
        exportfits(imagename=gal+'_cont_round_rebin.image',
                   fitsimage=gal+'_cont_round.fits',
                   velocity=True, overwrite=True, dropstokes=True, dropdeg=True)
        
        exportfits(imagename=gal+'_cont_rebin.pb',
                   fitsimage=gal+'_cont_pb.fits',
                   velocity=True, overwrite=True, dropstokes=True, dropdeg=True)
        
        exportfits(imagename=gal+'_chan0_co21_round_pbcor_rebin.image',
                   fitsimage=gal+'_chan0_co21_round_pbcor.fits',
                   velocity=True, overwrite=True, dropstokes=True, dropdeg=True)
        
        exportfits(imagename=gal+'_chan0_co21_round_rebin.image',
                   fitsimage=gal+'_chan0_co21_round.fits',
                   velocity=True, overwrite=True, dropstokes=True, dropdeg=True)
        
        exportfits(imagename=gal+'_chan0_co21_rebin.pb',
                   fitsimage=gal+'_chan0_co21_pb.fits',
                   velocity=True, overwrite=True, dropstokes=True, dropdeg=True)

        if has_c18o:

            exportfits(imagename=gal+'_chan0_c18o21_round_pbcor_rebin.image',
                       fitsimage=gal+'_chan0_c18o21_round_pbcor.fits',
                       velocity=True, overwrite=True, dropstokes=True, dropdeg=True)
            
            exportfits(imagename=gal+'_chan0_c18o21_round_rebin.image',
                       fitsimage=gal+'_chan0_c18o21_round.fits',
                       velocity=True, overwrite=True, dropstokes=True, dropdeg=True)
            
            exportfits(imagename=gal+'_chan0_c18o21_rebin.pb',
                       fitsimage=gal+'_chan0_c18o21_pb.fits',
                       velocity=True, overwrite=True, dropstokes=True, dropdeg=True)

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# TIMER
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

stop_time = time.time()

elapsed_time = (stop_time - start_time)/60.
print "This run took "+str(elapsed_time)+" minutes"
