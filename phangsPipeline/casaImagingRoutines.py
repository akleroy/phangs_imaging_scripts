"""
Standalone routines related to CASA imaging.
"""

#region Imports and definitions

import os
import numpy as np
import pyfits # CASA has pyfits, not astropy
import glob

import logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

# Analysis utilities
import analysisUtils as au

# CASA stuff
import casaStuff

from casaMaskingRoutines import signal_mask
from casaMaskingRoutines import stat_clean_cube

from clean_call import CleanCall

# Pipeline versionining
from pipelineVersion import version as pipeVer

#endregion

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# Routines to set up imaging
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

#region Setting up imaging

def estimate_cell_and_imsize(
    infile=None,    
    oversamp=5,
    forceSquare=False,
    ):
    """
    Pick a cell and image size for a measurement set. Requests an
    oversampling factor, which is by default 5. Will pick a good size
    for the FFT and will try to pick a round number for the cell size.
    """

    if os.path.isdir(infile) == False:
        logger.error('Error! The input file "'+infile+'" was not found!')
        return
    
    # These are the CASA-preferred sizes for fast FFTs

    valid_sizes = []
    for ii in range(10):
        for kk in range(3):
            for jj in range(3):
                valid_sizes.append(2**(ii+1)*5**(jj)*3**(kk))
    valid_sizes.sort()
    valid_sizes = np.array(valid_sizes)

    # Cell size implied by baseline distribution from analysis
    # utilities.

    au_cellsize, au_imsize, au_centralField = \
        au.pickCellSize(infile, imsize=True, npix=oversamp)
    xextent = au_cellsize*au_imsize[0]*1.2
    yextent = au_cellsize*au_imsize[1]*1.2

    # Make the cell size a nice round number

    if au_cellsize < 0.1:
        cell_size = au_cellsize
    if au_cellsize >= 0.1 and au_cellsize < 0.5:
        cell_size = np.floor(au_cellsize/0.05)*0.05
    if au_cellsize >= 0.5 and au_cellsize < 1.0:
        cell_size = np.floor(au_cellsize/0.1)*0.1
    if au_cellsize >= 1.0 and au_cellsize < 2.0:
        cell_size = np.floor(au_cellsize/0.25)*0.25
    if au_cellsize >= 2.0 and au_cellsize < 5.0:
        cell_size = np.floor(au_cellsize/0.5)*0.5
    if au_cellsize >= 5.0:
        cell_size = np.floor(au_cellsize/1.0)*0.5

    # Now make the image size a good number for the FFT

    need_cells_x = xextent / cell_size
    need_cells_y = yextent / cell_size

    cells_x = np.min(valid_sizes[valid_sizes > need_cells_x])
    cells_y = np.min(valid_sizes[valid_sizes > need_cells_y])

    # If requested, force the mosaic to be square. This avoids
    # pathologies in CASA versions 5.1 and 5.3.

    if forceSquare == True:
        if cells_y < cells_x:
            cells_y = cells_x
        if cells_x < cells_y:
            cells_x = cells_y

    image_size = [int(cells_x), int(cells_y)]
    cell_size_string = str(cell_size)+'arcsec'

    x_size_string = str(image_size[0])
    y_size_string = str(image_size[1])

    return cell_size_string, x_size_string, y_size_string

#endregion

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# Routines to set manipulate files associated with imaging
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

#region Input and output of imaging products

def wipe_imaging(
    image_root=None,
    ):
    """
    Wipe files associated with a cube or continuum imaging. Tries to
    delete all images and supporting products, including the output of
    any MFS imaging.
    """

    if image_root == None:
        return
    
    logger.debug('wipe_imaging')
    cmd_list = [
        'rm -rf '+image_root+'.image',
        'rm -rf '+image_root+'.alpha',
        'rm -rf '+image_root+'.beta',
        'rm -rf '+image_root+'.tt0',
        'rm -rf '+image_root+'.tt1',
        'rm -rf '+image_root+'.tt2',
        'rm -rf '+image_root+'.model',
        'rm -rf '+image_root+'.mask',
        'rm -rf '+image_root+'.pb',
        'rm -rf '+image_root+'.psf',
        'rm -rf '+image_root+'.residual',
        'rm -rf '+image_root+'.weight',
        'rm -rf '+image_root+'.sumwt',
        ]

    for this_cmd in cmd_list:
        logger.debug(this_cmd)
        os.system(this_cmd)

    return()

def save_copy_of_imaging(
    input_root=None,
    output_root=None):
    """
    Copy all of the files from a cube or continuum imaging output by
    clean to have a new root name. Most commonly used to make a backup
    copy of imaging output during iterative imaging (e.g., clean
    loops, shifting clean modes, selfcal, etc). Overwrites any
    previous imaging with that output name.
    """
    
    wipe_cube(output_root)
    
    logger.debug('save_copy_of_imaging')
    cmd_list = [
        'cp -r '+input_root+'.image '+output_root+'.image',
        'cp -r '+input_root+'.alpha '+output_root+'.alpha',
        'cp -r '+input_root+'.beta '+output_root+'.beta',
        'cp -r '+input_root+'.tt0 '+output_root+'.tt0',
        'cp -r '+input_root+'.tt1 '+output_root+'.tt1',
        'cp -r '+input_root+'.tt2 '+output_root+'.tt2',
        'cp -r '+input_root+'.model '+output_root+'.model',
        'cp -r '+input_root+'.mask '+output_root+'.mask',
        'cp -r '+input_root+'.pb '+output_root+'.pb',
        'cp -r '+input_root+'.psf '+output_root+'.psf',
        'cp -r '+input_root+'.residual '+output_root+'.residual',
        'cp -r '+input_root+'.weight '+output_root+'.weight',
        'cp -r '+input_root+'.sumwt '+output_root+'.sumwt',
        ]
    
    for this_cmd in cmd_list:
        logger.debug(this_cmd)
        os.system(this_cmd)

def replace_imaging_with_copy(
    input_root=None,
    output_root=None):
    """
    Replace a cube with a copy. Wipes any imaging associated with the
    target root root first.
    """
    
    wipe_cube(output_root)
    
    logger.debug('replace_imaging_with_copy')

    cmd_list = [
        'cp -r '+input_root+'.image '+output_root+'.image',
        'cp -r '+input_root+'.alpha '+output_root+'.alpha',
        'cp -r '+input_root+'.beta '+output_root+'.beta',
        'cp -r '+input_root+'.tt0 '+output_root+'.tt0',
        'cp -r '+input_root+'.tt1 '+output_root+'.tt1',
        'cp -r '+input_root+'.tt2 '+output_root+'.tt2',
        'cp -r '+input_root+'.mask '+output_root+'.mask',
        'cp -r '+input_root+'.pb '+output_root+'.pb',
        'cp -r '+input_root+'.psf '+output_root+'.psf',
        'cp -r '+input_root+'.residual '+output_root+'.residual',
        'cp -r '+input_root+'.weight '+output_root+'.weight',
        'cp -r '+input_root+'.sumwt '+output_root+'.sumwt',
        ]

    for this_cmd in cmd_list:
        logger.debug(this_cmd)
        os.system(this_cmd)

def export_imaging_to_fits(
    image_root=None,
    bitpix=-32,
    just_image=False):
    """
    Export the products associated with a CASA imaging run to FITS.
    """
    
    logger.debug('export_imaging_to_fits')

    ext_map = {
        '.image':'.fits',
        '.tt0':'.fits',
        '.tt1':'_tt1.fits',
        '.tt2':'_tt2.fits',
        '.alpha':'_alpha.fits',
        '.beta':'_beta.fits',
        '.mask':'_mask.fits',
        '.pb':'_pb.fits',
        '.psf':'_psf.fits',
        '.residual':'_residual.fits',
        '.weight':'_weight.fits',
        '.sumwt':'_sumwt.fits',
        }

    for this_ext in ext_map.keys():
        if just_image and ((this_ext != '.tt0') and this_ext != '.image'):
            continue

        this_casa_ext = this_ext
        this_fits_ext = ext_map[this_ext]

        casa_image = image_root + this_casa_ext
        if not os.pwd.isdir(casa_image):
            continue
        fits_image = image_root + this_fits_ext

        logger.debug('exportfits from '+casa_image+' to '+fits_image)
        casaStuff.exportfits(imagename=casa_image,
                             fitsimage=fits_image,
                             velocity=True, overwrite=True, dropstokes=True, 
                             dropdeg=True, bitpix=bitpix)
    
    return()

#endregion

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# Routines to actually execute the cleaning
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

#region clean call execution

def execute_clean_call(
    clean_call = None
    ):
    """
    Execute the clean call.
    """
    
    logger.debug('execute_clean_call')
    
    if not isinstance(clean_call, CleanCall):
        logger.error("Please input a valid clean call!")
        raise Exception("Please input a valid clean call!")
      
    if clean_call.vis == None:
        logger.info("No visibility. Returning.")
        return

    if os.path.isdir(clean_call.vis) == False:
        logger.info("Visibility file not found. Returning.")
        return
    
    if clean_call.cell_size == None or clean_call.image_size == None:
        logger.info("Estimating cell and image size.")
        cell_size, x_size, y_size = \
            estimate_cell_and_imsize(clean_call.vis, oversamp=5)
        clean_call.cell_size = cell_size
        clean_call.image_size = [x_size, y_size]

    if clean_call.logfile != None:
        oldlogfile = casaStuff.casalog.logfile()
        casaStuff.casalog.setlogfile(clean_call.logfile)

    if clean_call.reset:
        logger.info("Wiping previous versions of the cube.")
        wipe_cube(clean_call.image_root)
    
    logger.debug('Clean call: '+str(clean_call.kwargs()))
    casaStuff.tclean(**clean_call.kwargs_for_clean())

    if clean_call.logfile != None:
        casaStuff.casalog.setlogfile(oldlogfile)
    
def make_dirty_map(
    clean_call = None, 
    ):
    """
    Create a dirty map from a measurement set.
    """
    
    if not isinstance(clean_call, CleanCall):
        logger.error("Please input a valid clean call!")
        raise Exception("Please input a valid clean call!")
    
    clean_call.niter = 0
    clean_call.reset = True
    clean_call.usemask = 'pb'
    clean_call.logfile = clean_call.image_root+'_dirty.log'
    
    clean_call.calcres = True
    clean_call.calcpsf = True
    #clean_call.execute()
    execute_clean_call(clean_call)
    
    clean_call.reset = False
    clean_call.usemask = 'pb'
    clean_call.logfile = None
    
    save_copy_of_cube(
        input_root=clean_call.image_root,
        output_root=clean_call.image_root+'_dirty')










def clean_loop(
    clean_call = None,
    record_file=None,
    log_ext=None,
    delta_flux_threshold=0.02,    
    absolute_delta=True,
    absolute_threshold=None,
    snr_threshold=4.0,
    stop_at_negative=True,
    remask=False,
    max_loop = 20
    ):
    """
    Carry out an iterative clean until a convergence criteria is met.
    """
    
    if not isinstance(clean_call, CleanCall):
        logger.error("Please input a valid clean call!")
        raise Exception("Please input a valid clean call!")

    # Note the number of channels, which is used in setting the number
    # of iterations that we give to an individual clean call.

    vm = au.ValueMapping(clean_call.vis)
    nchan = vm.spwInfo[0]['numChannels']

    # Figure out the number of iterations we will use. Note that this
    # step is highly tunable, and can still be improved as we go
    # forward.

    base_niter = 10*nchan
    base_cycle_niter = 100
    loop = 1

    # Initialize our tracking of the flux in the model

    model_flux = 0.0

    # Open the text record if desired

    if record_file != None:
        f = open(record_file,'w')
        f.write("# column 1: loop type\n")
        f.write("# column 2: loop number\n")
        f.write("# column 3: supplied threshold\n")
        f.write("# column 4: model flux at end of this clean\n")
        f.write("# column 5: fractional change in flux (current-previous)/current\n")
        f.write("# column 6: number of iterations allocated (not necessarily used)\n")
        f.close()

    # Run the main loop

    proceed = True
    while proceed == True and loop <= max_loop:

        # Figure out how many iterations to give clean.

        if loop > 5:
            factor = 5
        else:
            factor = (loop-1)
        
        clean_call.niter = base_niter*(2**factor)
        clean_call.cycle_niter = base_cycle_niter*factor
        
        # Set the threshold for the clean call.

        if snr_threshold != None:
            resid_stats = stat_clean_cube(clean_call.image_root+'.residual')        
            current_noise = resid_stats['medabsdevmed'][0]/0.6745
            clean_call.threshold = str(current_noise*snr_threshold)+'Jy/beam'
        elif absolute_threshold != None:
            clean_call.threshold = absolute_threshold

        # If requested mask at each step (this is experimental, we're
        # seeing if it helps to avoid divergence during the deep
        # single scale clean.)

        if remask:
            logger.info("")
            logger.info("Remasking.")
            logger.info("")
            signal_mask(
                cube_root=clean_call.image_root,
                out_file=clean_call.image_root+'.mask',
                operation='AND',
                high_snr=4.0,
                low_snr=2.0,
                absolute=False)
            clean_call.usemask='user'

        logger.debug('clean_loop '+clean_call.image_root+' loop '+str(loop))
        
        # Set the log file

        if log_ext != None:
            clean_call.logfile = clean_call.image_root+"_loop_"+str(loop)+"_"+log_ext+".log"
        else:
            clean_call.logfile = None

        # Save the previous version of the file
        
        save_copy_of_cube(
            input_root=clean_call.image_root,
            output_root=clean_call.image_root+'_prev')

        # Execute the clean call.

        clean_call.reset = False
        #clean_call.execute()
        execute_clean_call(clean_call)

        clean_call.niter = 0
        clean_call.cycle_niter = 200

        # Record the new model flux and check for convergence. A nice
        # way to improve this would be to calculate the flux per
        # iteration.

        model_stats = stat_clean_cube(clean_call.image_root+'.model')

        prev_flux = model_flux
        model_flux = model_stats['sum'][0]

        delta_flux = (model_flux-prev_flux)/model_flux
        if absolute_delta:
            delta_flux = abs(delta_flux)

        if delta_flux_threshold >= 0.0:
            proceed = \
                (delta_flux > delta_flux_threshold)

        if stop_at_negative:
            if model_flux < 0.0:
                proceed = False
            
        # Print output
                
        logger.info("")
        logger.info("******************************")
        logger.info("CLEAN LOOP "+str(loop))
        logger.info("... threshold "+str(clean_call.threshold))
        logger.info("... old flux "+str(prev_flux))
        logger.info("... new flux "+str(model_flux))
        logger.info("... fractional change "+str(delta_flux)+" compare to stopping criterion of "+str(delta_flux_threshold))
        logger.info("... proceed? "+str(proceed))
        logger.info("******************************")
        logger.info("")

        # Record to log

        if record_file != None:
            line = 'LOOP '+str(loop)+ \
                ' '+clean_call.threshold+' '+str(model_flux)+ \
                ' '+str(delta_flux) + ' ' + str(clean_call.niter)+ '\n' 
            f = open(record_file,'a')
            f.write(line)
            f.close()

        if proceed == False:
            break
        loop += 1

    return



def singlescale_loop(
    clean_call = None,
    scales_as_angle=[],
    record_file=None,
    delta_flux_threshold=0.02,
    absolute_delta=True,
    absolute_threshold=None,
    snr_threshold=4.0,
    stop_at_negative=True,
    remask=False,
    max_loop = 20
    ):
    """
    Carry out an iterative multiscale clean loop.
    """
    
    # Check that we have a vile clean call

    if not isinstance(clean_call, CleanCall):
        logger.error("Please input a valid clean call!")
        raise Exception("Please input a valid clean call!")
        
    clean_call.deconvolver = 'hogbom'
    clean_call.calcres = False
    clean_call.calcpsf = False

    # Call the loop

    clean_loop(
        clean_call=clean_call,
        record_file=record_file,
        delta_flux_threshold=delta_flux_threshold,
        absolute_delta=absolute_delta,
        absolute_threshold=absolute_threshold,
        snr_threshold=snr_threshold,
        stop_at_negative=stop_at_negative,
        remask=remask,
        max_loop = max_loop, 
        log_ext = 'singlescale', #<TODO># set log_ext to output log files
        )

    # Save a copy

    save_copy_of_cube(
        input_root=clean_call.image_root,
        output_root=clean_call.image_root+'_singlescale')



def multiscale_loop(
    clean_call = None,
    record_file=None,
    delta_flux_threshold=0.02,
    absolute_delta=True,
    absolute_threshold=None,
    snr_threshold=4.0,
    stop_at_negative=True,
    max_loop = 20
    ):
    """
    Carry out an iterative multiscale clean loop.
    """
    
    # Check that we have a vile clean call

    if not isinstance(clean_call, CleanCall):
        logger.error("Please input a valid clean call!")
        raise Exception("Please input a valid clean call!")
    
    # Figure out the scales to use in pixel units

    cell_as_num = float((clean_call.cell_size.split('arcsec'))[0])
    scales_as_pix = []
    for scale in clean_call.scales_as_angle:
        scales_as_pix.append(int(scale/cell_as_num))
        
    clean_call.deconvolver = 'multiscale'
    clean_call.scales_as_pix = scales_as_pix
    clean_call.calcres = False
    clean_call.calcpsf = False

    logger.info("I will use the following scales: ")
    logger.info("... as pixels: " + str(clean_call.scales_as_pix))
    logger.info("... as arcseconds: " + str(clean_call.scales_as_angle))

    # Call the loop

    clean_loop(
        clean_call=clean_call,
        record_file=record_file,
        delta_flux_threshold=delta_flux_threshold,
        absolute_delta=absolute_delta,
        absolute_threshold=absolute_threshold,
        snr_threshold=snr_threshold,
        stop_at_negative=stop_at_negative,
        max_loop = max_loop, 
        log_ext = 'multiscale', #<TODO># set log_ext to output log files
        )

    # Save a copy

    save_copy_of_cube(
        input_root=clean_call.image_root,
        output_root=clean_call.image_root+'_multiscale')

#endregion






