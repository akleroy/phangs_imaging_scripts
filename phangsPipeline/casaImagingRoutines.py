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






#region Setting up imaging

def estimate_cell_and_imsize(
    in_file=None,    
    oversamp=5,
    forceSquare=False,
    ):
    """
    Pick a cell and image size for a measurement set. Requests an
    oversampling factor, which is by default 5. Will pick a good size
    for the FFT and will try to pick a round number for the cell size.
    """

    if os.path.isdir(in_file) == False:
        logger.error('Error! The input file "'+in_file+'" was not found!')
        return
    
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
        au.pickCellSize(in_file, imsize=True, npix=oversamp)
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

#region Input and output of imaging products

def wipe_imaging(
    image_root=None):
    """
    Wipe files associated with a cube.
    """
    if image_root == None:
        return
    
    logger.debug('wipe_imaging')
    logger.debug('rm -rf '+image_root+'.image')
    logger.debug('rm -rf '+image_root+'.model')
    logger.debug('rm -rf '+image_root+'.mask')
    logger.debug('rm -rf '+image_root+'.pb')
    logger.debug('rm -rf '+image_root+'.psf')
    logger.debug('rm -rf '+image_root+'.residual')
    logger.debug('rm -rf '+image_root+'.weight')
    logger.debug('rm -rf '+image_root+'.sumwt')
    
    os.system('rm -rf '+image_root+'.image')
    os.system('rm -rf '+image_root+'.model')
    os.system('rm -rf '+image_root+'.mask')
    os.system('rm -rf '+image_root+'.pb')
    os.system('rm -rf '+image_root+'.psf')
    os.system('rm -rf '+image_root+'.residual')
    os.system('rm -rf '+image_root+'.weight')
    os.system('rm -rf '+image_root+'.sumwt')

def save_copy_of_imaging(
    input_root=None,
    output_root=None):
    """
    Copy a cube to a new name. Used to make a backup copy. Overwrites
    the previous cube of that name.
    """
    
    wipe_cube(output_root)
    
    logger.debug('save_copy_of_imaging')
    logger.debug('cp -r '+input_root+'.image '+output_root+'.image')
    logger.debug('cp -r '+input_root+'.model '+output_root+'.model')
    logger.debug('cp -r '+input_root+'.mask '+output_root+'.mask')
    logger.debug('cp -r '+input_root+'.pb '+output_root+'.pb')
    logger.debug('cp -r '+input_root+'.psf '+output_root+'.psf')
    logger.debug('cp -r '+input_root+'.residual '+output_root+'.residual')
    logger.debug('cp -r '+input_root+'.weight '+output_root+'.weight')
    logger.debug('cp -r '+input_root+'.sumwt '+output_root+'.sumwt')
    
    os.system('cp -r '+input_root+'.image '+output_root+'.image')
    os.system('cp -r '+input_root+'.model '+output_root+'.model')
    os.system('cp -r '+input_root+'.mask '+output_root+'.mask')
    os.system('cp -r '+input_root+'.pb '+output_root+'.pb')
    os.system('cp -r '+input_root+'.psf '+output_root+'.psf')
    os.system('cp -r '+input_root+'.residual '+output_root+'.residual')
    os.system('cp -r '+input_root+'.weight '+output_root+'.weight')
    os.system('cp -r '+input_root+'.sumwt '+output_root+'.sumwt')

def replace_imaging_with_copy(
    to_root=None,
    from_root=None):
    """
    Replace a cube with a copy.
    """
    
    wipe_cube(to_root)
    
    logger.debug('replace_imaging_with_copy')
    logger.debug('cp -r '+from_root+'.image '+to_root+'.image')
    logger.debug('cp -r '+from_root+'.model '+to_root+'.model')
    logger.debug('cp -r '+from_root+'.mask '+to_root+'.mask')
    logger.debug('cp -r '+from_root+'.pb '+to_root+'.pb')
    logger.debug('cp -r '+from_root+'.psf '+to_root+'.psf')
    logger.debug('cp -r '+from_root+'.residual '+to_root+'.residual')
    logger.debug('cp -r '+from_root+'.weight '+to_root+'.weight')
    logger.debug('cp -r '+from_root+'.sumwt '+to_root+'.sumwt')

    os.system('cp -r '+from_root+'.image '+to_root+'.image')
    os.system('cp -r '+from_root+'.model '+to_root+'.model')
    os.system('cp -r '+from_root+'.mask '+to_root+'.mask')
    os.system('cp -r '+from_root+'.pb '+to_root+'.pb')
    os.system('cp -r '+from_root+'.psf '+to_root+'.psf')
    os.system('cp -r '+from_root+'.residual '+to_root+'.residual')
    os.system('cp -r '+from_root+'.weight '+to_root+'.weight')
    os.system('cp -r '+from_root+'.sumwt '+to_root+'.sumwt')

def export_imaging_to_fits(
    image_root=None,
    bitpix=-32):
    """
    Export the various products associated with a CASA cube to FITS.
    """
    
    logger.debug('export_imaging_to_fits')
    logger.debug('exportfits '+image_root+'.image '+image_root+'.fits')
    exportfits(imagename=image_root+'.image',
               fitsimage=image_root+'.fits',
               velocity=True, overwrite=True, dropstokes=True, 
               dropdeg=True, bitpix=bitpix)
    
    logger.debug('exportfits '+image_root+'.model '+image_root+'_model.fits')
    exportfits(imagename=image_root+'.model',
               fitsimage=image_root+'_model.fits',
               velocity=True, overwrite=True, dropstokes=True, 
               dropdeg=True, bitpix=bitpix)
    
    logger.debug('exportfits '+image_root+'.residual '+image_root+'_residual.fits')
    exportfits(imagename=image_root+'.residual',
               fitsimage=image_root+'_residual.fits',
               velocity=True, overwrite=True, dropstokes=True, 
               dropdeg=True, bitpix=bitpix)
    
    logger.debug('exportfits '+image_root+'.mask '+image_root+'_mask.fits')
    exportfits(imagename=image_root+'.mask',
               fitsimage=image_root+'_mask.fits',
               velocity=True, overwrite=True, dropstokes=True, 
               dropdeg=True, bitpix=bitpix)
    
    logger.debug('exportfits '+image_root+'.pb '+image_root+'_pb.fits')
    exportfits(imagename=image_root+'.pb',
               fitsimage=image_root+'_pb.fits',
               velocity=True, overwrite=True, dropstokes=True, 
               dropdeg=True, bitpix=bitpix)

    return


#endregion







# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# Routines to characterize and manipulate cubes
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

# 
# stat_clean_cube is now in casaMaskingRoutines.py
# 
#def stat_clean_cube(cube_file=None):
#    """
#    Calculate statistics for an image cube.
#    """
#    if cube_file == None:
#        logger.info("No cube file specified. Returning")
#        return
#    imstat_dict = imstat(cube_file)
#    
#    return imstat_dict


def save_copy_of_cube(
    input_root=None,
    output_root=None):
    """
    Copy a cube to a new name. Used to make a backup copy. Overwrites
    the previous cube of that name.
    """

    wipe_cube(output_root)
    
    logger.debug('save_copy_of_cube')
    logger.debug('cp -r '+input_root+'.image '+output_root+'.image')
    logger.debug('cp -r '+input_root+'.model '+output_root+'.model')
    logger.debug('cp -r '+input_root+'.mask '+output_root+'.mask')
    logger.debug('cp -r '+input_root+'.pb '+output_root+'.pb')
    logger.debug('cp -r '+input_root+'.psf '+output_root+'.psf')
    logger.debug('cp -r '+input_root+'.residual '+output_root+'.residual')
    logger.debug('cp -r '+input_root+'.psf '+output_root+'.weight')
    logger.debug('cp -r '+input_root+'.residual '+output_root+'.sumwt')
    
    os.system('cp -r '+input_root+'.image '+output_root+'.image')
    os.system('cp -r '+input_root+'.model '+output_root+'.model')
    os.system('cp -r '+input_root+'.mask '+output_root+'.mask')
    os.system('cp -r '+input_root+'.pb '+output_root+'.pb')
    os.system('cp -r '+input_root+'.psf '+output_root+'.psf')
    os.system('cp -r '+input_root+'.residual '+output_root+'.residual')
    os.system('cp -r '+input_root+'.psf '+output_root+'.weight')
    os.system('cp -r '+input_root+'.residual '+output_root+'.sumwt')


def wipe_cube(
    cube_root=None):
    """
    Wipe files associated with a cube.
    """
    if cube_root == None:
        return
    
    logger.debug('wipe_cube')
    logger.debug('rm -rf '+cube_root+'.image')
    logger.debug('rm -rf '+cube_root+'.model')
    logger.debug('rm -rf '+cube_root+'.mask')
    logger.debug('rm -rf '+cube_root+'.pb')
    logger.debug('rm -rf '+cube_root+'.psf')
    logger.debug('rm -rf '+cube_root+'.residual')
    logger.debug('rm -rf '+cube_root+'.weight')
    logger.debug('rm -rf '+cube_root+'.sumwt')
    
    os.system('rm -rf '+cube_root+'.image')
    os.system('rm -rf '+cube_root+'.model')
    os.system('rm -rf '+cube_root+'.mask')
    os.system('rm -rf '+cube_root+'.pb')
    os.system('rm -rf '+cube_root+'.psf')
    os.system('rm -rf '+cube_root+'.residual')
    os.system('rm -rf '+cube_root+'.weight')
    os.system('rm -rf '+cube_root+'.sumwt')


def replace_cube_with_copy(
    to_root=None,
    from_root=None):
    """
    Replace a cube with a copy.
    """
    
    wipe_cube(to_root)
    
    logger.debug('replace_cube_with_copy')
    logger.debug('cp -r '+from_root+'.image '+to_root+'.image')
    logger.debug('cp -r '+from_root+'.model '+to_root+'.model')
    logger.debug('cp -r '+from_root+'.mask '+to_root+'.mask')
    logger.debug('cp -r '+from_root+'.pb '+to_root+'.pb')
    logger.debug('cp -r '+from_root+'.psf '+to_root+'.psf')
    logger.debug('cp -r '+from_root+'.residual '+to_root+'.residual')
    logger.debug('cp -r '+from_root+'.psf '+to_root+'.weight')
    logger.debug('cp -r '+from_root+'.residual '+to_root+'.sumwt')
    
    os.system('cp -r '+from_root+'.image '+to_root+'.image')
    os.system('cp -r '+from_root+'.model '+to_root+'.model')
    os.system('cp -r '+from_root+'.mask '+to_root+'.mask')
    os.system('cp -r '+from_root+'.pb '+to_root+'.pb')
    os.system('cp -r '+from_root+'.psf '+to_root+'.psf')
    os.system('cp -r '+from_root+'.residual '+to_root+'.residual')
    os.system('cp -r '+from_root+'.psf '+to_root+'.weight')
    os.system('cp -r '+from_root+'.residual '+to_root+'.sumwt')


def export_to_fits(
    cube_root=None,
    bitpix=-32):
    """
    Export the various products associated with a CASA cube to FITS.
    """
    
    logger.debug('export_to_fits')
    logger.debug('exportfits '+cube_root+'.image '+cube_root+'.fits')
    casaStuff.exportfits(\
               imagename=cube_root+'.image',
               fitsimage=cube_root+'.fits',
               velocity=True, overwrite=True, dropstokes=True, 
               dropdeg=True, bitpix=bitpix)
    
    logger.debug('exportfits '+cube_root+'.model '+cube_root+'_model.fits')
    casaStuff.exportfits(\
               imagename=cube_root+'.model',
               fitsimage=cube_root+'_model.fits',
               velocity=True, overwrite=True, dropstokes=True, 
               dropdeg=True, bitpix=bitpix)
    
    logger.debug('exportfits '+cube_root+'.residual '+cube_root+'_residual.fits')
    casaStuff.exportfits(\
               imagename=cube_root+'.residual',
               fitsimage=cube_root+'_residual.fits',
               velocity=True, overwrite=True, dropstokes=True, 
               dropdeg=True, bitpix=bitpix)
    
    logger.debug('exportfits '+cube_root+'.mask '+cube_root+'_mask.fits')
    casaStuff.exportfits(\
               imagename=cube_root+'.mask',
               fitsimage=cube_root+'_mask.fits',
               velocity=True, overwrite=True, dropstokes=True, 
               dropdeg=True, bitpix=bitpix)
    
    logger.debug('exportfits '+cube_root+'.pb '+cube_root+'_pb.fits')
    casaStuff.exportfits(\
               imagename=cube_root+'.pb',
               fitsimage=cube_root+'_pb.fits',
               velocity=True, overwrite=True, dropstokes=True, 
               dropdeg=True, bitpix=bitpix)
    
    return










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
    casaStuff.tclean(**clean_call.kwargs())

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








