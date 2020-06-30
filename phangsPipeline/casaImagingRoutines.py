"""
Standalone routines related to CASA imaging.
"""

#region Imports and definitions

import os
import numpy as np
import pyfits # CASA has pyfits, not astropy
import glob, copy, inspect

import logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

# Analysis utilities
import analysisUtils as au

# CASA stuff
import casaStuff
import casaMaskingRoutines as cmr

# Clean call object
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
    force_square=False,
    ):
    """
    Pick a cell and image size for a measurement set. Requests an
    oversampling factor, which is by default 5. Will pick a good size
    for the FFT and will try to pick a round number for the cell
    size. Returns variables appropriate for cell= and imsize= in
    tclean.
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
        au.pickCellSize(infile, imsize=True,
                        npix=oversamp, intent='')
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

    if force_square == True:
        if cells_y < cells_x:
            cells_y = cells_x
        if cells_x < cells_y:
            cells_x = cells_y

    image_size = [int(cells_x), int(cells_y)]
    cell_size_string = str(cell_size)+'arcsec'

    return cell_size_string, image_size

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

    logger.debug('wipe_imaging under "'+os.getcwd()+'"')
    cmd_list = [
        'rm -rf '+image_root+'.image',
        'rm -rf '+image_root+'.model',
        'rm -rf '+image_root+'.mask',
        'rm -rf '+image_root+'.pb',
        'rm -rf '+image_root+'.psf',
        'rm -rf '+image_root+'.residual',
        'rm -rf '+image_root+'.weight',
        'rm -rf '+image_root+'.sumwt',
        'rm -rf '+image_root+'.alpha',
        'rm -rf '+image_root+'.alpha.error',
        'rm -rf '+image_root+'.beta',
        'rm -rf '+image_root+'.beta.error',
        'rm -rf '+image_root+'.image.tt0',
        'rm -rf '+image_root+'.image.tt1',
        'rm -rf '+image_root+'.image.tt2',
        'rm -rf '+image_root+'.model.tt0',
        'rm -rf '+image_root+'.model.tt1',
        'rm -rf '+image_root+'.model.tt2',
        'rm -rf '+image_root+'.mask.tt0',
        'rm -rf '+image_root+'.mask.tt1',
        'rm -rf '+image_root+'.mask.tt2',
        'rm -rf '+image_root+'.pb.tt0',
        'rm -rf '+image_root+'.pb.tt1',
        'rm -rf '+image_root+'.pb.tt2',
        'rm -rf '+image_root+'.psf.tt0',
        'rm -rf '+image_root+'.psf.tt1',
        'rm -rf '+image_root+'.psf.tt2',
        'rm -rf '+image_root+'.residual.tt0',
        'rm -rf '+image_root+'.residual.tt1',
        'rm -rf '+image_root+'.residual.tt2',
        'rm -rf '+image_root+'.weight.tt0',
        'rm -rf '+image_root+'.weight.tt1',
        'rm -rf '+image_root+'.weight.tt2',
        'rm -rf '+image_root+'.sumwt.tt0',
        'rm -rf '+image_root+'.sumwt.tt1',
        'rm -rf '+image_root+'.sumwt.tt2',
        ]

    for this_cmd in cmd_list:
        logger.debug(this_cmd)
        os.system(this_cmd+' 2>/dev/null')

    return()

def copy_imaging(
    input_root=None,
    output_root=None,
    wipe_first=True):
    """
    Copy all of the files from a cube or continuum imaging output by
    clean to have a new root name. Most commonly used to make a backup
    copy of imaging output during iterative imaging (e.g., clean
    loops, shifting clean modes, selfcal, etc). Overwrites any
    previous imaging with that output name.
    """

    if wipe_first:
        wipe_imaging(output_root)

    logger.debug('Copying imaging from root '+input_root+' to root '+output_root)
    cmd_list = [
        'cp -r '+input_root+'.image '+output_root+'.image',
        'cp -r '+input_root+'.model '+output_root+'.model',
        'cp -r '+input_root+'.residual '+output_root+'.residual',
        'cp -r '+input_root+'.mask '+output_root+'.mask',
        'cp -r '+input_root+'.pb '+output_root+'.pb',
        'cp -r '+input_root+'.psf '+output_root+'.psf',
        'cp -r '+input_root+'.weight '+output_root+'.weight',
        'cp -r '+input_root+'.sumwt '+output_root+'.sumwt',
        'cp -r '+input_root+'.alpha '+output_root+'.alpha',
        'cp -r '+input_root+'.alpha.error '+output_root+'.alpha.error',
        'cp -r '+input_root+'.beta '+output_root+'.beta',
        'cp -r '+input_root+'.beta.error '+output_root+'.beta.error',
        'cp -r '+input_root+'.image.tt0 '+output_root+'.image.tt0',
        'cp -r '+input_root+'.image.tt1 '+output_root+'.image.tt1',
        'cp -r '+input_root+'.image.tt2 '+output_root+'.image.tt2',
        'cp -r '+input_root+'.model.tt0 '+output_root+'.model.tt0',
        'cp -r '+input_root+'.model.tt1 '+output_root+'.model.tt1',
        'cp -r '+input_root+'.model.tt2 '+output_root+'.model.tt2',
        'cp -r '+input_root+'.residual.tt0 '+output_root+'.residual.tt0',
        'cp -r '+input_root+'.residual.tt1 '+output_root+'.residual.tt1',
        'cp -r '+input_root+'.residual.tt2 '+output_root+'.residual.tt2',
        'cp -r '+input_root+'.mask.tt0 '+output_root+'.mask.tt0',
        'cp -r '+input_root+'.mask.tt1 '+output_root+'.mask.tt1',
        'cp -r '+input_root+'.mask.tt2 '+output_root+'.mask.tt2',
        'cp -r '+input_root+'.pb.tt0 '+output_root+'.pb.tt0',
        'cp -r '+input_root+'.pb.tt1 '+output_root+'.pb.tt1',
        'cp -r '+input_root+'.pb.tt2 '+output_root+'.pb.tt2',
        'cp -r '+input_root+'.psf.tt0 '+output_root+'.psf.tt0',
        'cp -r '+input_root+'.psf.tt1 '+output_root+'.psf.tt1',
        'cp -r '+input_root+'.psf.tt2 '+output_root+'.psf.tt2',
        'cp -r '+input_root+'.weight.tt0 '+output_root+'.weight.tt0',
        'cp -r '+input_root+'.weight.tt1 '+output_root+'.weight.tt1',
        'cp -r '+input_root+'.weight.tt2 '+output_root+'.weight.tt2',
        'cp -r '+input_root+'.sumwt.tt0 '+output_root+'.sumwt.tt0',
        'cp -r '+input_root+'.sumwt.tt1 '+output_root+'.sumwt.tt1',
        'cp -r '+input_root+'.sumwt.tt2 '+output_root+'.sumwt.tt2',
        ]

    for this_cmd in cmd_list:
        logger.debug(this_cmd)
        os.system(this_cmd+' 2>/dev/null')

def export_imaging_to_fits(
    image_root=None,
    bitpix=-32,
    just_image=False):
    """
    Export the products associated with a CASA imaging run to FITS.
    """

    ext_map = {
        '.alpha':'_alpha.fits',
        '.alpha.error':'_alpha_error.fits',
        '.beta':'_beta.fits',
        '.beta.error':'_beta_error.fits',
        '.image.tt0':'.fits',
        '.image.tt1':'_tt1.fits',
        '.image.tt2':'_tt2.fits',
        '.model.tt0':'_model.fits',
        '.model.tt1':'_model_tt1.fits',
        '.model.tt2':'_model_tt2.fits',
        '.residual.tt0':'_residual.fits',
        '.residual.tt1':'_residual_tt1.fits',
        '.residual.tt2':'_residual_tt2.fits',
        '.mask.tt0':'_mask.fits',
        '.mask.tt1':'_mask_tt1.fits',
        '.mask.tt2':'_mask_tt2.fits',
        '.pb.tt0':'_pb.fits',
        '.pb.tt1':'_pb_tt1.fits',
        '.pb.tt2':'_pb_tt2.fits',
        '.psf.tt0':'_psf.fits',
        '.psf.tt1':'_psf_tt1.fits',
        '.psf.tt2':'_psf_tt2.fits',
        '.weight.tt0':'_weight.fits',
        '.weight.tt1':'_weight_tt1.fits',
        '.weight.tt2':'_weight_tt2.fits',
        '.image':'.fits',
        '.model':'_model.fits',
        '.residual':'_residual.fits',
        '.mask':'_mask.fits',
        '.pb':'_pb.fits',
        '.psf':'_psf.fits',
        '.weight':'_weight.fits',
        }

    for this_ext in ext_map.keys():
        if just_image and ((this_ext != '.tt0') and this_ext != '.image'):
            continue

        this_casa_ext = this_ext
        this_fits_ext = ext_map[this_ext]

        casa_image = image_root + this_casa_ext
        if not os.path.isdir(casa_image):
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
# Execute a clean call
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

#region clean call execution

def execute_clean_call(
    clean_call = None,
    reset = False,
    ):
    """
    Execute a clean call object, optionally deleting previous versions
    of the imaging first.

    Set reset to wipe any existing output data.
    """

    if not isinstance(clean_call, CleanCall):
        logger.error("Please input a valid clean call!")
        raise Exception("Please input a valid clean call!")

    if not clean_call.has_param('vis'):
        logger.error("No visibility defined in clean_call. Returning.")
        return

    if not clean_call.has_param('imagename'):
        logger.error("No imagename defined in clean_call. Returning.")
        return

    #<TODO><DEBUG><DL># this will not overwrite existing data
    #if os.path.isdir(clean_call.get_param('imagename')+'.image') and not os.path.isdir(clean_call.get_param('imagename')+'.image'+'.touch'):
    #    logger.info('Found existing data "'+clean_call.get_param('imagename')+'.image'+'". Will not overwrite.')
    #    return

    if not os.path.isdir(clean_call.get_param('vis')):
        logger.error("Visibility file not found: "+clean_call.get_param('vis'))
        return

    if clean_call.logfile != None:
        oldlogfile = casaStuff.casalog.logfile()
        casaStuff.casalog.setlogfile(clean_call.logfile)

    if reset:
        logger.debug("Wiping previous versions of the cube.")
        wipe_imaging(clean_call.get_param('imagename'))


    # a simple way to slightly solve the compatible issue is to check
    # the list of expected_kwargs and only return keys inside it.
    logger.debug("Running CASA "+str(clean_call))
    clean_kwargs = clean_call.kwargs_for_clean()
    expected_kwargs = inspect.getargspec(casaStuff.tclean)[0]
    active_kwargs = {} # kwarg dict
    if expected_kwargs is not None:
        missing_kwargs = [] # list
        unused_kwargs = [] # list
        for k in expected_kwargs:
            if k in clean_kwargs:
                active_kwargs[k] = clean_kwargs[k]
            else:
                missing_kwargs.append(k)
        for k in clean_kwargs:
            if not (k in expected_kwargs):
                unused_kwargs.append(k)
        if len(unused_kwargs) > 0:
            logger.warning('Unused key arguments for clean: '+str(unused_kwargs))
        if len(missing_kwargs) > 0:
            logger.warning('Missing key arguments for clean: '+str(missing_kwargs)+'. Caution that CASA will use some default values depending on the CASA version.')
    else:
        active_kwargs = copy.deepcopy(clean_kwargs)

    # Force pbmask == pblimit
    # For cubes with an empty channel, pbmask < pblimit triggers
    # NaNs in cleaning as of CASA 5.6.1.
    # Force pbmask >= pblimit to avoid this.
    if 'pblimit' in active_kwargs and 'pbmask' in active_kwargs:
        if active_kwargs['pblimit'] > active_kwargs['pbmask']:
            active_kwargs['pbmask'] = active_kwargs['pblimit']

    #os.mkdir(clean_call.get_param('imagename')+'.image'+'.touch') #<TODO><DEBUG><DL>#
    casaStuff.tclean(**active_kwargs)
    #os.rmdir(clean_call.get_param('imagename')+'.image'+'.touch') #<TODO><DEBUG><DL>#

    if clean_call.logfile != None:
        casaStuff.casalog.setlogfile(oldlogfile)

    return

#endregion

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# Run a clean call with NITER=0 to make a dirty image
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

def make_dirty_image(
    clean_call = None,
    ):
    """
    Create a dirty image using the provided clean call. Forces number
    of iterations to zero before excuting the clean call and enforces
    psf and residual caculation but otherwise leaves the clean_call
    unchanged. Making a dirty image also forces a reset, wiping any
    previous version of the imaging. Avoids mutating the clean_call.
    """

    if not isinstance(clean_call, CleanCall):
        logger.error("Please input a valid clean call!")
        raise Exception("Please input a valid clean call!")

    dirty_clean_call = copy.deepcopy(clean_call)

    dirty_clean_call.set_param('niter', 0)
    dirty_clean_call.set_param('calcres',True)
    dirty_clean_call.set_param('calcpsf',True)

    execute_clean_call(dirty_clean_call, reset=True)

    return()


# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# Repeated run a clean call, looking for convergence
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

def eval_niter(
    loopnum=1,
    baseval=10,
    model='geometric',
    factor=2.0,
    saturation=1000,
    other_input=None,
    ):
    """
    Helper function to evaluate the number of iterations.
    """
    niter = None

    # Fixed number of iterations

    if model.lower() == 'fixed':
        niter = baseval

    # A geometric model starts at the base value and scales by
    # factor each time.

    if model.lower() == 'geometric':
        niter = baseval*factor**(loopnum)

    # A linear model starts at the base value and increases by
    # baseval*(factor*loopnum) each time.

    if model.lower() == 'linear':
        niter = baseval*(1.0+factor*loopnum)

    # Experimental/untested: a sequence of iterations
    if model.lower() == 'sequence':
        if loopnum >= len(other_input):
            index = len(other_input)-1
        else:
            index = loopnum
        niter = other_input[loopnum]

    # Experimental/untested: an expression to be 'exec'ed
    if model.lower() == 'expr':
        exec(other_input)

    # Cap at the saturation value
    if saturation is not None:
        if niter >= saturation:
            niter = saturation

    return int(niter)

def clean_loop(
    clean_call=None,
    record_file=None,
    suffix='',
    log_ext=None,
    niter_base_perchan = 10,
    niter_growth_model = 'geometric',
    niter_growth_factor = 2.0,
    niter_saturation_perchan = 1000,
    niter_other_input = None,
    cycleniter_base = 100,
    cycleniter_growth_model='linear',
    cycleniter_growth_factor = 1.0,
    cycleniter_saturation_value = 1000,
    cycleniter_other_input = None,
    threshold_type = 'snr',
    threshold_value = 4.0,
    min_loops = 0,
    max_loops = 20,
    max_total_niter = None,
    convergence_fracflux=0.02,
    convergence_totalflux=None,
    convergence_fluxperniter=None,
    use_absolute_delta=True,
    stop_at_negative=True,
    remask_each_loop=False,
    force_dirty_image=False,
    ):
    """
    Carry out an iterative clean until a convergence criteria is
    met. The loop releases progressively more iterations to the
    provided clean call and checks for convergence after each call to
    clean.

    Optional parameters control the growth of major cycle
    (cycle_niter) and total (niter) iterations/components, the
    implementation of absolute or signal-to-noise thresholds, and the
    total number of loops allowed.

    Convergence is checked by looking at the fractional change in the
    integrated model flux during successive iterations. With a strong
    threshold and enough iterations released, this can practically
    reduce to approaching the signal to noise threshold. The user also
    supplies a maximum number of loops, which prevent running away in
    cases of divergence.

    Future improvements could allow the execution of arbitrary
    functions to characterize the data and check for convergence or to
    manipulate the clean mask. Right now "remasking" and "stat_cube"
    are the only option implemented in this step.
    """

    # Some definitions and error checking

    if not isinstance(clean_call, CleanCall):
        logger.error("Please input a valid clean call!")
        raise Exception("Please input a valid clean call!")

    if suffix == '':
        if clean_call.has_param('specmode'):
            if clean_call.get_param('specmode') == 'mfs':
                suffix = '.tt0'

    valid_model_types = ['fixed', 'geometric','linear','sequence','expr']
    for growth_model in [niter_growth_model.lower(), cycleniter_growth_model.lower()]:
        if growth_model not in valid_model_types:
            logger.warning("Growth model not recognized: ", growth_model)
            return()

    valid_threshold_types = ['snr','absolute']
    if threshold_type.lower() not in valid_threshold_types:
        logger.warning("Threshold type not recognized: ", threshold_type)
        return()

    # Check if a residual image exists. If not, then build the dirty
    # image. Also build the dirty image if the flag to
    # force_dirt_image is set to True.

    missing_image = True
    if os.path.isdir(clean_call.get_param('imagename')+'.residual'+suffix):
        missing_image = False

    if missing_image or force_dirty_image:
        make_dirty_image(clean_call)

    # Copy the clean call so we can manipulate it without changing the
    # input version call.

    working_call = copy.deepcopy(clean_call)
    working_call.set_param('calcres',False)
    working_call.set_param('calcpsf',False)

    # Note the number of channels, which is used in setting the number
    # of iterations that we give to an individual clean call.

    vm = au.ValueMapping(working_call.get_param('vis'))
    nchan = vm.spwInfo[0]['numChannels']

    # Create a text record of progress through successive clean calls.

    record = []
    record.append("loopnum, deconvolver, niter, cycleniter, threshold, noise, model_flux, frac_delta_flux\n")
    record.append("# column 1: Loop number.\n")
    record.append("# column 2: Deconvolver used in clean.\n")
    record.append("# column 3: Allocated number of iterations.\n")
    record.append("# column 4: Cycleniter used to force major cycles.\n")
    record.append("# column 5: Threshold supplied to clean.\n")
    record.append("# column 6: Noise level measured in residuals.\n")
    record.append("# column 7: Integrated model flux.\n")
    record.append("# column 8: Fractional change in flux from previous loop.\n")

    # Initialize the loop counter and our tracking of the flux in the
    # model (which we use to estimate convergence).

    loop = 0
    cumulative_niter = 0
    previous_flux = 0.0
    current_flux = 0.0

    # Run the main loop

    proceed = True
    while proceed == True:

        # Calculate the number of total iterations (niter) and
        # iterations per major cycle (cycleniter) released to clean
        # during this call.

        niter = eval_niter(loopnum = loop,
                           baseval = niter_base_perchan*nchan,
                           model=niter_growth_model, factor=niter_growth_factor,
                           saturation = niter_saturation_perchan*nchan,
                           other_input=niter_other_input)

        cycleniter = eval_niter(loopnum = loop, baseval = cycleniter_base,
                                model=cycleniter_growth_model, factor=cycleniter_growth_factor,
                                saturation = cycleniter_saturation_value,
                                other_input=cycleniter_other_input)

        working_call.set_param('niter', niter, nowarning=True)
        working_call.set_param('cycleniter', cycleniter, nowarning=True)

        cumulative_niter = cumulative_niter + niter

        # Calculate the current noise in the residual image. Don't
        # exclude the masked region from the noise calculation but do
        # turn on iterative noise estimation (Chauvenet+m.a.d. using 5
        # iterations should be quite robust).

        current_noise = cmr.noise_for_cube(
            infile=working_call.get_param('imagename')+'.residual'+suffix,
            method='chauvmad', niter=5)

        # Set the threshold for the clean call. Clean expects a value
        # in Jy/beam. Switch on the threshold type to make the string
        # and attach it to the clean call.

        if threshold_type == 'snr':
            threshold_string = str(current_noise*threshold_value)+'Jy/beam'
        elif threshold_type == 'absolute':
            if type(threshold_value) == type(0.0):
                threshold_string = str(threshold_value)+'Jy/beam'
            else:
                threshold_string = threshold_value
        else:
            threshold_string = '0.0Jy/beam'

        working_call.set_param('threshold', threshold_string, nowarning=True)

        # If requested mask at each step (this is experimental, we're
        # seeing if it helps to avoid divergence during the deep
        # single scale clean.)

        if remask_each_loop:
            logger.info("")
            logger.info("Remasking.")
            logger.info("")
            signal_mask(
                cube_root=working_call.get_param('imagename'),
                out_file=working_call.get_param('imagename')+'.mask'+suffix,
                suffix_in=suffix,
                suffix_out=suffix,
                operation='AND',
                high_snr=4.0,
                low_snr=2.0,
                absolute=False)
            working_call.usemask='user'

        # Set the log file (revisit this)

        if log_ext is not None:
            working_call.logfile = working_call.get_param('imagename')+"_loop_"+str(loop)+"_"+log_ext+".log"
        else:
            working_call.logfile = None

        # Save the previous version of the imaging for comparison

        copy_imaging(
            input_root=working_call.get_param('imagename'),
            output_root=working_call.get_param('imagename')+'_prev')

        # Execute the clean call.

        execute_clean_call(working_call)

        # Calculate the new model flux and the change relative to the
        # previous step, normalized by current flux and by iterations.

        model_stats = cmr.stat_cube(working_call.get_param('imagename')+'.model'+suffix)

        previous_flux = current_flux
        current_flux = model_stats['sum'][0]

        delta_flux = (current_flux-previous_flux)
        if use_absolute_delta:
            delta_flux = abs(delta_flux)

        flux_per_iter = delta_flux / niter
        frac_delta_flux = delta_flux / previous_flux

        # Check whether the model flux convergence criteria is met

        if convergence_fracflux is not None:
            if frac_delta_flux < convergence_fracflux:
                proceed = False

        if convergence_totalflux is not None:
            if delta_flux < convergence_totalflux:
                proceed = False

        if convergence_fluxperniter is not None:
            if flux_per_iter < convergence_fluxperniter:
                proceed = False

        if max_total_niter is not None:
            if cumulative_niter >= max_total_niter:
                proceed = False

        # If requested, stop if the integrated model flux becomes
        # negative.

        if stop_at_negative:
            if current_flux < 0.0:
                proceed = False

        # Enforce minimum and maximum limits on number of loops. These
        # override other convergence criteria.

        if loop >= max_loops:
            proceed = False
        if loop < min_loops:
            proceed = True

        # Generate a record line and print the current status to the screen

        this_record = ''
        this_record += str(loop)+', '
        this_record += str(working_call.get_param('deconvolver'))+', '
        this_record += str(working_call.get_param('niter'))+', '
        this_record += str(working_call.get_param('cycleniter'))+', '
        this_record += str(working_call.get_param('threshold'))+', '
        this_record += str(current_noise)+'Jy/beam, '
        this_record += str(current_flux)+'Jy*chan, '
        this_record += str(frac_delta_flux)+''
        this_record += '\n'

        # Print the current record to the screen

        record.append(this_record)
        for line in record:
            print(line)

        logger.info("... proceeding? "+str(proceed))

        if proceed == False:
            break

        loop += 1

    # ... if requested also write this to a file.

    if record_file != None:
        f = open(record_file,'w')
        f.writelines(record)
        f.close()

    return()

#endregion

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# Evaluate the output of imaging
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

def calc_residual_statistics(
    resid_name=None,
    mask_name=None,
    ):
    """
    """
    
    if os.path.isdir(resid_name) == False:
        logger.error('Error! The input file "'+resid_name+'" was not found!')
        return

    if os.path.isdir(mask_name) == False:
        logger.error('Error! The input file "'+mask_name+'" was not found!')
        return
    
    myia = au.createCasaTool(casaStuff.iatool)

    myia.open(mask_name)
    mask = myia.getchunk()    
    myia.close()

    myia.open(resid_name)
    resid = myia.getchunk()    
    myia.close()

    vec = resid[((mask == 1)*np.isfinite(resid))]
    del mask
    del resid

    current_noise = cmr.noise_for_cube(
        infile=resid_name,
        method='chauvmad', niter=5)
    
    out_dict = {
        'cubename':resid_name,
        'maskname':mask_name,
        'max':np.max(vec),
        'p99':np.percentile(vec,99),
        'p95':np.percentile(vec,95),
        'p90':np.percentile(vec,90),
        'noise':current_noise,
        }
    
    return(out_dict)
    
    
