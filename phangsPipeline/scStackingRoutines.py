import numpy as np
import scipy.ndimage as nd
import astropy.units as u
import astropy.wcs as wcs
from astropy.io import fits
from spectral_cube import SpectralCube
from reproject import reproject_interp
import logging

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

def channelShiftVec(x, ChanShift):
    # Shift an array of spectra (x) by a set number of Channels (array)
    ftx = np.fft.fft(x, axis=0)
    m = np.fft.fftfreq(x.shape[0])
    phase = np.exp(-2 * np.pi * m[:, np.newaxis]
                   * 1j * ChanShift[np.newaxis, :])
    x2 = np.real(np.fft.ifft(ftx * phase, axis=0))
    return(x2)

def ShuffleCube(DataCube, vfield, chunk=1000):
    """
    Shuffles cube so that the velocity appearing in the vfield is set 
    to the middle channel and velocity vfield.
    
    Parameters
    ----------
    DataCube : SpectralCube
        The original spectral cube with spatial dimensions Nx, Ny and spectral dimension Nv
    vfield : 2D numpy.ndarray
        A 2D map of the vfield velocities for the lines to stack of dimensions Nx, Ny.
        Note that DataCube and vfield map must have equivalent (but not necessarily equal) 
        spectral units (e.g., km/s and m/s)
    
    Keywords
    --------
    chunk : int
        Number of data points to include in a chunk for processing.
    
    Returns
    -------
    OutCube : np.array
        Output SpectralCube cube shuffled so that the emission is at 0.
 
    """

    spaxis = DataCube.spectral_axis
    y, x = np.where(np.isfinite(vfield))
    vfield = vfield.to(spaxis.unit)
    # v0 = spaxis[len(spaxis) // 2]
    # newspaxis = spaxis - v0
    relative_channel = np.arange(len(spaxis)) - (len(spaxis) // 2)
    vfields = vfield[y, x]
    sortindex = np.argsort(spaxis)
    channel_shift = -1 * np.interp(vfields, spaxis[sortindex],
                                   np.array(relative_channel[sortindex], dtype=float))
    NewCube = np.empty(DataCube.shape)
    NewCube.fill(np.nan)
    nchunk = (len(x) // chunk)
    for thisx, thisy, thisshift in zip(np.array_split(x, nchunk), 
                                       np.array_split(y, nchunk),
                                       np.array_split(channel_shift, nchunk)):
        spectrum = DataCube.filled_data[:, thisy, thisx].value
        baddata = ~np.isfinite(spectrum)
        shifted_spectrum = channelShiftVec(np.nan_to_num(spectrum),
                                           np.atleast_1d(thisshift))
        # shifted_mask = channelShiftVec(baddata, np.atleast_1d(thisshift))
        shifted_spectrum[baddata>0] = np.nan
        NewCube[:, thisy, thisx] = shifted_spectrum
    hdr = DataCube.header
    hdr['CRVAL3'] = 0.0
    hdr['CRPIX3'] = len(spaxis) // 2 + 1
    newwcs = wcs.WCS(hdr)
    NewCube *= u.Unit(hdr['BUNIT'])
    return(SpectralCube(NewCube, newwcs, header=hdr))


def recipe_phangs_vfield(
    template_vfield,
    list_of_vfields=None,
    outfile=None,
    overwrite=False):
    """
    Task to create combined velocity field from a list of velocity maps in
    a hierarchical manner (from first to last).

    Parameters:

    -----------

    template_vfield : string or fits.hdu

        The reference velocity filed that holds the target WCS. The 
        other maps will be reprojected onto this one. This map is 
        included in the final output.

    Keywords:
    ---------

    list_of_vfields : list or list of fits.hdu

        List of velocity fields or fits.hdus. These will be reprojected
        onto the template vfield and then combined via hierarchical addition.

    outfile : string
        Filename where the vfield will be written. The vfield maps is alsoreturned.

    """

    # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
    # Error checking and work out inputs
    # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

    # check input vfield map
    if (type(template_vfield) == fits.hdu.hdulist.HDUList):
        if (type(template_vfield.data) is np.array) & (np.ndim(template_vfield.data) == 2):
            hdr_template = template_vfield.header
            template_vfield = template_vfield.data * u.Unit(hdr_template['BUNIT'])
        else:
            logger.error("Input velocity vfield map must be a 2D np.array of dimensions Nx, Ny.")
    elif type(template_vfield) == type("hello"):
        template_vfield, hdr_template = fits.getdata(template_vfield, header=True)
        template_vfield *= u.Unit(hdr_template['BUNIT'])
    else:
        logger.error("Input velocity vfield map must be a 2D np.array of dimensions Nx, Ny.")

    # initialise combined velocity field
    vfield_combined = np.copy(template_vfield)
    vfield_combined[np.isnan(vfield_combined)] = 0

    # loop over list of velocity fields
    for this_vfield in list_of_vfields:

        # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
        # Error checking and work out inputs
        # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

        # check input vfield map
        if (type(this_vfield) == fits.hdu.hdulist.HDUList):
            if (type(this_vfield.data) is np.array) & (np.ndim(this_vfield.data) == 2):
                this_vfield = this_vfield.data * u.Unit(this_vfield.header['BUNIT'])
            else:
                logger.error("Input velocity vfield map must be a 2D np.array of dimensions Nx, Ny.")
        elif type(this_vfield) == type("hello"):
            this_vfield, hdr_this_vfield = fits.getdata(this_vfield, header=True)
            this_vfield *= u.Unit(hdr_this_vfield['BUNIT'])
        else:
            logger.error("Input velocity vfield map must be a 2D np.array of dimensions Nx, Ny.")

        # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
        # Reproject velocity field to template grid
        # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

        # regrid velocity field if grids do not match
        if (template_vfield.shape[0] != np.shape(this_vfield)[0]) | (template_vfield.shape[1] != np.shape(this_vfield)[1]):
            logger.warning('Regridding velocity field with interpolated reprojection.')
            this_vfield, _ = reproject_interp((this_vfield, hdr_this_vfield), hdr_template)
            this_vfield *= u.Unit(hdr_this_vfield['BUNIT'])

        # check that grids match after regridding
        if template_vfield.shape[0] != np.shape(this_vfield)[0]:
            logger.error("Input velocity vfield map must be a 2D np.array of dimensions Nx, Ny, but Nx does not match.")
        if template_vfield.shape[1] != np.shape(this_vfield)[1]:
            logger.error("Input velocity vfield map must be a 2D np.array of dimensions Nx, Ny, but Ny does not match.")

        # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
        # Merge velocity fields
        # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

        # fill empty area with current vfield
        mom1_temp = np.copy(this_vfield)
        mom1_temp[vfield_combined != 0] = 0
        vfield_combined += mom1_temp
        vfield_combined[np.isnan(vfield_combined)] = 0

    # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
    # Write or return as requested
    # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

    # prepate velocity field for output
    vfield_combined[vfield_combined == 0] = np.nan

    # Write to disk, if desired
    if outfile is not None:
        fits.writeto(outfile, vfield_combined.value, hdr_template, overwrite=overwrite)

    # return combined velocity field
    return(vfield_combined)

def recipe_shuffle_cube(
    incube,
    invfield,
    outfile=None,
    return_spectral_cube=False,
    overwrite=False):
    """
    Task to create velocity-shuffled cubes via input velocity vfield map.

    Parameters:
    -----------

    incube : string or SpectralCube
        The cube to be masked.

    invfield : 2D numpy.ndarray
        A 2D map of the vfield velocities for the lines to stack of dimensions Nx, Ny.
        Note that DataCube and vfield map must have equivalent (but not necessarily equal) 
        spectral units (e.g., km/s and m/s)
    

    Keywords:
    ---------

    outfile : string
        Filename where the shuffled cube will be written. The shuffled cube is also
        returned.

    """

    # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
    # Error checking and work out inputs
    # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

    # check input cube
    if type(incube) is SpectralCube:
        cube = incube
    elif type(incube) == type("hello"):
        cube = SpectralCube.read(incube)
    else:
        logger.error("Input cube must be a SpectralCube object or a filename.")

    cube.allow_huge_operations = True

    # check input vfield map
    if (type(invfield) is np.array) & (np.ndim(invfield) == 2):
        vfield = invfield
    elif type(invfield) == type("hello"):
        vfield, hdr_map = fits.getdata(invfield, header=True)
        vfield *= u.Unit(hdr_map['BUNIT'])
    else:
        logger.error("Input velocity vfield map must be a 2D np.array of dimensions Nx, Ny.")

    # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
    # Reproject velocity field to cube grid
    # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

    # regrid velocity field if grids do not match
    if (cube.shape[1] != np.shape(vfield)[0]) | (cube.shape[2] != np.shape(vfield)[1]):
        logger.warning('Regridding velocity field with interpolated reprojection.')
        vfield, _ = reproject_interp((vfield, hdr_map), cube.moment0().header)
        vfield *= u.Unit(hdr_map['BUNIT'])

    # check that grids match after regridding
    if cube.shape[1] != np.shape(vfield)[0]:
        logger.error("Input velocity vfield map must be a 2D np.array of dimensions Nx, Ny, but Nx does not match.")
    if cube.shape[2] != np.shape(vfield)[1]:
        logger.error("Input velocity vfield map must be a 2D np.array of dimensions Nx, Ny, but Ny does not match.")

    # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
    # Run the shuffling
    # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

    # run shuffling
    shuffled_cube = ShuffleCube(cube, vfield, chunk=1000)    

    # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
    # Write or return as requested
    # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

    # prepare cube for output
    # convert data type to 32-bit float and assign WCS and header
    shuffled_cube = SpectralCube(shuffled_cube.filled_data[:].astype(np.float32), wcs=shuffled_cube.wcs, header=shuffled_cube.header)

    # Write to disk, if desired
    if outfile is not None:
        shuffled_cube.write(outfile, overwrite=overwrite)

    # return cube
    if return_spectral_cube:
        return(shuffled_cube)
    else:
        return(shuffled_cube.filled_data[:].value)


def BinByMask(DataCube, mask, vfield, weight_map=None):
    """
    Bin a data cube by a label mask, aligning the data to a common vfield.  Returns an array.

    Parameters
    ----------
    DataCube : SpectralCube
        The original spectral cube with spatial dimensions Nx, Ny and spectral dimension Nv
    Mask : 2D numpy.ndarray
        A 2D map containing boolean values with True indicate where the spectra should be aggregated.
    vfield : 2D numpy.ndarray
        A 2D map of the vfield velocities for the lines to stack of dimensions Nx, Ny.
        Note that DataCube and vfield map must have equivalent spectral units (e.g., km/s)
    weight_map : 2D numpy.ndarray
        Map containing the weight values to be used in averaging
    Returns
    -------
    Spectrum : np.array
        Spectrum of average over mask.
    """
    spaxis = DataCube.spectral_axis.value
    y, x = np.where(mask)
    v0 = spaxis[len(spaxis) // 2] * DataCube.spectral_axis.unit
    relative_channel = np.arange(len(spaxis)) - (len(spaxis) // 2)
    vfields = vfield[y, x].to(DataCube.spectral_axis.unit).value
    sortindex = np.argsort(spaxis)
    channel_shift = -1 * np.interp(vfields, spaxis[sortindex],
                                   np.array(relative_channel[sortindex], dtype=float))
    spectra = DataCube.filled_data[:, y, x].value
    shifted_spectra = channelShiftVec(spectra, channel_shift)
    if weight_map is not None:
        wts = weight_map[y, x]
    else:
        wts = np.ones_like(channel_shift)
    accum_spectrum = np.nansum(wts[np.newaxis, :]
                               * shifted_spectra, axis=1) / np.nansum(wts)
    shifted_spaxis = (DataCube.spectral_axis - v0)
    return(accum_spectrum, shifted_spaxis)


def BinByLabel(DataCube, LabelMap, vfield,
               weight_map=None,
               background_labels=[0]):
    """
    Bin a data cube by a label mask, aligning the data to a common vfield.

    Parameters
    ----------
    DataCube : SpectralCube
        The original spectral cube with spatial dimensions Nx, Ny and spectral dimension Nv
    LabelMap : 2D numpy.ndarray
        A 2D map containing integer labels for each pixel into objects defining the stacking.
    vfield : 2D numpy.ndarray
        A 2D map of the vfield velocities for the lines to stack of dimensions Nx, Ny.
        Note that DataCube and vfield map must have equivalent spectral units to DataCube
    background_labels : list
        List of values in the label map that correspond to background objects and should not
        be processed with the stacking. 

    Returns
    -------
    output_list : list of dict
        List of dict where each entry contains the stacked spectrum for a given label
    unique_labels = array of unique labels in same order as output list
    """
    UniqLabels = np.unique(LabelMap)
    output_list = []
    for ThisLabel in UniqLabels:
        if ThisLabel not in background_labels:
            thisspec, spaxis = BinByMask(DataCube,
                                         (LabelMap == ThisLabel),
                                         vfield,
                                         weight_map=weight_map)
            output_list += [{'label': ThisLabel,
                             'spectrum': thisspec,
                             'spectral_axis': spaxis}]
    return(output_list, UniqLabels)
