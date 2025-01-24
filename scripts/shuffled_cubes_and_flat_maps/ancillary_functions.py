######################################################################################################
# PHANGS-ALMA (DR5) shuffled cubes pipeline, ancillary scripts
# Date: Jan 23, 2025
# Authors: Lukas Neumann, lukas.neumann@eso.org
######################################################################################################

import numpy as np
from astropy.io import fits   # to load fits files
from astropy.wcs import WCS  # for coordinates
from matplotlib.path import Path  # create paths


def shuffle_cube(cube, vaxis, vfield):
    
    """
    Takes a cube, its velocity axis and a velocity field and shuffles the cube in 
    integer channels according to the velocity field.
    """
    
    # copy cube
    cube_shuffled = np.copy(cube)
    
    # define reference channel = the zero-velocity channel
    n_ch = len(vaxis)  # number of channels
    ch_med = n_ch//2  # reference channel
    v_med = vaxis[ch_med]  # reference velocity
    ch_axis = np.arange(n_ch) - ch_med
    
    # create cube with velocity values
    vaxis_cube = np.zeros_like(cube, dtype=float)
    for i in range(len(vaxis)):
        vaxis_cube[i,:,:] = vaxis[i]
    
    # resample velocity field to velocity axis
    vfield_ch = np.argmin(np.abs(vfield - vaxis_cube), axis=0)  # channel indeces of velocity field
    vfield_ch_delta = vfield_ch - ch_med  # shift by reference channel

    # loop over velocity shuffling array
    for ch_roll in ch_axis:
        
        # get spaxels where map matches shuffling value
        idx_x = np.where(vfield_ch_delta == ch_roll)[0]
        idx_y = np.where(vfield_ch_delta == ch_roll)[1]
        
        # shuffle cube
        cube_shuffled[:,idx_x, idx_y] = np.roll(cube[:,idx_x, idx_y], -ch_roll, axis=0)

    return cube_shuffled


def get_major_axis_bins(ctr_ra, ctr_dec, posang, header, bin_width_major, bin_width_minor, rgal_axis_length):

    """
    Takes a galaxy's fits file and its coordinates and return a PV diagram array.
    """

    # make 2D header
    del header['*3*']
    header['NAXIS'] = 2
    header['WCSAXES'] = 2

    # get coordinatesspectrum
    wcs = WCS(header)
    
    # get map dimensions
    n_ra = header['NAXIS1']
    n_dec = header['NAXIS2']

    # convert centre to pixel coordinates
    ctr_ra_pix, ctr_dec_pix = wcs.wcs_world2pix([[ctr_ra, ctr_dec]], 1)[0]  

    # convert axis length from degree to pixel coordinates
    _, axis_max_pix = wcs.wcs_world2pix([[ctr_ra, ctr_dec]], 1)[0]  
    _, axis_ctr_pix = wcs.wcs_world2pix([[ctr_ra, ctr_dec + rgal_axis_length]], 1)[0]  
    rgal_axis_length_pix = np.abs(axis_max_pix - axis_ctr_pix)

    # make axis oriented with north-south axis
    dec_min = wcs.wcs_pix2world([[ctr_ra_pix, ctr_dec_pix - rgal_axis_length_pix]], 1)[0][1]
    dec_max = wcs.wcs_pix2world([[ctr_ra_pix, ctr_dec_pix + rgal_axis_length_pix]], 1)[0][1]
    major_dec = np.concatenate((np.arange(ctr_dec-bin_width_major/2, dec_min, -bin_width_major), np.arange(ctr_dec+bin_width_major/2, dec_max, bin_width_major)))
    major_dec.sort()
    major_ra = np.full_like(major_dec, ctr_ra)
    axis = np.column_stack((major_ra, major_dec))

    # get major axis position (relative to centre)
    bin_position_list = major_dec - ctr_dec
    
    # conver to relative coordinates
    axis[:,0] -= ctr_ra
    axis[:,1] -= ctr_dec
    
    # rotate according to position angle
    theta = - np.deg2rad(posang)
    rotation_matrix = np.array([[np.cos(theta), -np.sin(theta)], 
                                [np.sin(theta),  np.cos(theta)]])
    major_axis = np.full_like(axis, np.nan)
    for i in range(axis.shape[0]):
        major_axis[i,:] = np.dot(rotation_matrix, axis[i,:])
    
    # convert back to absolute coordinates
    major_axis[:,0] += ctr_ra
    major_axis[:,1] += ctr_dec
    
    # convert to world coordinates
    major_axis_pix = wcs.wcs_world2pix(major_axis, 1)
    
    ########################################
    # define size of bins
    size_x = bin_width_minor
    size_y = bin_width_major

    # make list for bins
    bin_pixels_list = []
    
    # create bins and select pixel inside bin
    for x, y in zip(major_axis[:,0], major_axis[:,1]):
    
        # create vertices of the bin
        verts = np.array([
            [x - size_x/2, y - size_y/2],  # bottom left
            [x - size_x/2, y + size_y/2],  # top left
            [x + size_x/2, y + size_y/2],  # top right
            [x + size_x/2, y - size_y/2],  # bottom right
            [x - size_x/2, y - size_y/2],  # bottom left
        ])
        
        # conver to relative coordinates
        verts[:,0] -= x
        verts[:,1] -= y
        
        # rotate according to position angle
        verts_rot = np.full_like(verts, np.nan)
        for i in range(verts.shape[0]):
            verts_rot[i,:] = np.dot(rotation_matrix, verts[i,:])
        
        # convert back to absolute coordinates
        verts_rot[:,0] += x
        verts_rot[:,1] += y
        
        # convert to pixel coordinates
        verts_world = wcs.wcs_world2pix(verts_rot, 1)
                
        # define code to link vertices
        codes = [
            Path.MOVETO,
            Path.LINETO,
            Path.LINETO,
            Path.LINETO,
            Path.CLOSEPOLY,
        ]
        
        # build hexagon path
        path = Path(verts_world, codes)
            
        # get pixels inside bin
        ra_v, dec_v = np.meshgrid(np.arange(n_ra), np.arange(n_dec))  # create meshgrid for pixel indeces in x and y
        pixels = np.column_stack((ra_v.flatten(), dec_v.flatten()))  # create list of pixels (x,y)
        bin_mask = path.contains_points(pixels)  # create mask to select pixels inside bin
        bin_pixels = pixels[bin_mask]  # select pixels of map inside bin

        # append to mask
        bin_pixels_list.append(bin_pixels)

    return bin_position_list, bin_pixels_list



def get_bin_spectra(cube, major_axis_pixels):

    """
    Takes a cube and its major axis and returns binned spectra along the major axis
    of the galaxy.
    """

    # make list for average bin spectra
    bin_spectra = []
    
    # loop over bins
    for bin_pixels in major_axis_pixels:
    
        try:
            # select spectra inside bin
            spectra = cube[:, bin_pixels[:,0], bin_pixels[:,1]]
            
            # compute average spectrum inside bin
            spectrum = np.nanmean(spectra, axis=1)
            
        except:
            spectrum = np.full_like(cube[:,0,0], np.nan)
        
        # append to list
        bin_spectra.append(spectrum)

    return bin_spectra



def get_pv_data(fits_cube, ctr_ra, ctr_dec, posang, bin_width_major=None, bin_width_minor=None, bin_width_velocity=None, rgal_axis_length=None, fits_cube_mask=None):
    
    """
    Takes a galaxy's fits cube and its coordinates and return a PV diagram array.
    """
    
    # get data and header
    cube, header = fits.getdata(fits_cube, header=True)

    if fits_cube_mask != None:
        cube_mask = fits.getdata(fits_cube_mask)
        cube[cube_mask == 0] = np.nan
    
    # build velocity axis
    n_ch = header['NAXIS3'] # number of channels
    v_ch = header['CDELT3'] # channel width (m/s, with sign)
    v_ref = header['CRVAL3'] # referecne channel (m/s)
    ch_ref = header['CRPIX3'] # reference channel
    v_ch0 = v_ref - (ch_ref-1)*v_ch # velocity of first channel
    vaxis_kms = np.linspace(v_ch0, v_ch0+(n_ch-1)*v_ch, n_ch) * 1e-3  # in km/s
    # vaxis_kms = np.arange(v_ch0, v_ch0+n_ch*v_ch, v_ch) * 1e-3  # in km/s

    # position bin width
    if bin_width_major == None:
        bin_width_major = header['BMAJ']
    if bin_width_minor == None:
        bin_width_minor = min(header['CDELT1']*header['NAXIS1'], header['CDELT2']*header['NAXIS2'])/2

    # major axis length
    if rgal_axis_length == None:

        # make 2D header
        del header['*3*']
        header['NAXIS'] = 2
        header['WCSAXES'] = 2
        
        # get coordinates
        wcs = WCS(header)
        
        # get map dimensions
        n_ra = header['NAXIS1']
        n_dec = header['NAXIS2']
        n_max = max(n_ra, n_dec)  # extend of the map in any direction
    
        # convert centre to pixel coordinatesx
        ctr_ra_pix, ctr_dec_pix = wcs.wcs_world2pix([[ctr_ra, ctr_dec]], 1)[0]  

        # compute radial extend of major axis
        dec_max = wcs.wcs_pix2world([[ctr_ra_pix, n_max]], 1)[0][1]
        rgal_axis_length = np.abs(dec_max - ctr_dec)  # degrees



    # get major axis bins
    position, major_axis_pixels = get_major_axis_bins(ctr_ra, ctr_dec, posang, header, bin_width_major, bin_width_minor, rgal_axis_length)
    
    # get major axis spectra
    spectra = get_bin_spectra(cube, major_axis_pixels)

    if bin_width_velocity == None:
        # velocity resolution is native, so skip velocity binning
        pv_data = np.transpose(np.array(spectra))
        velocity = vaxis_kms
        
    else:
        # compute velocity bins
        vaxis_sign = 1 if vaxis_kms[0] < vaxis_kms[-1] else -1
        vaxis_bin_edges = np.arange(vaxis_kms[0]-0.5*bin_width_velocity, vaxis_kms[-1]+vaxis_sign*2*bin_width_velocity, vaxis_sign*bin_width_velocity)
    
        # compute velocity bin centres
        velocity = vaxis_bin_edges[:-1] + (vaxis_bin_edges[2]-vaxis_bin_edges[1])/2 
        
        # get bin indeces
        bin_indeces = np.digitize(vaxis_kms, bins=vaxis_bin_edges, right=True)
        
        # get dimensions of pv plot
        n_pos = len(position)
        n_vel = len(velocity)
        
        # creata pv-data array
        pv_data = np.ones([n_vel, n_pos]) * np.nan
        
        # loop over spectra
        id_spec = 0
        for spec in spectra:
    
            # make lists to store binned data
            vaxis_bin_list = []  # velocity bin mean
            spec_bin_list = []  # intensity bin mean
        
            # loop over bin indeces
            for bin_idx in range(n_vel):
        
                # select data in bin
                bin_mask = bin_indeces == bin_idx
        
                # compute bin mean
                vaxis_bin = np.nanmean(vaxis_kms[bin_mask])
                
                # check if spectrum is empty
                if sum(np.isnan(spec)) == len(spec):
                    spec_bin = np.full_like(vaxis_bin, np.nan)
                else:
                    spec_bin = np.nanmean(spec[bin_mask])
        
                # append to list
                vaxis_bin_list.append(vaxis_bin)
                spec_bin_list.append(spec_bin)
        
            # cast list to numpy array
            vaxis_bin = np.array(vaxis_bin_list)
            spec_bin = np.array(spec_bin_list)
    
            # add spectrum to pv array
            pv_data[:, id_spec] = spec_bin
            id_spec += 1

    # find position indices that contain only nans
    id_del = []  # position indeces that contain only nan values
    id_spec = 0
    for spec in pv_data.T:
        if sum(np.isnan(spec)) == len(spec):
            id_del.append(id_spec)
        id_spec += 1
    # trim data
    pv_data = np.delete(pv_data, id_del, axis=1)
    position = np.delete(position, id_del)
    
    return pv_data, position, velocity
