# Imports. This is a python, not a CASA module.

import os, glob
import numpy as np

import spectral_cube as sc

from astropy import units as u
from astropy.wcs import WCS
from astropy.stats import mad_std

from utilsImages import *

def make_linear_mosaic(
        input_file_list,
        reponse_file_list = None,
        response_type = 'pb',
        global_weight_list = None,
        global_weight_from_mad = True,
        target_header = None,
        template_file = None,
        target_beam = None,        
        skip_convolve = False
        tol_convolve = 0.05,):
    
    """Make a linear mosaic from a list of input files.
    
    """

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Work out the target and initialize output
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    if target_header is None and template_file is not None:
        template_sc = sc.SpectralCube.read(template_file)
        target_header = target_sc.header
    elif target_header is None:
        pass

    # Initialize output arrays to hold the sum of weight*data and weights
    
    output_shape = (target_header['NAXIS3'], target_header['NAXIS2'], target_header['NAXIS1'])        
    sum_of_data = np.array(output_shape)
    sum_of_weights = np.array(output_shape)

    # If needed, find the common beam for the stack of input images
    
    if target_beam is None and skip_convolve == False:
        
        target_beam = common_beam_from_flist(
            input_file_list, round_beam=True)
        
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Loop over input
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    for ii, this_file in input_file_list:
        
        # Read the data
        this_sc = sc.SpectralCube.read(this_file)
        this_sc.allow_huge_operations = True
        
        # Detect if the data is a cube or an image

        print("Need cube vs. image detection.")

        #this_tp = sc.SpectralCube.read(this_file)
        #rest_freq_line = 230.538*u.GHz
        #this_tp_vel = this_tp.with_spectral_unit(
        #u.m/u.s, velocity_convention='radio', rest_value=rest_freq_line)

        # ......................................................
        # Convolve the data to target resolution
        # ......................................................
        
        if (skip_convolve == False) and (target_beam is not None):

            proceed_with_convolve = True
            
            if tol_convolve > 0.0:
                
                proceed_with_convolve = \
                    beams_match(target_beam, this_sc.beam, tol=tol_convolve,
                                check_pa = True, pa_tol = 5.*u.deg)

            if proceed_with_convolve:
                this_sc = this_sc.convolve_to(target_beam)
        
        # ...........................
        # Response and weight
        # ...........................
        
        # Get the response for the image. This is a spatially variable
        # value that weights the mosaic at each position. The optimal
        # weighting scheme goes as 1/noise**2 at each location. We
        # allow three input conventions:

        # 'pb' : The response is a "primary beam" of the type output for
        # interferometer data. The weight is proportional to
        # (primary_beam)**2 because the primary beam is proportional
        # to 1/noise.

        # 'noise' : The response is a "noise cube" or image giving the
        # estimated noise at each location in the data set. Then the
        # weight is proportional to 1/noise**2.

        # 'weight' : the response is a positional weight and directly
        # read as such.

        # If nothing is supplied then default to a response of 1.0 at
        # each location.
        
        if response_file_list is not None:
            
            this_response_file = response_file_list[ii]
            this_response_sc = sc.SpectralCube.read(this_response_file)

            if response_type == 'pb':
                weight_sc = this_response_sc**2

            if response_type == 'noise':
                weight_sc = 1./this_response_sc**2
            
            if response_type == 'weight':
                weight_sc = this_response_sc

        else:

            weight_sc = np.isfinite(this_sc)                        
            
        # Get the global weight for this image. This is a single value
        # per image (with the response file handling the spatial
        # variations).

        # If no global weight list is supplied AND the response type
        # is 'pb' then the global weights are calculated by estimating
        # the noise from a mad-based standard deviation of the data
        # which is converted to a weight by 1/noise**2. This is the
        # desired behavior if the response is normalized, e.g., as
        # expected for a primary beam response input.

        # If no weight list is supplied and this option is disabled
        # then the program will weight each image equally.
        
        if global_weight_list is not None:
            this_global_weight = global_weight_list[ii]
        else:
            
            if response_type == 'pb' and global_weight_from_mad:
                this_noise = mad_std(this_data.filled(fill_value=np.nan))
                this_global_weight = 1./this_noise**2
            else:
                this_global_weight = 1.0
                
        # Combine the response and global weight to give the scaling
        # for the mosaic.

        this_weight_image = this_resolved_weight * this_global_weight
        this_weight_image.allow_huge_operations = True
        
        # ...........................
        # Reproject to the output grid
        # ...........................

        if working_with_cube:            
            this_sc_vel_grid = this_sc.spectral_interpolate(
                target_spec_axis, fill_value=np.nan)
            this_weight_vel_grid = this_weight_image.spectral_interpolate(
                target_spec_axis, fill_value=np.nan)
        else:
            this_sc_vel_grid = this_sc
            this_weight_vel_grid = this_weight_image
            
        this_sc_output_grid = this_sc_vel_grid.reproject(
            target_header, order='bilinear')
        this_weight_output_grid = this_weight_vel_grid.reproject(
            target_header, order='bilinear')

        this_sc_output_grid.allow_huge_operations = True
        this_weight_output_grid.allow_huge_operations = True        
        
        # ...........................
        # Compile output
        # ...........................

        this_sc_values = this_sc_output_grid.filled(fill_value=np.nan)
        this_weight_values = this_weight_output_grid(fill_value=np.nan)

        mask = np.where(np.isfinite(this_sc_values))
        sum_of_data[mask] = sum_of_data[mask] + this_sc_values[mask]*this_weight_values[mask]
        sum_of_weight[mask] = sum_of_weight[mask] + this_weight_values[mask]
    
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Clean up, write output, and return
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    # Divide the sum by the weights
    output_data = sum_of_data / sum_of_weight

    # Convert to spectral cubes
    output_sc = sc.SpectralCube(output_data, wcs=WCS(target_header), beam=target_beam)
    output_weight_sc = sc.SpectralCube(sum_of_weight, wcs=WCS(target_header), beam=target_beam)

    # Write output
    if outfile is not None:
        output_sc.write(outfile, overwrite=overwrite)
    
    if outfile_weight is not None:
        output_weight_sc.write(outfile_weight, overwrite=overwrite)    

    # Return the mosaicked cube and weight image
    return(output_sc, output_weight_sc)
