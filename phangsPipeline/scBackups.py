
def write_moment1_hybrid(
        cube, rms=None, channel_correlation=None,
        outfile=None, errorfile=None,
        overwrite=True, unit=None,
        return_products=True,
        strict_vfield=None,
        broad_vfield=None,
        broad_signal=None,
        vfield_prior=None,
        vfield_prior_res=None,
        vfield_reject_thresh='30km/s',
        mom0_thresh_for_mom1=2.0,
        context=None):
    """Write out moment1 map using combination of other moment maps.
    This is a secondary moment that needs to be calculated in the
    context of other moments.
    
    Keywords:
    ---------
    
    cube : SpectralCube
        Included to keep same call signature but not used
    
    outfile : str
        File name of output file
        
    errorfile : str
        File name of map for the uncertainty
        
    rms : SpectralCube
        Included to keep the same call signature but not used.
        
    channel_correlation : np.array
        Included to keep the same call signature but not used.

    overwrite : bool
        Set to True (the default) to overwrite existing maps if present. 
        
    unit : astropy.Unit
        Preferred unit for moment masks
        
    return_products : bool
        Return products calculated in the map

    strict_vfield : str
        Moment tag for velocity field to be used as a high confidence map

    broad_vfield : str
        Moment tag for velocity field for low confidence map

    broad_signal : str
        Moment tag to be used as an estimate of the signal for a S/N
        cut on where the broad_vfield is valid.  Also finds a noise
        estimate of the same and uses this for the Noise component

    vfield_prior : str
        Moment tag for low-resolution prior map of velocity field

    vfield_prior_res : str
        Resolution tag for low-resolution prior map of velocity field

    vfield_reject_thresh : astropy.units.Quantity
        The maximum difference between the broad field and the prior
        field in units that can convert to that of the velocity field.

     mom0_thresh_for_mom1 : float
        S/N threshold for using a broad_vfield estimate in the map
    """


    # Resolution to work with
    resname = context['res_tag']
    if resname is None:
        resname = ''

    # The threshold for outlier rejection from the prior velocity field

    vfield_reject_thresh = u.Quantity(vfield_reject_thresh)

    # The root for the maps

    moment_root = utilsFilenames.get_cube_filename(
        target=context['target'], config=context['config'],
        product=context['product'],
        ext=resname + context['extra_ext'])
    moment_root = moment_root.replace('.fits','')

    # This strict moment1 field will remain in place no matter what

    strict_moment1_name = ''.join([context['indir'],
                                   moment_root,
                                   context['allmoments'][strict_vfield]['ext'],
                                   '.fits'])
    
    mom1strict = convert_and_reproject(strict_moment1_name, unit=unit)

    strict_moment1_err_name = ''.join([context['indir'],
                                       moment_root,
                                       context['allmoments'][strict_vfield]['ext_error'],
                                       '.fits'])
    
    mom1strict_error = convert_and_reproject(strict_moment1_err_name, unit=unit)

    # This broad moment 0 map will be used to help prune faint emission

    broad_moment0_name = ''.join([context['indir'],
                                  moment_root,
                                  context['allmoments'][broad_signal]['ext'],
                                  '.fits'])

    mom0broad = convert_and_reproject(broad_moment0_name, template=mom1strict)

    broad_moment0_err_name = ''.join([context['indir'],
                                      moment_root,
                                      context['allmoments'][broad_signal]['ext_error'],
                                      '.fits'])

    mom0broad_error = convert_and_reproject(broad_moment0_err_name, template=mom1strict)

    # This broad moment 1 map will be used as a candidate velocity field

    broad_moment1_name = ''.join([context['indir'],
                                  moment_root,
                                  context['allmoments'][broad_vfield]['ext'],
                                  '.fits'])

    mom1broad = convert_and_reproject(broad_moment1_name, template=mom1strict,
                                      unit=unit)

    broad_moment1_err_name = ''.join([context['indir'],
                                      moment_root,
                                      context['allmoments'][broad_vfield]['ext_error'],
                                      '.fits'])
    
    mom1broad_error = convert_and_reproject(broad_moment1_err_name, template=mom1strict,
                                            unit=unit)    

    # This prior velocity field will be used to reject outliers

    resname = vfield_prior_res

    # AKL - need to make the config tunable and separate from the input maps here

    moment_root = utilsFilenames.get_cube_filename(
        target=context['target'], config=context['config'],
        product=context['product'],
        ext=resname + context['extra_ext'])

    moment_root = moment_root.replace('.fits','')

    prior_moment1_name = ''.join([context['indir'],
                                  moment_root,
                                  context['allmoments'][vfield_prior]['ext'],
                                  '.fits'])

    
    mom1prior = convert_and_reproject(prior_moment1_name, template=mom1strict,
                                      unit=unit)

    # Now hybridize

    # ... start with the high quality strict mask
    mom1hybrid = mom1strict.value

    # ... candidate entries are places with a broad value
    valid_broad_mom1 = np.isfinite(mom1broad.value)

    # ... but not any strict value
    valid_broad_mom1[np.isfinite(mom1strict)] = False

    # If thresholding on intensity, apply that

    if mom0broad_error is not None:
        valid_broad_mom1 *= (mom0broad.value
                             > (mom0_thresh_for_mom1
                             * mom0broad_error.value))

    # If thresholding relative to prior field, apply that
        
    if mom1prior is not None:
        valid_broad_mom1 = (valid_broad_mom1 *
                            (np.abs(mom1broad - mom1prior)
                             < vfield_reject_thresh)
                            )
    
    # Fill in the still-valid locations in the hybrid

    mom1hybrid[valid_broad_mom1] = (mom1broad.value)[valid_broad_mom1]
    mom1hybrid = u.Quantity(mom1hybrid, unit)
    if unit is not None:
        mom1hybrid = mom1hybrid.to(unit)
    
    # Attach to WCS

    mom1hybrid_proj = Projection(mom1hybrid,
                                 wcs=mom1strict.wcs,
                                 header=mom1strict.header,
                                 meta=mom1strict.meta)

    # Write

    if outfile is not None:
        mom1hybrid_proj.write(outfile,
                              overwrite=overwrite)

    # Propagate errors from the input map to an error map

    mom1hybrid_error = None
    
    if (type(mom1broad_error) is Projection and 
        type(mom1strict_error) is Projection):
        mom1hybrid_error = mom1broad_error
        mom1hybrid_error[~np.isfinite(mom1hybrid.value)] = np.nan
        strictvals = np.isfinite(mom1strict_error.value)
        mom1hybrid_error[strictvals] = mom1strict_error[strictvals]
        if unit is not None:
            mom1hybrid_error = mom1hybrid_error.to(unit)
        mom1hybrid_error_proj = Projection(mom1hybrid_error,
                                           wcs=mom1strict.wcs,
                                           header=mom1strict.header,
                                           meta=mom1strict.meta)
        if errorfile is not None:
            mom1hybrid_error_proj = update_metadata(mom1hybrid_error_proj,
                                                    cube, error=True)
            mom1hybrid_error_proj.write(errorfile,
                                        overwrite=overwrite)
    
    if return_products and mom1hybrid_error_proj is not None:
        return(mom1hybrid_proj, mom1hybrid_error_proj)
    elif return_products and mom1hybrid_error_proj is None:
        return(mom1hybrid_proj)

    
def old_write_moment1_hybrid(cube,
                             broad_mask=None,
                             moment1_prior=None,
                             order='bilinear',
                             outfile=None,
                             errorfile=None,
                             rms=None,
                             channel_correlation=None,
                             overwrite=True,
                             vfield_reject_thresh=30 * u.km / u.s,
                             mom0_thresh_for_mom1=2.0,
                             unit=None,
                             return_products=False):
    """
    Writes a moment 1 map 
    
    Parameters:
    -----------
    
    cube : SpectralCube
        SpectralCube of original data with strict masking applied
    
    Keywords:
    ---------
    
    broad_mask : SpectralCube or np.array
        Array with same shape as the input SpectralCube to be used 
        as the broad (permissive) mask
        
    moment1_prior : FITS filename or Projection
        FITS filename or Projection containting the velocity field prior
        
    order : str
        Specifies the order of interpolation to be used for aligning spectral 
        cubes to each other from 'nearest-neighbor', 'bilinear', 
        'biquadratic', 'bicubic'. Defaults to 'bilinear'.
    
    errorfile : str
        File name of map for the uncertainty
        
    rms : SpectralCube
        Root-mean-square estimate of the error.  This must have an estimate
        the noise level at all positions where there is signal, and only at 
        those positions.
        
    channel_correlation : np.array
        One-dimensional array containing the channel-to-channel 
        normalize correlation coefficients
        
    overwrite : bool
        Set to True (the default) to overwrite existing maps if present. 
        
    unit : astropy.Unit
        Preferred unit for moment masks
        
    vfield_reject_thresh : astropy.Quantity
        Velocity range beyond which deviations from a prior velocity 
        field are rejected.  Default 30 km/s
    
    mom0_thresh_for_mom1 : int
        Signal-to-noise ratio in a moment-0 to accept a measurement 
        of a moment1 map.
        
    return_products : bool
        Return products calculated in the map
    """

    (mom1strict, 
     mom1strict_error) = write_moment1(cube, rms=rms,
                                       channel_correlation=channel_correlation,
                                       unit=unit,
                                       return_products=True)

    spaxis = cube.spectral_axis.value

    if moment1_prior is not None:
        if type(moment1_prior) is Projection:
            mom1prior = moment1_prior
        elif type(moment1_prior) is str:
            hdu_list = fits.open(moment1_prior)
            mom1prior = Projection.from_hdu(hdu_list[0])
        mom1prior = mom1prior.to(cube.spectral_axis.unit)
        mom1prior = mom1prior.reproject(mom1strict.header, order=order)
    else:
        mom1prior = None

    if type(broad_mask) is SpectralCube:
        strict_mask = SpectralCube(cube.mask.include(),
                                   wcs=cube.wcs,
                                   header=cube.header)
        hybrid_mask = hybridize_mask(strict_mask,
                                     broad_mask,
                                     return_cube=False)
        broad_cube = cube.with_mask(hybrid_mask,
                                    inherit_mask=False)
        
    elif type(broad_mask) is str:
        broad_mask = SpectralCube.read(broad_mask)
        strict_mask = SpectralCube(cube.mask.include(),
                                   wcs=cube.wcs,
                                   header=cube.header)
        hybrid_mask = hybridize_mask(strict_mask,
                                     broad_mask,
                                     return_cube=False)
        broad_cube = cube.with_mask(hybrid_mask,
                                    inherit_mask=False)
        
    elif type(broad_mask) is np.ndarray:
        broad_cube = cube.with_mask(broad_mask.astype(bool),
                                    inherit_mask=False)

    (mom0broad,
     mom0broad_error) = write_moment0(broad_cube, rms=rms,
                                      channel_correlation=channel_correlation,
                                      return_products=True)
     
    (mom1broad,
     mom1broad_error) = write_moment1(broad_cube, rms=rms,
                                      channel_correlation=channel_correlation,
                                      unit=unit,
                                      return_products=True)
    
    mom1hybrid = mom1strict.value
    valid_broad_mom1 = np.isfinite(mom1broad.value)
    valid_broad_mom1[np.isfinite(mom1strict)] = False

    if mom0broad_error is not None:
        valid_broad_mom1 *= (mom0broad.value
                             > (mom0_thresh_for_mom1
                             * mom0broad_error.value))
        
    if mom1prior is not None:
        valid_broad_mom1 = (valid_broad_mom1 *
                            (np.abs(mom1broad - mom1prior)
                             < vfield_reject_thresh)
                            )
    
    mom1hybrid[valid_broad_mom1] = (mom1broad.value)[valid_broad_mom1]
    mom1hybrid = u.Quantity(mom1hybrid, cube.spectral_axis.unit)
    if unit is not None:
        mom1hybrid = mom1hybrid.to(unit)
    
    mom1hybrid_proj = Projection(mom1hybrid,
                                 wcs=mom1strict.wcs,
                                 header=mom1strict.header,
                                 meta=mom1strict.meta)
    if outfile is not None:
        mom1hybrid_proj = update_metadata(mom1hybrid_proj, cube)
        mom1hybrid_proj.write(outfile,
                              overwrite=overwrite)
    mom1hybrid_error = None
    
    if (type(mom1broad_error) is Projection and 
        type(mom1strict_error) is Projection):
        mom1hybrid_error = mom1broad_error
        mom1hybrid_error[~np.isfinite(mom1hybrid.value)] = np.nan
        strictvals = np.isfinite(mom1strict_error.value)
        mom1hybrid_error[strictvals] = mom1strict_error[strictvals]
        if unit is not None:
            mom1hybrid_error = mom1hybrid_error.to(unit)
        mom1hybrid_error_proj = Projection(mom1hybrid_error,
                                           wcs=mom1strict.wcs,
                                           header=mom1strict.header,
                                           meta=mom1strict.meta)
        if errorfile is not None:
            mom1hybrid_error_proj = update_metadata(mom1hybrid_error_proj,
                                                    cube, error=True)
            mom1hybrid_error_proj.write(errorfile,
                                        overwrite=overwrite)
    
    if return_products and mom1hybrid_error_proj is not None:
        return(mom1hybrid_proj, mom1hybrid_error_proj)
    elif return_products and mom1hybrid_error_proj is None:
        return(mom1hybrid_proj)
