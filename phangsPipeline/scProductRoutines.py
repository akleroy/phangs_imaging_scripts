from spectral_cube import SpectralCube, Projection
import astropy.units as u
# from pipelineVersion import version as pipeVer
from scMaskingRoutines import noise_cube, simple_mask
import numpy as np

import logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


def write_moment0(cube,
                  outfile=None,
                  errorfile=None,
                  rms=None, 
                  channel_correlation=None,
                  overwrite=True,
                  unit=None):
    """
    Write out moment0 map for a SpectralCube
    
    Keywords:
    ---------
    
    cube : SpectralCube
        (Masked) spectral cube to write a moment0 map
    
    outfile : str
        File name of output file
        
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
    
    """

    mom0 = cube.moment0()
    if unit is not None:
        mom0 = mom0.to(unit)
    mom0.write(outfile, overwrite=overwrite)

    if errorfile is not None and rms is None:
        logger.error("Moment 0 error requested but no RMS provided")

    if rms is not None and errorfile is not None:

        if channel_correlation is None:
            channel_correlation = np.array([1])

        mom0err = np.empty(mom0.shape)
        mom0err.fill(np.nan)
  
        for x, y, slc in cube._iter_rays(0):
            mask = np.squeeze(cube._mask.include(view=slc))
            if not mask.any():
                continue
            index = np.where(mask)[0]
            rms_spec = rms.flattened(slc).value
            spec = cube.flattened(slc).value
            
            covar = build_covariance(spectrum=spec,
                                     rms=rms_spec,
                                     channel_correlation=channel_correlation,
                                     index=index)

            mom0err = (np.sum(covar**2))**0.5
        mom0err = u.Quantity(mom0err, mom0.unit, copy=False)
        if unit is not None:
            mom0err = mom0err.to(unit)
        mom0projection = Projection(mom0err,
                                    wcs=mom0.wcs,
                                    header=mom0.header)
        
        mom0projection.write(errorfile, overwrite=overwrite)
    #     return(mom0, mom0err)
    # else:
    #     return(mom0, None)


def write_moment1(cube,
                  outfile=None,
                  errorfile=None,
                  rms=None,
                  channel_correlation=None,
                  overwrite=True,
                  unit=None):
    """
    Write out moment1 map for a SpectralCube
    
    Keywords:
    ---------
    
    cube : SpectralCube
        (Masked) spectral cube to write a moment1 map
    
    outfile : str
        File name of output file
        
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
    
    """

    mom1 = cube.moment1()
    if unit is not None:
        mom1 = mom1.to(unit)
        spaxis = cube.spectral_axis.to(unit).value
    else:
        spaxis = cube.spectral_axis.value
    mom1.write(outfile, overwrite=True)

    if errorfile is not None and rms is None:
        logger.error("Moment 1 error requested but no RMS provided")

    if rms is not None and errorfile is not None:

        if channel_correlation is None:
            channel_correlation = np.array([1])

        mom1err = np.empty(mom1.shape)
        mom1err.fill(np.nan)
        # Ensure the same mask applied to both.
        rms = rms.with_mask(cube._mask, inherit_mask=False)
        for x, y, slc in cube._iter_rays(0):
            mask = np.squeeze(cube._mask.include(view=slc))
            if not mask.any():
                continue
            index = np.where(mask)[0]
            rms_spec = rms.flattened(slc).value
            spec = cube.flattened(slc).value
            covar = build_covariance(spectrum=spec,
                                     rms=rms_spec,
                                     channel_correlation=channel_correlation,
                                     index=index)

            vval = spaxis[index]
            sum_T = np.sum(spec)
            sum_vT = np.sum(spec * vval)
            
            jacobian = vval / sum_T - sum_vT / sum_T**2
            mom1err[x, y] = np.dot(
                np.dot(jacobian[np.newaxis, :], covar),
                jacobian[:, np.newaxis])**0.5
        mom1err = u.Quantity(mom1err, mom1.unit, copy=False)
        if unit is not None:
            mom1err = mom1err.to(unit)
        mom1projection = Projection(mom1err,
                                    wcs=mom1.wcs,
                                    header=mom1.header)

        mom1projection.write(errorfile, overwrite=overwrite)
    #     return(mom1, mom1err)
    # else:
    #     return(mom1, None)


def write_moment2(cube,
                  outfile=None,
                  errorfile=None,
                  rms=None,
                  channel_correlation=None,
                  overwrite=True,
                  unit=None):
    """
    Write out linewidth (moment2-based) map for a SpectralCube
    
    Keywords:
    ---------
    
    cube : SpectralCube
        (Masked) spectral cube to write a Moment 2 (velocity 
        dispersion) map
    
    outfile : str
        File name of output file
        
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
    
    """

    mom2 = cube.linewidth_sigma()
    if unit is not None:
        mom2 = mom2.to(unit)
        spaxis = cube.spectral_axis.to(unit).value
    else:
        spaxis = cube.spectral_axis.value
    mom2.write(outfile, overwrite=True)

    if errorfile is not None and rms is None:
        logger.error("Moment 2 error requested but no RMS provided")

    if rms is not None and errorfile is not None:

        if channel_correlation is None:
            channel_correlation = np.array([1])

        mom2err = np.empty(mom2.shape)
        mom2err.fill(np.nan)

        for x, y, slc in cube._iter_rays(0):
            mask = np.squeeze(cube._mask.include(view=slc))
            if not mask.any():
                continue
            index = np.where(mask)[0]
            rms_spec = rms.flattened(slc).value
            spec = cube.flattened(slc).value
            covar = build_covariance(spectrum=spec,
                                     rms=rms_spec,
                                     channel_correlation=channel_correlation,
                                     index=index)

            vval = spaxis[index]
            sum_T = np.sum(spec)
            sum_vT = np.sum(spec * vval)
            vbar = sum_vT / sum_T
            vdisp = (vval - vbar)**2
            wtvdisp = np.sum(spec * vdisp)
            # Dear future self: There is no crossterm (error term from 
            # vbar) since dispersion is at a minimum around vbar
            jacobian = (vdisp / sum_T
                        - wtvdisp / sum_T**2)
            mom2err[x, y] = np.dot(
                np.dot(jacobian[np.newaxis, :], covar),
                jacobian[:, np.newaxis])**0.5
        mom2err = u.Quantity(mom2err, mom2.unit, copy=False)
        if unit is not None:
            mom2err = mom2err.to(unit)
        mom2projection = Projection(mom2err,
                                    wcs=mom2.wcs,
                                    header=mom2.header)

        mom2projection.write(errorfile, overwrite=overwrite)
    #     return(mom2, mom2err)
    # else:
    #     return(mom2, None)


def write_ew(cube,
             outfile=None,
             errorfile=None,
             rms=None,
             channel_correlation=None,
             overwrite=True,
             unit=None):
    """
    Write out linewidth (equivalent-width-based) map for a SpectralCube
    
    Keywords:
    ---------
    
    cube : SpectralCube
        (Masked) spectral cube to write a EW-based velocity dispersion map
    
    outfile : str
        File name of output file
        
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
    
    """

    maxmap = cube.max(axis=0)
    mom0 = cube.moment0()
    sqrt2pi = np.sqrt(2 * np.pi)
    sigma_ew = mom0 / maxmap / sqrt2pi

    if unit is not None:
        sigma_ew = sigma_ew.to(unit)
        spaxis = cube.spectral_axis.to(unit).value
    else:
        spaxis = cube.spectral_axis.value
    sigma_ew.write(outfile, overwrite=True)

    if errorfile is not None and rms is None:
        logger.error("Equivalent width error requested but no RMS provided")

    if rms is not None and errorfile is not None:
        argmaxmap = cube.argmax(axis=0)
        rms_at_max = np.take_along_axis(
            rms.filled_data[:], 
            argmaxmap[np.newaxis, :, :], 0).value
        
        if channel_correlation is None:
            channel_correlation = np.array([1])

        sigma_ew_err = np.empty(sigma_ew.shape)
        sigma_ew_err.fill(np.nan)


        for x, y, slc in cube._iter_rays(0):
            mask = np.squeeze(cube._mask.include(view=slc))
            if not mask.any():
                continue
            index = np.where(mask)[0]
            rms_spec = rms.flattened(slc).value
            spec = cube.flattened(slc).value
            covar = build_covariance(spectrum=spec,
                                     rms=rms_spec,
                                     channel_correlation=channel_correlation,
                                     index=index)
            sigma_ew_err[x, y] = (np.sum(covar**2) 
                                  + (sigma_ew[x, y].value**2 
                                     * rms_at_max[0, x, y]**2 
                                     / maxmap[x, y].value**2))**0.5

        sigma_ew_err = u.Quantity(sigma_ew_err, sigma_ew.unit, copy=False)
        if unit is not None:
            sigma_ew_err = sigma_ew_err.to(unit)
        sigma_ewerr_projection = Projection(sigma_ew_err,
                                         wcs=sigma_ew.wcs,
                                         header=sigma_ew.header)

        sigma_ewerr_projection.write(errorfile, overwrite=overwrite)


def write_tmax(cubein,
               outfile=None,
               errorfile=None,
               rms=None,
               channel_correlation=None,
               overwrite=True,
               unit=None,
               window=None):
    """
    Write out Tmax map for a SpectralCube
    
    Keywords:
    ---------
    
    cube : SpectralCube
        (Masked) spectral cube to write a Tmax map
    
    outfile : str
        File name of output file
        
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
        
    window : astropy.Quantity
        Spectral window over which the data should be smoothed
    """
    if type(window) is u.Quantity:
        from astropy.convolution import Box1DKernel
        channel_width = np.abs(cubein.spectral_axis[0]
                               - cubein.spectral_axis[1])
        nChan = (window / channel_width).to(u.dimensionless_unscaled).value
        if nChan > 1:
            cube = cubein.spectral_smooth(Box1DKernel(nChan))
        else:
            cube = cubein
    else:
        cube = cubein
    maxmap = cube.max(axis=0)


    maxmap.write(outfile, overwrite=True)

    if errorfile is not None and rms is None:
        logger.error("Tmax error requested but no RMS provided")

    if rms is not None and errorfile is not None:
        argmaxmap = cube.argmax(axis=0)
        rms_at_max = np.take_along_axis(
            rms.filled_data[:],
            argmaxmap[np.newaxis, :, :], 0).value
        rms_at_max = np.squeeze(rms_at_max)
        rms_at_max = u.Quantity(rms_at_max, maxmap.unit, copy=False)
        if unit is not None:
            rms_at_max = rms_at_max.to(unit)
        rms_map_projection = Projection(rms_at_max,
                                        wcs=maxmap.wcs,
                                        header=maxmap.header)

        rms_map_projection.write(errorfile, overwrite=overwrite)


def write_vmax(cubein,
               outfile=None,
               errorfile=None,
               rms=None,
               channel_correlation=None,
               overwrite=True,
               unit=None,
               window=None):
    """
    Write out velocity map at max brightness temp for a SpectralCube
    
    Keywords:
    ---------
    
    cube : SpectralCube
        (Masked) spectral cube to write a Vmax map
    
    outfile : str
        File name of output file
        
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
    
    window : astropy.Quantity
        Spectral window over which the data should be smoothed
    """
    if type(window) is u.Quantity:
        from astropy.convolution import Box1DKernel
        channel_width = np.abs(cubein.spectral_axis[0]
                               - cubein.spectral_axis[1])
        nChan = (window / channel_width).to(u.dimensionless_unscaled).value
        if nChan > 1:
            cube = cubein.spectral_smooth(Box1DKernel(nChan))
        else:
            cube = cubein
    else:
        cube = cubein
    maxmap = cube.max(axis=0)
    argmaxmap = cube.argmax(axis=0)
    vmaxmap = np.take_along_axis(cube.spectral_axis[:, np.newaxis, np.newaxis],
                                 argmaxmap[np.newaxis, :, :], 0)
    vmaxmap = np.squeeze(vmaxmap)
    if unit is not None:
        vmaxmap = vmaxmap.to(unit)
    vmaxmap[~np.isfinite(maxmap)] = np.nan
    vmaxmap_projection = Projection(vmaxmap,
                                    wcs=maxmap.wcs,
                                    header=maxmap.header)
    vmaxmap_projection.write(outfile, overwrite=True)

    if errorfile is not None and rms is None:
        logger.error("Moment 2 error requested but no RMS provided")

    if rms is not None and errorfile is not None:
        channel_width = np.abs(cube.spectral_axis[0] 
                               - cube.spectral_axis[1])
        vmaxerror = np.empty(maxmap.shape)
        vmaxerror.fill(np.nan)
        vmaxerror[np.isfinite(maxmap)] = channel_width.value
        vmaxerror = u.Quantity(vmaxerror, channel_width.unit, copy=False)
        if unit is not None:
            vmaxerror = vmaxerror.to(unit)
        vmaxerror_projection = Projection(vmaxerror,
                                          wcs=maxmap.wcs,
                                          header=maxmap.header)

        vmaxerror_projection.write(errorfile, overwrite=overwrite)


def write_vquad(cubein,
                outfile=None,
                errorfile=None,
                rms=None,
                channel_correlation=None,
                overwrite=True,
                unit=None,
                window=None,
                maxshift=0.5):
    """
    Write out velocity map at max brightness temp for a 
    SpectralCube using the quadratic peak interpolation
    
    Keywords:
    ---------
    
    cube : SpectralCube
        (Masked) spectral cube to write a Vmax map
    
    outfile : str
        File name of output file
        
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
    
    window : astropy.Quantity
        Spectral window over which the data should be smoothed
    
    maxshift : np.float
        Maximum number of channels that the algorithm can shift the 
        peak estimator (default = 0.5).  Set to None to suppress clipping.
    """
    
    from scipy.interpolate import interp1d
    
    if type(window) is u.Quantity:
        from astropy.convolution import Box1DKernel
        channel_width = np.abs(cubein.spectral_axis[0]
                               - cubein.spectral_axis[1])
        nChan = (window / channel_width).to(u.dimensionless_unscaled).value
        if nChan > 1:
            cube = cubein.spectral_smooth(Box1DKernel(nChan))
        else:
            cube = cubein
    else:
        cube = cubein
        
    if unit is not None:
        spaxis = cube.spectral_axis.to(unit).value
    else:
        spaxis = cube.spectral_axis.value
        unit = cube.spectral_axis.unit
    pixinterp = interp1d(np.arange(spaxis.size),
                         spaxis)
    maxmap = cube.max(axis=0)
    argmaxmap = cube.argmax(axis=0)
    Tup = np.take_along_axis(cube.filled_data[:],
                             argmaxmap[np.newaxis, :, :]+1, 0).value
    Tdown = np.take_along_axis(cube.filled_data[:],
                               argmaxmap[np.newaxis, :, :]-1, 0).value
    Tup = np.squeeze(np.nan_to_num(Tup))
    Tup[Tup < 0] = 0
    Tdown = np.squeeze(np.nan_to_num(Tdown))
    Tdown[Tdown < 0] = 0

    delta = -1 * ((Tup - Tdown) / (Tup + Tdown - 2 * maxmap.value))
    if maxshift is not None:
        delta = np.clip(delta, -maxshift, maxshift)
    peakchan = argmaxmap + delta

    vmaxmap = np.empty(maxmap.shape)
    vmaxmap.fill(np.nan)
    good = np.isfinite(maxmap)
    vmaxmap[good] = u.Quantity(pixinterp(peakchan[good]),
                               unit)

    vmaxmap[~np.isfinite(maxmap)] = np.nan
    vmaxmap_projection = Projection(vmaxmap,
                                    wcs=maxmap.wcs,
                                    header=maxmap.header)
    vmaxmap_projection.write(outfile, overwrite=overwrite)

    if errorfile is not None and rms is None:
        logger.error("Vquad error requested but no RMS provided")

    if rms is not None and errorfile is not None:
        channel_width = np.abs(np.median(cube.spectral_axis[1:] 
                                         - cube.spectral_axis[:-1])).to(unit)
        RMSup = np.take_along_axis(rms.filled_data[:],
                                   argmaxmap[np.newaxis, :, :]+1,
                                   0).value
        RMSdown = np.take_along_axis(rms.filled_data[:],
                                     argmaxmap[np.newaxis, :, :]-1,
                                     0).value
        RMSmax = np.take_along_axis(rms.filled_data[:],
                                    argmaxmap[np.newaxis, :, :],
                                    0).value
        denom = (Tup + Tdown - 2 * maxmap.value)
        j1 = (1/denom - (Tup - Tdown) / denom**2)
        j2 = (2 * (Tup - Tdown) / denom**2)
        j3 = (-1/denom - (Tup - Tdown) / denom**2)
        jacobian = np.r_[j1[np.newaxis, :, :],
                         j2[np.newaxis, :, :],
                         j3[np.newaxis, :, :]]
        if ((channel_correlation is None) 
            or len(channel_correlation) == 1):
            covar = np.r_[RMSup,
                          RMSmax,
                          RMSdown]
            covar = covar**2
            error = np.einsum('i...,i...', jacobian**2, covar)
        else:
            if len(channel_correlation) == 2:
                ccor = np.r_[channel_correlation,
                             np.array([0])]
            else:
                ccor = channel_correlation[0:3]
            corrmat = ccor[np.array([[0, 1, 2],
                                     [1, 0, 1],
                                     [2, 1, 0]])]
            rmsvec = np.r_[RMSup,
                           RMSmax,
                           RMSdown]
            # They never should have taught me how to do this.
            covar  = np.einsum('ij,ilm,jlm->ijlm', corrmat,
                               rmsvec, rmsvec)
            error = np.einsum('ilm,jlm,ijlm->lm',
                              jacobian, jacobian, covar)
        if maxshift is not None:
            error = np.clip(error, -maxshift, maxshift)
        
        vquaderror = error * channel_width
        vquaderror_projection = Projection(vquaderror,
                                           wcs=maxmap.wcs,
                                           header=maxmap.header)
        vquaderror_projection.write(errorfile, overwrite=overwrite)


def build_covariance(spectrum=None,
                     rms=None,
                     channel_correlation=None,
                     index=None):
    """
    Build a covariance matrix from a channel_correlation vector
    
    Keywords:
    ---------
    
    spectrum : np.array
        One-dimensional array of spectrum values
        
    rms : np.array
        One-dimensional array containing the root-mean-squared error
        estimate for the values in the spectrum
    
    channel_correlation : np.array
        One-dimensional array containing the channel-to-channel 
        normalize correlation coefficients
    
    index : np.array
        Integer array indicating the spectral indices of the data 
        in the original cube
    """
    
    # Note that this assumes you are masking out values to make sure 
    # arrays stay the same shape as the input

    if index is None:
        index = np.arange(len(spectrum))
    if channel_correlation is None:
        return(np.diag(rms**2))
    if len(channel_correlation) == 1:
        return(np.diag(rms**2))
    distance = np.abs(index[:, np.newaxis] - index[np.newaxis, :])
    covar = rms[:, np.newaxis] * rms[np.newaxis, :]
    maxdist = len(channel_correlation)
    covar[distance >= maxdist] = 0
    covar[distance < maxdist] *= channel_correlation[distance[distance 
                                                              < maxdist]]
    return(covar)    

def calculate_channel_correlation(cube, length=1):
    raise NotImplementedError


def write_moment1_advanced(cube,
                           broad_mask=None,
                           moment1_prior=None,
                           velocity_convention='radio',
                           order='bilinear',
                           outfile=None,
                           errorfile=None,
                           rms=None,
                           channel_correlation=None,
                           overwrite=True,
                           vfield_reject_thresh=30 * u.km / u.s,
                           mom0_thresh_for_mom1=2.0,
                           unit=None):
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
    
    """

    mom1strict = cube.moment1()
    if unit is not None:
        mom1strict = mom1strict.to(unit)
        spaxis = cube.spectral_axis.to(unit).value
    else:
        spaxis = cube.spectral_axis.value

    if moment1_prior is not None:
        if type(moment1_prior) is Projection:
            mom1prior = moment1_prior
        elif type(moment1_prior) is str:
            hdu_list = fits.open(moment1_prior)
            mom1prior = Projection.from_hdu(hdulist[0])
        if unit is not None:
            mom1prior.to(unit)
        mom1prior = mom1prior.reproject(mom1strict.header, order=order)
    else:
        mom1prior = None

    
    if type(broad_mask) is SpectralCube:
        broad_cube = cube.with_mask(
            broad_mask.filled_data[:].value.astype(np.bool),
            inherit_mask=False)
    elif type(broad_mask) is np.ndarray:
        broad_cube = cube.with_mask(broad_mask.astype(np.bool),
                                    inherit_mask=False)
        
    mom0broad = broad_cube.moment0()
    mom1broad = broad_cube.moment1()
    mom0err = None
    mom1err = None
    # Calculate errors in mom0broad, mom1broad
    if rms is not None and errorfile is not None:
        if channel_correlation is None:
            channel_correlation = np.array([1])

        mom0err = np.empty(mom0broad.shape)
        mom0err.fill(np.nan)
        mom1err = np.empty(mom1broad.shape)
        mom1err.fill(np.nan)

        for x, y, slc in cube._iter_rays(0):
            mask = np.squeeze(cube._mask.include(view=slc))
            if not mask.any():
                continue
            index = np.where(mask)[0]
            rms_spec = rms.flattened(slc).value
            spec = cube.flattened(slc).value

            covar = build_covariance(spectrum=spec,
                                     rms=rms_spec,
                                     channel_correlation=channel_correlation,
                                     index=index)
            vval = spaxis[index]
            sum_T = np.sum(spec)
            sum_vT = np.sum(spec * vval)

            jacobian = vval / sum_T - sum_vT / sum_T**2
            mom1err[x, y] = np.dot(
                np.dot(jacobian[np.newaxis, :], covar),
                jacobian[:, np.newaxis])**0.5
            mom0err = (np.sum(covar**2))**0.5
        mom0err = u.Quantity(mom0err, mom0broad.unit, copy=False)
        mom1err = u.Quantity(mom1err, mom1broad.unit, copy=False)
        if unit is not None:
            mom1err = mom1err.to(unit)
    
    mom1hybrid = mom1strict.value
    valid_broad_mom1 = np.isfinite(mom1broad.value)
    valid_broad_mom1[np.isfinite(mom1strict)] = False
    
    if mom0err is not None:
        valid_broad_mom1 = (mom0.value 
                            > (mom0_thresh_for_mom1 * mom0err.value))
    if mom1prior is not None:
        valid_broad_mom1 = (valid_broad_mom1 * 
                            (np.abs(mom1broad - mom1prior) 
                             < vfield_reject_thresh)
                            )
    mom1hybrid[valid_broad_mom1] = (mom1broad.value)[valid_broad_mom1]
    return(mom1hybrid)
        # mom0err_proj = Projection(mom0err,
        #                           wcs=mom0_broad.wcs,
        #                           header=mom0_broad.header)
        # mom1err_proj = Projection(mom1err,
        #                           wcs=mom1_broad.wcs,
        #                           header=mom1_broad.header)
        
