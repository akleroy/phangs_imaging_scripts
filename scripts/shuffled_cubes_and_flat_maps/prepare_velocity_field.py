######################################################################################################
# PHANGS-ALMA (DR5) shuffled cubes pipeline, step 0: prepare velocity field
# Date: Jan 23, 2025
# Authors: Lukas Neumann, lukas.neumann@eso.org
######################################################################################################

# libraries
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from astropy.convolution import convolve, Gaussian2DKernel
from reproject import reproject_interp  # reprojection to target coordinate system
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import warnings
warnings.filterwarnings("ignore")

######################################################################################################
# CONFIG
######################################################################################################

# directories
wd = '/Users/lneumann/Documents/'  # working directory

# data directories
co_cubes_dir = wd + 'Data/PHANGS/alma/cubes/'  # PHANGS-ALMA CO cubes: native resolution cubes (to be shuffled)
co_maps_dir = wd + 'Data/PHANGS/alma/maps/'  # PHANGS-ALMA CO maps:CO moment-1 maps + err (native and 15" res)
ha_maps_dir = wd + 'Data/PHANGS/muse/maps/' # PHANGS-MUSE Halpha maps: if availabe, Halpha moment-1 maps + err (native and 15" res)
hi_maps_dir = wd + 'Data/PHANGS/hi/maps/'  # if availabe,  HI moment-1 maps + err
model_maps_dir = wd + 'Data/PHANGS/velocity_fields/smoothvfieldmodels2024/'  # if available, modelled velocity field

# diretory to save velocity field
vel_maps_dir = wd + 'Data/PHANGS/alma/shuffled_cubes/velocity_field/'  # combined velocity fields

# diagnostic plots
plot_dir = wd + 'Products/PHANGS/shuffled_cubes/diagnostics/'

# PHANGS sample table (full path)
sample_table_dir = wd + 'Data/PHANGS/tables/'  # directory
fpath_sample_table = sample_table_dir + 'phangs_sample_table_v1p6.csv'  # path

######################################################################################################
# target lists

# phangs alma sample
galaxies_alma = ['circinus','ic1954','ic5273','ic5332','ngc0247','ngc0253','ngc0300','ngc0628','ngc0685','ngc1068',
                 'ngc1087','ngc1097','ngc1300','ngc1317','ngc1365','ngc1385','ngc1433','ngc1511','ngc1512','ngc1546',
                 'ngc1559','ngc1566','ngc1637','ngc1672','ngc1792','ngc1809','ngc2090','ngc2283','ngc2566','ngc2775',
                 'ngc2835','ngc2903','ngc2997','ngc3059','ngc3137','ngc3239','ngc3351','ngc3489','ngc3507','ngc3511',
                 'ngc3521','ngc3596','ngc3599','ngc3621','ngc3626','ngc3627','ngc4207','ngc4254','ngc4293','ngc4298',
                 'ngc4303','ngc4321','ngc4424','ngc4457','ngc4459','ngc4476','ngc4477','ngc4496a','ngc4535','ngc4536',
                 'ngc4540','ngc4548','ngc4569','ngc4571','ngc4579','ngc4596','ngc4654','ngc4689','ngc4694','ngc4731',
                 'ngc4781','ngc4826','ngc4941','ngc4945','ngc4951','ngc5042','ngc5068','ngc5128','ngc5134','ngc5236',
                 'ngc5248','ngc5530','ngc5643','ngc6300','ngc6744','ngc7456','ngc7496','ngc7743','ngc7793']

# galaxies without 12m coverage
galaxies_alma_7m = ['circinus', 'ngc0247', 'ngc0253', 'ngc0300', 'ngc1068', 'ngc4945', 'ngc5128', 'ngc7793']

# phangs alma sample
# comment: note that these targets are not in the PHANGS sample table (v1.6) which is required for the diagnostics plots 
galaxies_alma_postlp = ['ngc1187', 'ngc1808', 'ngc3673', 'ngc4123', 'ngc4981', 'ngc5101']

# phangs muse sample
galaxies_muse = ['ic5332', 'ngc0628', 'ngc1087', 'ngc1300', 'ngc1365', 
                 'ngc1385', 'ngc1433', 'ngc1512', 'ngc1566', 'ngc1672', 
                 'ngc2835', 'ngc3351', 'ngc3627', 'ngc4254', 'ngc4303', 
                 'ngc4321', 'ngc4535', 'ngc5068', 'ngc7496']
muse_res_list = [0.87, 0.92, 0.92, 0.89, 1.15, 
                 0.77, 0.91, 1.25, 0.80, 0.96, 
                 1.15, 1.05, 1.05, 0.89, 0.78, 
                 1.16, 0.56, 1.04, 0.89]

# phangs muse sample extension
galaxies_muse_extended = ['circinus','ic1954','ic5273','ngc1097','ngc1317',
                          'ngc1808','ngc2775','ngc2903','ngc3239','ngc3489',
                          'ngc3521','ngc3596','ngc3599','ngc3626','ngc4424',
                          'ngc4457','ngc4496a','ngc4548','ngc4596','ngc4694',
                          'ngc4731','ngc4781','ngc4941','ngc4945','ngc5248',
                          'ngc7456','ngc7743']

######################################################################################################

# select sample
galaxies = galaxies_alma

# toggle diagnostic plot
diagnostics = True


######################################################################################################
# PIPELINE (0): VELOCITY FIELD
######################################################################################################

print('[INFO] Step 0:\tInitialise velocity field script.')

for galaxy in galaxies:

    print(f'[INFO] {galaxy.upper()}: Start producing combined velocity field.')

    # ALMA array and sample
    if galaxy in galaxies_alma_7m:
        array = '7m+tp'
    else:
        array = '12m+7m+tp'
        
    if galaxy in galaxies_alma_postlp:
        line = 'co21lores'
    else:
        line = 'co21'

    # load PHANGS table
    phangs_table = pd.read_csv(fpath_sample_table, comment='#', skiprows=1)
    
    # get systemic velocity of the galaxy
    if galaxy == 'circinus':
        # comment: circinus has an alternative name in sample table
        phangs_table_galaxy = phangs_table[phangs_table['name'] == 'eso097-013']
    else:
        phangs_table_galaxy = phangs_table[phangs_table['name'] == galaxy]
    v_lsr = float(phangs_table_galaxy['orient_vlsr'].iloc[0])  # systemic velocity [km/s]
        
    ###################################################
    # ALMA CO
    ###################################################

    # load phangs-alma CO mom1 (12m+7m+tp or 12m+tp = native res)
    fname_co_mom1 = co_maps_dir + f'{galaxy}_{array}_{line}_strict_mom1.fits'
    co_mom1, hdr_co_mom1 = fits.getdata(fname_co_mom1, header=True)
    wcs_co = WCS(hdr_co_mom1)
    
    # load phangs-alma CO mom1 (7m+tp = 15" res)
    fname_co_mom1_15as = co_maps_dir + f'{galaxy}_7m+tp_{line}_15as_strict_mom1.fits'
    co_mom1_15as, hdr_co_mom1_15as = fits.getdata(fname_co_mom1_15as, header=True)
    
    # reproject to PHANGS-ALMA CO grid
    co_mom1_15as_repr, _ = reproject_interp((co_mom1_15as, hdr_co_mom1_15as), hdr_co_mom1)
    
    ###################################################
    # MUSE Halpha
    ###################################################

    # load phangs-muse at native res
    if galaxy in galaxies_muse:
        # MUSE maps filename (multiextension file)
        res_muse = muse_res_list[galaxies_muse.index(galaxy)]
        fname_halpha_mom1 = ha_maps_dir + f'{galaxy.upper()}-{res_muse}asec_MAPS.fits'
        # extract Halpha mom1 map
        hdu_halpha = fits.open(fname_halpha_mom1)
        halpha_mom1 = hdu_halpha[97].data # Halpha mom1
        hdr_halpha_mom1 = hdu_halpha[97].header # header
        halpha_mom1_err = hdu_halpha[98].data # error
        halpha_mom1[halpha_mom1_err > 10] = np.nan # clip low significant data
        halpha = True # activate Halpha key
        
    # load phangs-muse at native res (extended sample)
    elif galaxy in galaxies_muse_extended:
        # MUSE maps filename (multiextension file)
        res_muse = 1  # TBD: dig up MUSE resolution
        fname_halpha_mom1 = ha_maps_dir + f'{galaxy.upper()}_MAPS.fits'
        # extract Halpha mom1 map
        hdu_halpha = fits.open(fname_halpha_mom1)
        halpha_mom1 = hdu_halpha[95].data # Halpha mom1
        hdr_halpha_mom1 = hdu_halpha[95].header # header
        halpha_mom1_err = hdu_halpha[96].data # error
        halpha_mom1[halpha_mom1_err > 10] = np.nan # clip low significant data
        halpha = True # activate Halpha key

    else:
        # Halpha data NOT available
        print(f'[INFO] {galaxy.upper()}: MUSE Halpha data not available!')
        halpha_mom1_repr = np.full_like(co_mom1, np.nan) # dummy map for velocity field combination
        halpha = False # deactivate Halpha key

    # convolve and reproject Halpha map to PHANGS-ALMA CO grid
    if halpha:
        # convolve to PHANGS-ALMA CO resolution
        res_alma = hdr_co_mom1['BMAJ'] * 3600
        # check that MUSE resolution is higher than CO resolution
        if res_muse < res_alma:
            # build convolution kernel
            try:
                deg_per_px = abs(hdr_halpha_mom1['CD1_1']) 
            except:
                deg_per_px = abs(hdr_halpha_mom1['PC1_1']) 
            as_per_px = deg_per_px * 3600
            current_res_px = res_muse / as_per_px
            target_res_px = res_alma / as_per_px
            kernel_fwhm_px = np.sqrt(target_res_px**2 - current_res_px**2)
            kernel_std_px = kernel_fwhm_px / (2*np.sqrt(2*np.log(2)))
            kernel = Gaussian2DKernel(kernel_std_px)
        
            # convolve to target resolution
            halpha_mom1 =  convolve(halpha_mom1, kernel, boundary='extend', preserve_nan=True)
            halpha_mom1[halpha_mom1==0] = np.nan
        else:
            print(f'[INFO] {galaxy.upper()}: MUSE Halpha resolution is coarser than CO resolution!')
        
        # reproject to PHANGS-ALMA CO grid
        halpha_mom1_repr, _ = reproject_interp((halpha_mom1, hdr_halpha_mom1), hdr_co_mom1)
        
        # correct for velocity offset
        # comment: there is an offset between the LOS velocities of ALMA and MUSE observations
        # comment: this offsets comes from different definitions of the LSR velocities
        # comment: here, we apply a brute-force subtraction of the offset assuming that the medians agree
        halpha_mom1_repr -= np.nanmedian(halpha_mom1_repr - co_mom1)

        
    ###################################################
    # load phangs-muse at 15" res
    # comment: 15" maps do not exist for the extended sample (yet)
    if galaxy in galaxies_muse:
        # MUSE maps filename (multiextension file)
        fname_halpha_mom1_15as = ha_maps_dir + f'{galaxy.upper()}-15asec_MAPS.fits'
        # extract Halpha mom1 map
        hdu_halpha_15as = fits.open(fname_halpha_mom1_15as)
        halpha_mom1_15as = hdu_halpha_15as[97].data # Halpha mom1
        hdr_halpha_mom1_15as = hdu_halpha_15as[97].header # header
        halpha_mom1_15as_err = hdu_halpha_15as[98].data # error
        halpha_mom1_15as[halpha_mom1_15as_err > 2] = np.nan
        
        # reproject to PHANGS-ALMA CO grid
        halpha_mom1_15as_repr, _ = reproject_interp((halpha_mom1_15as, hdr_halpha_mom1_15as), hdr_co_mom1)
        
        # correct for offset
        halpha_mom1_15as_repr -= np.nanmedian(halpha_mom1_15as_repr - co_mom1)
    
        # Halpha data available
        halpha_15as = True # activate Halpha key
    else:
        print(f'[INFO] {galaxy.upper()}: MUSE Halpha data @ 15 arcsec not available!')
        # Halpha data NOT available
        halpha_mom1_15as_repr = np.full_like(co_mom1, np.nan) # dummy map for velocity field combination
        halpha_15as = False # deactivate Halpha key


    ###################################################
    # HI
    ###################################################
    try:
        # load HI mom1 map
        fname_hi_mom1 = hi_maps_dir + f'{galaxy.upper()}_21cm_strictmask_mom1.fits'
        hi_mom1, hdr_hi_mom1 = fits.getdata(fname_hi_mom1, header=True)
        
        # reproject to PHANGS-ALMA CO grid
        hi_mom1_repr, _ = reproject_interp((hi_mom1, hdr_hi_mom1), hdr_co_mom1)

        if galaxy == 'ngc4535':
            hi_mom1_repr -= np.nanmedian(hi_mom1_repr - co_mom1_15as_repr) # correct for offset (specific to this galaxy)
        
        # HI data available
        hi = True # activate HI key
        
    except:
        # HI data NOT available
        print(f'[INFO] {galaxy.upper()}: HI 21cm data not available!')
        hi_mom1_repr = np.full_like(co_mom1, np.nan) # dummy map for velocity field combination
        hi = False # deactivate HI key

    
    ###################################################
    # Smoothed rotation curves (Lang, Meid et al. 2020)
    ###################################################
    try:
        # load rotation curve
        fname_model = model_maps_dir + f'{galaxy}legendreLOSmodelFULL.fits'
        model_mom1, hdr_model = fits.getdata(fname_model, header=True)
        
        # reproject to PHANGS-ALMA CO grid
        model_mom1_repr, _ = reproject_interp((model_mom1, hdr_model), hdr_co_mom1)

        # add systemic velocity
        model_mom1_repr += v_lsr
        
        # Model available
        model = True # activate model key
        
    except:
        # Model NOT available
        print(f'[INFO] {galaxy.upper()}: Model rotation curve not available!')
        model_mom1_repr = np.full_like(co_mom1, np.nan) # dummy map for velocity field combination
        model = False # deactivate model key
        
    
    ###################################################
    # Combined velocity field
    ###################################################
    
    # start with CO map
    vel_field = np.copy(co_mom1)  # copy CO mom1
    vel_field[np.isnan(vel_field)] = 0
    tracers = 'CO'
    
    # add Halpha map
    if halpha:
        mom1_temp = np.copy(halpha_mom1_repr)
        mom1_temp[vel_field != 0] = 0
        vel_field += mom1_temp
        vel_field[np.isnan(vel_field)] = 0
        tracers += '+Ha'
    
    # add smoothed CO map
    mom1_temp = np.copy(co_mom1_15as_repr)
    mom1_temp[vel_field != 0] = 0
    vel_field += mom1_temp
    vel_field[np.isnan(vel_field)] = 0

    # add smoothed Halpha map
    if halpha_15as:
        mom1_temp = np.copy(halpha_mom1_15as_repr)
        mom1_temp[vel_field != 0] = 0
        vel_field += mom1_temp
        vel_field[np.isnan(vel_field)] = 0
    
    # add HI map
    if hi:
        mom1_temp = np.copy(hi_mom1_repr)
        mom1_temp[vel_field != 0] = 0
        vel_field += mom1_temp
        vel_field[np.isnan(vel_field)] = 0
        tracers += '+HI'

    # add model map
    if model:
        mom1_temp = np.copy(model_mom1_repr)
        mom1_temp[vel_field != 0] = 0
        vel_field += mom1_temp
        tracers += '+model'
        
    vel_field[vel_field == 0] = np.nan
    
    # save as fits file
    hdr_co_mom1['OBJECT'] = galaxy
    hdr_co_mom1['LINES'] = tracers
    fits.writeto(vel_maps_dir + f'{galaxy}_combined_velocity_field.fits', vel_field, hdr_co_mom1, overwrite=True)

    print(f'[INFO] {galaxy.upper()}: Velocity field includes: {tracers}.')
    print(f'[INFO] {galaxy.upper()}: Velocity field; successfully finished.')

    ######################################################################################################

    # diagnostics
    if diagnostics:

        ###################################################
        # Filling factors
        ###################################################
        
        # compute CO filling factor (native res)
        # load CO mom0 map (used to compute number of pixels in FOV)
        fname_co_mom0 = co_maps_dir + f'{galaxy}_{array}_{line}_broad_mom0.fits'
        co_mom0 = fits.getdata(fname_co_mom0, header=False)
        co_cutout = np.copy(co_mom1)  # copy CO mom1 map
        co_cutout[np.isnan(co_mom0)] = np.nan
        co_filling = np.sum(~np.isnan(co_cutout.flatten())) / np.sum(~np.isnan(co_mom0.flatten()))  

        # compute Halpha filling factor (native res + 15 arcsecond) 
        if halpha:
            
            halpha_cutout = np.copy(halpha_mom1_repr)
            halpha_cutout[np.isnan(co_mom0)] = np.nan
            halpha_filling = np.sum(~np.isnan(halpha_cutout.flatten())) / np.sum(~np.isnan(co_mom0.flatten()))
            halpha_15as_cutout = np.copy(halpha_mom1_15as_repr)
            halpha_15as_cutout[np.isnan(co_mom0)] = np.nan
            halpha_15as_filling = np.sum(~np.isnan(halpha_15as_cutout.flatten())) / np.sum(~np.isnan(co_mom0.flatten()))
        else:
            halpha_filling = np.nan
            halpha_15as_filling = np.nan
            
        # compute CO filling factor (15 arcsecond)
        co_15as_cutout = np.copy(co_mom1_15as_repr)
        co_15as_cutout[np.isnan(co_mom0)] = np.nan
        co_15as_filling = np.sum(~np.isnan(co_15as_cutout.flatten())) / np.sum(~np.isnan(co_mom0.flatten()))
        
        # compute HI filling factor
        if hi:
            hi_cutout = np.copy(hi_mom1_repr)
            hi_cutout[np.isnan(co_mom0)] = np.nan
            hi_filling = np.sum(~np.isnan(hi_cutout.flatten())) / np.sum(~np.isnan(co_mom0.flatten()))
        else:
            hi_filling = np.nan
        
        # compute model filling factor
        if model:
            model_cutout = np.copy(model_mom1_repr)
            model_cutout[np.isnan(co_mom0)] = np.nan
            model_filling = np.sum(~np.isnan(model_cutout.flatten())) / np.sum(~np.isnan(co_mom0.flatten()))
        else:
            model_filling = np.nan


        ###################################################
        # PLOTTING
        ###################################################

        # create figure
        fig = plt.figure(figsize=(15, 14))
        fig.subplots_adjust(wspace=0.3, hspace=0.3)

        # matplotlib settings
        plt.rcParams['xtick.direction'] = 'in'
        plt.rcParams['ytick.direction'] = 'in'
        
        # figure dimensions
        nrows = 6
        ncols = 6
        
        # maps to be plotted
        data_list = [co_mom1, halpha_mom1_repr, co_mom1_15as_repr, halpha_mom1_15as_repr, hi_mom1_repr, model_mom1_repr]
        labels = ['CO', r'H$\alpha$', r'CO$_{15^{\prime\prime}}$', r'H$\alpha_{15^{\prime\prime}}$', 'HI', 'Model']
        filling_factors = [co_filling, halpha_filling, co_15as_filling, halpha_15as_filling, hi_filling, model_filling]
        map_indeces = [1,8,15,22,29,36]

        # velocity range (for all panels)
        vmin = np.nanmin(co_mom1)
        vmax= np.nanmax(co_mom1)
        
        # iterate over maps
        for data, idx, label in zip(data_list, map_indeces, labels):
            ax = plt.subplot(nrows, ncols, idx, projection=wcs_co) # create subplot
            im = ax.imshow(data, transform=ax.get_transform(wcs_co), cmap='rainbow', vmin=vmin, vmax=vmax)
            cax = ax.inset_axes([0, 1, 1, 0.05], transform=ax.transAxes, zorder=5)
            cbar = plt.colorbar(im, cax=cax, orientation='horizontal')
            cax.xaxis.set_ticks_position('top')
            ax.set_title(label+' vel. [km/s]', size=14, pad=10)
            # axis labels
            ax.set_xlabel('R.A. [J2000]')
            ax.set_ylabel('Dec. [J2000]', labelpad=-0.1)
            ax.tick_params(labelsize=6)
            edge_px=50
            ax.set_xlim(0-edge_px, hdr_co_mom1['NAXIS1']+edge_px)
            ax.set_ylim(0-edge_px, hdr_co_mom1['NAXIS2']+edge_px)
            
            # FOV contour
            data_fov, header_fov = fits.getdata(fname_co_mom0, header=True)
            frame_width = 100
            data_fov = np.vstack([np.full([frame_width, np.shape(data_fov)[1]], np.nan), data_fov, np.full([frame_width, np.shape(data_fov)[1]], np.nan)]) # expand 2D array with a frame of nans
            data_fov = np.vstack([np.full([frame_width, np.shape(data_fov)[0]], np.nan), data_fov.T, np.full([frame_width, np.shape(data_fov)[0]], np.nan)]).T
            header_fov['CRPIX1'] += frame_width
            header_fov['CRPIX2'] += frame_width
            wcs_fov = WCS(header_fov)
            data_fov[np.isnan(data_fov)] = -1
            data_fov[data_fov != -1] = 1
            # roll = 12
            # mask = (np.roll(data_fov,roll,axis=0)==1) | (np.roll(data_fov,-roll,axis=0)==1) | (np.roll(data_fov,roll,axis=1)==1) | (np.roll(data_fov,-roll,axis=1)==1)
            # data_fov[mask] = 1
            # mask = (np.roll(data_fov,roll,axis=0)==1) & (np.roll(data_fov,-roll,axis=0)==1) & (np.roll(data_fov,roll,axis=1)==1) & (np.roll(data_fov,-roll,axis=1)==1)
            # data_fov[~mask] = -1
            ax.contour(data_fov, levels=[0], colors='k', linewidths=1, linestyles='solid', transform=ax.get_transform(wcs_fov))
            
        
        ########################################
        # scatter plots
        xdata_list = [co_mom1, 
                      co_mom1, halpha_mom1_repr, 
                      co_mom1, halpha_mom1_repr, co_mom1_15as_repr,
                      co_mom1, halpha_mom1_repr, co_mom1_15as_repr, halpha_mom1_15as_repr,
                      co_mom1, halpha_mom1_repr, co_mom1_15as_repr, halpha_mom1_15as_repr, hi_mom1_repr]
        ydata_list = [halpha_mom1_repr, 
                      co_mom1_15as_repr, co_mom1_15as_repr, 
                      halpha_mom1_15as_repr, halpha_mom1_15as_repr, halpha_mom1_15as_repr, 
                      hi_mom1_repr, hi_mom1_repr, hi_mom1_repr, hi_mom1_repr,
                      model_mom1_repr, model_mom1_repr, model_mom1_repr, model_mom1_repr, model_mom1_repr]
        xlabels = ['CO', 
                   'CO', r'H$\alpha$', 
                   'CO', r'H$\alpha$', r'CO$_{15^{\prime\prime}}$',
                   'CO', r'H$\alpha$', r'CO$_{15^{\prime\prime}}$', r'H$\alpha_{15^{\prime\prime}}$',
                   'CO', r'H$\alpha$', r'CO$_{15^{\prime\prime}}$', r'H$\alpha_{15^{\prime\prime}}$', 'Model']
        ylabels = [r'H$\alpha$', 
                   r'CO$_{15^{\prime\prime}}$', r'CO$_{15^{\prime\prime}}$', 
                   r'H$\alpha_{15^{\prime\prime}}$', r'H$\alpha_{15^{\prime\prime}}$', r'H$\alpha_{15^{\prime\prime}}$',
                   'HI', 'HI', 'HI', 'HI',
                   'Model', 'Model', 'Model', 'Model', 'Model']
        
        # plot indeces
        scatter_indeces = [7, 13,14, 19,20,21, 25,26,27,28, 31,32,33,34,35]
        
        # iterate over 
        for xdata, ydata, idx, xlabel, ylabel in zip(xdata_list, ydata_list, scatter_indeces, xlabels, ylabels):
            ax = plt.subplot(nrows, ncols, idx)  # create subplot
            ax.scatter(xdata, ydata, c='darkblue', s=2, zorder=1)  # plot data
            ax.scatter(xdata, ydata, c='dodgerblue', s=1, zorder=1)
            # axis labels
            if idx > ncols*(nrows-1):
                ax.set_xlabel(xlabel+' vel. [km/s]', size=14)
            else:
                ax.set_xlabel(' ')
            if idx % ncols == 1:    
                ax.set_ylabel(ylabel+' vel. [km/s]', size=14)
            else:
                ax.set_ylabel(' ')
            # plot range
            ax.set_xlim(vmin-20, vmax+20)
            ax.set_ylim(vmin-20, vmax+20)
            # 1-1 relation
            ax.plot(np.linspace(vmin-20, vmax+20, 10), np.linspace(vmin-20, vmax+20, 10), lw=2, ls='-', c='k', zorder=5)
            # 20 km/s deviation
            ax.plot(np.linspace(vmin-20, vmax+20, 10), np.linspace(vmin-20, vmax+20, 10)+20, ls='dashed', c='k', zorder=5)
            ax.plot(np.linspace(vmin-20, vmax+20, 10), np.linspace(vmin-20, vmax+20, 10)-20, ls='dashed', c='k', zorder=5)
            # 50 km/s deviation
            ax.plot(np.linspace(vmin-20, vmax+20, 10), np.linspace(vmin-20, vmax+20, 10)+50, ls='dotted', c='k', zorder=5)
            ax.plot(np.linspace(vmin-20, vmax+20, 10), np.linspace(vmin-20, vmax+20, 10)-50, ls='dotted', c='k', zorder=5)
        
        
        ########################################
        # combined velocity field
        ax = plt.subplot(nrows, ncols, ncols, projection=wcs_co) # create subplot
        im = ax.imshow(vel_field, transform=ax.get_transform(wcs_co), cmap='rainbow', vmin=vmin, vmax=vmax)
        cax = ax.inset_axes([0, 1, 1, 0.05], transform=ax.transAxes, zorder=5)
        cbar = plt.colorbar(im, cax=cax, orientation='horizontal')
        cax.xaxis.set_ticks_position('top')
        ax.set_title('Combined vel. [km/s]', size=16, pad=10)
        # axis labels
        ax.set_xlabel('R.A. [J2000]')
        ax.set_ylabel('Dec. [J2000]', labelpad=-0.1)
        ax.tick_params(labelsize=6)
        # ax.text(0.95, 0.95, 'combined', ha='right', va='top', transform=ax.transAxes, bbox=box, zorder=10, size=16)
        ax.set_xlim(0-edge_px, hdr_co_mom1['NAXIS1']+edge_px)
        ax.set_ylim(0-edge_px, hdr_co_mom1['NAXIS2']+edge_px)
        
        # FOV contour
        ax.contour(data_fov, levels=[0], colors='k', linewidths=1, linestyles='solid', transform=ax.get_transform(wcs_fov))
        

        # arrows
        transFigure = fig.transFigure.inverted()
        coord1_list = [[0.07,0.1], [0.25, 0.81], [0.38, 0.7], [0.52, 0.58], [0.63, 0.53], [0.72, 0.4], [0.85, 0.28]]
        coord2_list = [[0.07,1]] + [[0.85, 0.81]]*6

        if not model:
            del coord1_list[6]
            del coord2_list[6]
        
        if not hi:
            del coord1_list[5]
            del coord2_list[5]
        
        if not halpha_15as:
            del coord1_list[4]
            del coord2_list[4]
            
        if not halpha:
            del coord1_list[2]
            del coord2_list[2]

        for coord1, coord2 in zip(coord1_list, coord2_list):
            arrow = mpatches.FancyArrowPatch(
                coord1,  # posA
                coord2,  # posB
                shrinkA=0,  # so tail is exactly on posA (default shrink is 2)
                shrinkB=120,  # so head is exactly on posB (default shrink is 2)
                transform=fig.transFigure,
                color="black",
                arrowstyle="-|>",  # "normal" arrow
                mutation_scale=30,  # controls arrow head size
                linewidth=3)
            fig.patches.append(arrow)
        
        # table with filling factors
        textsize = 14
        y_off = 0.015
        
        # show filling factors
        plt.text(0.46, 0.93+0.25*y_off, 'data', ha='center', va='top', transform=fig.transFigure, size=textsize)
        plt.text(0.54, 0.93+0.25*y_off, 'filling factor', ha='center', va='top', transform=fig.transFigure, size=textsize)
        for label, filling_factor in zip(labels, filling_factors):
            idx_label = labels.index(label)
            plt.text(0.46, 0.93-(idx_label+1)*y_off, label, ha='center', va='top', transform=fig.transFigure, size=textsize)
            plt.text(0.54, 0.93-(idx_label+1)*y_off, '%.2f' % filling_factor, ha='center', va='top', transform=fig.transFigure, size=textsize)
        
        ########################################
        # figure title
        plt.text(0.05, 0.5, 'Priority', ha='center', va='center', transform=fig.transFigure, size=24, rotation=90)
        fig.suptitle(galaxy.upper(), size=30)
        
        # save figure
        plotname = f'{galaxy}_velocity_field'
        plt.savefig(plot_dir + plotname + '.png', facecolor='white', transparent=False, bbox_inches='tight')

print('[INFO] Step 0:\tSuccessfully finished.')