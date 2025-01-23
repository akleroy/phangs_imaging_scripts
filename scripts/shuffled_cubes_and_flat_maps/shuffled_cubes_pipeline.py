######################################################################################################
# PHANGS-ALMA (DR5) shuffled cubes pipeline
# Date: Jan 23, 2025
# Authors: Lukas Neumann, lukas.neumann@eso.org
######################################################################################################

# libraries
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm  # for colour stretch
from matplotlib.ticker import ScalarFormatter  # ticklabel format
import matplotlib.gridspec as gridspec
from ancillary_functions import shuffle_cube, get_pv_data
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
vel_maps_dir = wd + 'Data/PHANGS/alma/shuffled_cubes/velocity_field/'  # combined velocity fields

# pipeline products
shuffled_cubes_dir = wd + 'Data/PHANGS/alma/shuffled_cubes/cubes/'  # shuffled cubes
flat_maps_dir = wd + 'Data/PHANGS/alma/shuffled_cubes/maps/'  # flat maps

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
                          'ngc2775','ngc2903','ngc3239','ngc3489','ngc3521',
                          'ngc3596','ngc3599','ngc3626','ngc4424','ngc4457',
                          'ngc4496a','ngc4548','ngc4596','ngc4694','ngc4731',
                          'ngc4781','ngc4941','ngc4945','ngc5248','ngc7456',
                          'ngc7743']

######################################################################################################
# velocity windows
# fixed windows
v_width_50kms = 50 # km/s (plus/minus 25 km/s velocity range)
v_width_100kms = 100 # km/s (plus/minus 50 km/s velocity range)
# narrow window (adapted to each galaxy)
v_width_narrow_list_alma = [70, 30, 30, 20, 30, 40, 20, 40, 40, 40, 
                            40, 60, 30, 50, 200, 60, 30, 80, 30, 60, 
                            60, 40, 30, 40, 60, 40, 40, 40, 40, 30, 
                            30, 50, 50, 40, 40, 30, 30, 60, 40, 40, 
                            60, 40, 30, 60, 40, 60, 50, 60, 60, 50, 
                            60, 50, 50, 60, 60, 50, 40, 40, 50, 40, 
                            50, 40, 40, 30, 50, 50, 50, 40, 40, 50, 
                            50, 50, 30, 80, 40, 40, 30, 50, 50, 50, 
                            50, 40, 50, 50, 20, 30, 30, 50, 30]
# wide window (fixed across galaxies)
v_width_wide = 100 # km/s

######################################################################################################

# select sample
galaxies = galaxies_alma
v_width_narrow_list = v_width_narrow_list_alma

# toggle pipeline step
run_shuffling = True
run_flatmaps = True

# toggle diagnostic plot
diagnostics = True


######################################################################################################
# PIPELINE
######################################################################################################

print('[INFO] Start pipeline.')

if run_shuffling:

    ######################################################################################################
    # PIPELINE (1): SHUFFLING
    ######################################################################################################

    print('[INFO] Step 1:\tInitialise shuffled cubes script.')

    # list of cubes to be shuffled
    cube_files = ['.fits', '_noise.fits', '_broadmask.fits', '_strictmask.fits']

    # loop over galaxies
    for galaxy in galaxies:

        print(f'[INFO] {galaxy.upper()}: Start producing shuffled cubes.')

        # loop over cube files
        for cube_file in cube_files:

            # ALMA array
            if galaxy in galaxies_alma_7m:
                array = '7m+tp'
            else:
                array = '12m+7m+tp'

            # ALMA line (accounting for post-LP sample)
            if galaxy in galaxies_alma_postlp:
                line = 'co21lores'
            else:
                line = 'co21'
            
            ###################################################
            # LOAD DATA
            ###################################################

            # load cube to be shuffled (DR4)
            fname_cube = co_cubes_dir + f'{galaxy}_{array}_{line}{cube_file}'
            cube, hdr_cube = fits.getdata(fname_cube, header=True)
            
            # load cube header
            fname_co_cube = co_cubes_dir + f'{galaxy}_{array}_{line}.fits'
            hdr_co_cube = fits.getheader(fname_co_cube)
            
            # load velocity field
            fname_vfield = vel_maps_dir + f'{galaxy}_combined_velocity_field.fits'
            vfield, hdr_vfield = fits.getdata(fname_vfield, header=True)


            ###################################################
            # MASKING (only for NGC 6744)
            ###################################################

            if galaxy == 'ngc6744':

                # trim edges using loaded mask
                fname_mask = co_maps_dir + f'{galaxy}_{array}_{line}_mask_edges.fits'
                co_mask_2d = fits.getdata(fname_mask)

                # create cube mask
                co_mask_3d = np.copy(cube)
                for i in range(np.shape(co_mask_3d)[0]):
                    co_mask_3d[i,:,:] = co_mask_2d

                # apply mask to cube
                if cube_file in ['_broadmask.fits', '_strictmask.fits']:
                    cube[co_mask_3d == 0] = 0
                else:
                    cube[co_mask_3d == 0] = np.nan 
            
            ###################################################
            # SHUFFLING
            ###################################################
            
            # get velocity axis
            n_ch = hdr_co_cube['NAXIS3'] # number of channels
            v_ch = hdr_co_cube['CDELT3'] # channel width (m/s, with sign)
            v_ref = hdr_co_cube['CRVAL3'] # reference channel (m/s)
            ch_ref = hdr_co_cube['CRPIX3'] # reference channel
            v_ch0 = v_ref - (ch_ref-1)*v_ch # velocity of first channel
            vaxis = np.linspace(v_ch0, v_ch0+(n_ch-1)*v_ch, n_ch)  # velocity axis
            if hdr_co_cube['CUNIT3'] == 'm/s':
                vaxis *= 1e-3  # convert to km/s
            else:
                print('[WARNING] Velocity axis units not specified, assuming km/s')
            
            # shuffle cube
            cube_shuffled = shuffle_cube(cube, vaxis, vfield)
        
            # edit header (shuffle velocity axis)
            hdr_cube['CRVAL3'] = 0
            hdr_cube['CRPIX3'] = hdr_cube['NAXIS3']//2 + 1  
            
            # save shuffled cube
            fname_cube_shuffled = shuffled_cubes_dir + f'{galaxy}_{array}_{line}_shuffled{cube_file}'
            fits.writeto(fname_cube_shuffled, cube_shuffled, hdr_cube, overwrite=True)

        ######################################################################################################

        if diagnostics:

            # load PHANGS table
            phangs_table = pd.read_csv(fpath_sample_table, comment='#', skiprows=1)
            
            # get position angle and centre coordinates of galaxy
            if galaxy == 'circinus':
                phangs_table_galaxy = phangs_table[phangs_table['name'] == 'eso097-013']  # select row in table
            else:
                phangs_table_galaxy = phangs_table[phangs_table['name'] == galaxy]  # select row in table
            ctr_ra = float(phangs_table_galaxy['orient_ra'].iloc[0])  # get galaxy centre right ascension
            ctr_dec = float(phangs_table_galaxy['orient_dec'].iloc[0])  # get galaxy centre declination
            posang = float(phangs_table_galaxy['orient_posang'].iloc[0])  # get position angle

            # filenames
            co_cube_fname = co_cubes_dir + f'{galaxy}_{array}_{line}.fits'
            co_cube_shuffled_fname = shuffled_cubes_dir + f'{galaxy}_{array}_{line}_shuffled.fits'
            fname_list = [co_cube_fname, co_cube_shuffled_fname]
            labels = ['original cube', 'shuffled cube']
            
            # labels
            box = dict(boxstyle='round', facecolor='k', linewidth=0.5, alpha=0.8, pad=0.2)  # set box
            
            # create figure
            fig = plt.figure(figsize=(12, 3))
            fig.subplots_adjust(wspace=0.3)
            
            # iterate over cubes
            for fname, label in zip(fname_list, labels):
            
                # get pv-diagram data for CO cube
                pv_data, position, velocity = get_pv_data(fits_cube = fname, 
                                                        ctr_ra=ctr_ra, ctr_dec=ctr_dec, posang=posang, 
                                                        # bin_width_major=10/3600, # units are degree (so 10 arcseconds is 10/3600 degree, native resolution is around 1 arcsecond)
                                                        # bin_width_minor=100/3600, 
                                                        # bin_width_velocity=10,  # units are km/s (channel width of the data is 2.5 km/s)
                                                        # rgal_axis_length = 200/3600  # units are degree
                                                        # fits_cube_mask = data_dir + fname_cube_mask
                                                        )  
                
                ########################################
                
                # create subplot
                n = fname_list.index(fname)
                ax = plt.subplot(1, 3, n+1)
                
                # plot PV diagram
                vmax = np.nanmax(pv_data)
                im = ax.imshow(pv_data, # plot pv data as 2d array
                            cmap='gist_heat',  # colormap
                            norm=SymLogNorm(linthresh=np.nanpercentile(pv_data, 90), vmin=0, vmax=vmax), # color stretch
                            extent=[position[0]*3600, position[-1]*3600, velocity[0], velocity[-1]],  # set axis values according to position, velocity values
                            aspect='auto'  # stretch pixels according to figure size
                            )
                
                # colourbar
                vmin = 0.01
                cbar_ticks = np.logspace(np.log10(vmin), np.log10(vmax), 5)
                cbar_ticks = [0] + list(cbar_ticks[:-1]) + [vmax]
                cax = ax.inset_axes([0, 1, 1, 0.05], transform=ax.transAxes, zorder=5)
                cbar = plt.colorbar(im, cax=cax, orientation='horizontal', ticks=cbar_ticks)
                cax.xaxis.set_ticks_position('top')
                cbar.ax.set_title('$T_\mathrm{CO(2-1)}$ [K]', pad=10, size=12)
            
                # colourbar ticks
                formatter = ScalarFormatter()
                formatter.set_scientific(False)
                cax.xaxis.set_major_formatter(formatter)
                cbar.ax.tick_params(which='minor', top=False, bottom=False)
                cbar.ax.tick_params(which='minor', labeltop=False, labelbottom=False)    
                cbar.ax.set_xticklabels(['{:.2f}'.format(x) for x in cbar_ticks])
                cbar.ax.tick_params(labelsize=6)

                # axis labels
                ax.set_xlabel('Position along major axis [arcsecond]')
                ax.set_ylabel('Velocity [km/s]')
                ax.tick_params(axis='both', which='both', colors='white', labelcolor='k')
                ax.tick_params(labelsize=6)
            
                # indicate reference velocity
                if n == 1:
                    ax.axhline(0, ls='dashed', c='white')
                    
                # labels
                ax.text(0.97, 0.97, label, ha='right', va='top', transform=ax.transAxes, c='white', bbox=box, zorder=10, size=8)

            #############
            # CHANNEL MAP
            #############

            # get channel map at 0 km/s
            co_cube, hdr_co_cube = fits.getdata(co_cube_shuffled_fname, header=True)
            n_ch = hdr_co_cube['NAXIS3'] # number of channels
            v_ch = hdr_co_cube['CDELT3'] # channel width (m/s, with sign)
            v_ref = hdr_co_cube['CRVAL3'] # reference channel (m/s)
            ch_ref = hdr_co_cube['CRPIX3'] # reference channel
            v_ch0 = v_ref - (ch_ref-1)*v_ch # velocity of first channel
            vaxis = np.linspace(v_ch0, v_ch0+(n_ch-1)*v_ch, n_ch)
            id_ch0 = np.where(vaxis==0)[0][0]
            T_ch0 = co_cube[id_ch0]
            fname_co_mom0 = co_maps_dir + f'{galaxy}_{array}_{line}_strict_mom0.fits'
            hdr_T_ch0 = fits.getheader(fname_co_mom0)
            wcs = WCS(hdr_T_ch0)


            ########################################
            # combined velocity field
            ax = plt.subplot(1, 3, 3, projection=wcs)
            vmax = np.nanmax(T_ch0)
            im = ax.imshow(T_ch0, 
                        transform=ax.get_transform(wcs), cmap='gist_heat', 
                        norm=SymLogNorm(linthresh=np.nanpercentile(T_ch0, 90), vmin=0, vmax=vmax), # color stretch
                        )
            vmin = 0.1
            cbar_ticks = np.logspace(np.log10(vmin), np.log10(vmax), 5)
            cbar_ticks = [0] + list(cbar_ticks[:-1]) + [vmax]
            cax = ax.inset_axes([0, 1, 1, 0.05], transform=ax.transAxes, zorder=5)
            cbar = plt.colorbar(im, cax=cax, orientation='horizontal', ticks=cbar_ticks)
            cax.xaxis.set_ticks_position('top')
            cbar.ax.set_title('$T_\mathrm{CO(2-1)}$ [K]', pad=10, size=12)
            cbar.ax.tick_params(labelsize=6)

            formatter = ScalarFormatter()
            formatter.set_scientific(False)
            cax.xaxis.set_major_formatter(formatter)
            cbar.ax.tick_params(which='minor', top=False, bottom=False)
            cbar.ax.tick_params(which='minor', labeltop=False, labelbottom=False)    
            cbar.ax.set_xticklabels(['{:.1f}'.format(x) for x in cbar_ticks])

            ax.set_ylabel('Dec. [J2000]', labelpad=-0.1)
            ax.set_xlabel('R.A. [J2000]')
            ax.tick_params(labelsize=6)
            ax.set_facecolor('grey')
            ax.tick_params(axis='both', which='both', colors='white', labelcolor='k')

            # labels
            ax.text(0.97, 0.97, r'channel map at 0 km/s', ha='right', va='top', transform=ax.transAxes, c='white', bbox=box, zorder=10, size=8)
            
            # figure title
            fig.suptitle(galaxy.upper(), size=30, y=1.25)

            ########################################
            # save figure
            plotname = f'{galaxy}_pv_diagrams_and_0kms_channel_map'
            plt.savefig(plot_dir + plotname + '.png', facecolor='white', transparent=False, bbox_inches='tight')
            plt.close()

        print(f'[INFO] {galaxy.upper()}: Shuffled cubes; successfully finished.')

    print('[INFO] Step 1:\tSuccessfully finished.') 


if run_flatmaps:

    ######################################################################################################
    # PIPELINE (2): FLAT MAPS
    ######################################################################################################
    print('[INFO] Step 2:\tInitialise flat maps script.')

    # loop over galaxies
    for galaxy, v_width_narrow in zip(galaxies, v_width_narrow_list):

        print(f'[INFO] {galaxy.upper()}: Start producing flat maps.')

        ########################################
        # LOAD DATA
        ########################################

        # ALMA array
        if galaxy in galaxies_alma_7m:
            array = '7m+tp'
        else:
            array = '12m+7m+tp'
            
        # ALMA line (accounting for post-LP sample)
        if galaxy in galaxies_alma_postlp:
            line = 'co21lores'
        else:
            line = 'co21'

        # CO shuffled cube
        fname_co_cube = shuffled_cubes_dir + f'{galaxy}_{array}_{line}_shuffled.fits'
        co_cube, hdr_co_cube = fits.getdata(fname_co_cube, header=True)
        
        # CO shuffled noise cube
        fname_co_cube_noise = shuffled_cubes_dir + f'{galaxy}_{array}_{line}_shuffled_noise.fits'
        co_cube_noise, hdr_co_cube_noise = fits.getdata(fname_co_cube_noise, header=True)

        # CO shuffled broad mask
        fname_co_broadmask = shuffled_cubes_dir + f'{galaxy}_{array}_{line}_shuffled_broadmask.fits'
        co_broadmask, hdr_co_broadmask = fits.getdata(fname_co_broadmask, header=True)

        # CO shuffled strict mask
        fname_co_strictmask = shuffled_cubes_dir + f'{galaxy}_{array}_{line}_shuffled_strictmask.fits'
        co_strictmask, hdr_co_strictmask = fits.getdata(fname_co_strictmask, header=True)
        
        # CO mom0 (DR4)
        fname_co_mom0 = co_maps_dir + f'{galaxy}_{array}_{line}_broad_mom0.fits'
        co_mom0, hdr_co_mom0 = fits.getdata(fname_co_mom0, header=True)

        # get velocity axis
        n_ch = hdr_co_cube['NAXIS3'] # number of channels
        v_ch = hdr_co_cube['CDELT3'] # channel width (m/s, with sign)
        v_ref = hdr_co_cube['CRVAL3'] # reference channel (m/s)
        ch_ref = hdr_co_cube['CRPIX3'] # reference channel
        v_ch0 = v_ref - (ch_ref-1)*v_ch # velocity of first channel
        vaxis = np.linspace(v_ch0, v_ch0+(n_ch-1)*v_ch, n_ch) * 1e-3 # velocity axis [km/s]

        if diagnostics:

            ######################################################################################################
            # VISULISE VELOCITY WINDOWS
            ######################################################################################################
                
            # load PHANGS table
            phangs_table = pd.read_csv(fpath_sample_table, comment='#', skiprows=1)
            
            # get position angle and centre coordinates of galaxy
            if galaxy == 'circinus':
                phangs_table_galaxy = phangs_table[phangs_table['name'] == 'eso097-013']  # select row in table
            else:
                phangs_table_galaxy = phangs_table[phangs_table['name'] == galaxy]  # select row in table
            ctr_ra = float(phangs_table_galaxy['orient_ra'].iloc[0])  # get galaxy centre right ascension
            ctr_dec = float(phangs_table_galaxy['orient_dec'].iloc[0])  # get galaxy centre declination
            posang = float(phangs_table_galaxy['orient_posang'].iloc[0])  # get position angle
                
            # filenames
            fname = shuffled_cubes_dir + f'{galaxy}_{array}_{line}_shuffled.fits'
            label = 'shuffled cube'
            
            # labels
            box = dict(boxstyle='round', facecolor='k', linewidth=0.5, alpha=0.8, pad=0.2)  # set box

            # create figure
            fig = plt.figure(figsize=(6, 6))
            fig.subplots_adjust(wspace=0.3, hspace=0.5)

            # get pv-diagram data for CO cube
            pv_data, position, velocity = get_pv_data(fits_cube=fname, ctr_ra=ctr_ra, ctr_dec=ctr_dec, posang=posang)  
            
            ########################################
            
            # create subplot
            ax = fig.add_subplot(111)
            
            # plot PV diagram
            vmax = np.nanmax(pv_data)
            im = ax.imshow(pv_data, # plot pv data as 2d array
                        cmap='gist_rainbow_r',  # colormap
                        vmin=0, vmax=1e-2, # color stretch
                        extent=[position[0]*3600, position[-1]*3600, velocity[0], velocity[-1]],  # set axis values according to position, velocity values
                        aspect='auto'  # stretch pixels according to figure size
                        )
            
            # colourbar
            vmin = 0.01
            cbar_ticks = np.logspace(np.log10(vmin), np.log10(vmax), 5)
            cbar_ticks = [0] + list(cbar_ticks[:-1]) + [vmax]
            cax = ax.inset_axes([0, 1, 1, 0.05], transform=ax.transAxes, zorder=5)
            cbar = plt.colorbar(im, cax=cax, orientation='horizontal', ticks=cbar_ticks)
            cax.xaxis.set_ticks_position('top')
            cbar.ax.set_title('$T_\mathrm{CO(2-1)}$ [K]', pad=10, size=12)
            # axis settings
            formatter = ScalarFormatter()
            formatter.set_scientific(False)
            cax.xaxis.set_major_formatter(formatter)
            cbar.ax.tick_params(which='minor', top=False, bottom=False)
            cbar.ax.tick_params(which='minor', labeltop=False, labelbottom=False)    
            cbar.ax.set_xticklabels(['{:.2f}'.format(x) for x in cbar_ticks])

            # axis labels
            ax.set_xlabel('Position along major axis [arcsecond]', size=12)
            ax.set_ylabel('Velocity [km/s]', size=12)
            ax.tick_params(axis='both', which='both', colors='white', labelcolor='k')
            ax.tick_params(labelsize=8)
            ax.set_yticks(np.linspace(-200, 200, 41), minor=True)
            ax.set_ylim(-200,200)
            ax.grid(axis='y', which='both', ls='dotted', c='k')

            # indicate varying velocity window
            ax.axhline(v_width_narrow/2, ls='solid', c='k', lw=2)
            ax.axhline(-v_width_narrow/2, ls='solid', c='k', lw=2)

            # indicate fixed velocity window
            ax.axhline(v_width_50kms/2, ls='dashed', c='k', lw=2)
            ax.axhline(-v_width_50kms/2, ls='dashed', c='k', lw=2)

            # indicate fixed velocity window
            ax.axhline(v_width_100kms/2, ls='dotted', c='k', lw=2)
            ax.axhline(-v_width_100kms/2, ls='dotted', c='k', lw=2)
                
            # labels
            ax.text(0.97, 0.97, label, ha='right', va='top', transform=ax.transAxes, c='white', bbox=box, zorder=10, size=8)

            # title
            fig.suptitle(galaxy.upper(), size=30, y=1.1)
            
            ########################################
            # save figure
            plotname = f'{galaxy}_pv_diagram_velocity_window'
            plt.savefig(plot_dir + plotname + '.png', facecolor='white', transparent=False, bbox_inches='tight') 
            plt.close()


        ######################################################################################################
        # (2.1) FIXED WINDOW FLAT MAPS
        ######################################################################################################

        # loop over velocity windows
        for v_width in [v_width_50kms, v_width_100kms]:
        
            ########################################
            # FLAT MASK
            ########################################
        
            # mask for velocity integration
            flat_mask = np.zeros_like(co_cube)
            mask_vrange = np.where(np.abs(vaxis) <= v_width)[0]
            flat_mask[mask_vrange,:,:] = 1
            
            ########################################
            # FLAT MAP
            ########################################
            
            # flat map (mom0)
            co_mom0_flat = np.nansum(co_cube*flat_mask, axis=0) * np.abs(v_ch) * 1e-3  # K km/s
            co_mom0_flat[np.isnan(co_mom0)] = np.nan
        
            # uncertainty (emom0)
            co_emom0_flat = np.sqrt(np.nansum((co_cube_noise*flat_mask)**2, axis=0)) * np.abs(v_ch) * 1e-3  # K km/s
            co_emom0_flat[np.isnan(co_mom0)] = np.nan
        
            # save maps
            fits.writeto(flat_maps_dir + f'{galaxy}_{array}_{line}_flat_{v_width}kms_mom0.fits', co_mom0_flat, hdr_co_mom0, overwrite=True)
            fits.writeto(flat_maps_dir + f'{galaxy}_{array}_{line}_flat_{v_width}kms_emom0.fits', co_emom0_flat, hdr_co_mom0, overwrite=True)

        print(f'[INFO] {galaxy.upper()}: Step 2.1: Fixed window flat maps; successfully finished.')

        ######################################################################################################
        # (2.2) BROAD FLAT MAPS
        ######################################################################################################

        ########################################
        # NARROW BROAD MASK
        ########################################
        
        # mask for velocity integration
        flat_mask = np.zeros_like(co_cube)  # initiate cube
        mask_vrange = np.where(np.abs(vaxis) <= v_width_narrow/2)[0]  # index within fixed velocity window
        flat_mask[mask_vrange,:,:] = 1  # set velocity window to one
        flat_mask[co_broadmask == 1] = 1  # extend mask to shuffled strict mask to capture signal outside of the fixed window
        flat_mask = flat_mask.astype(np.uint8)
        
        ########################################
        # NARROW BROAD FLAT MAP
        ########################################
        
        # flat map (mom0)
        co_mom0_flat = np.nansum(co_cube*flat_mask, axis=0) * np.abs(v_ch) * 1e-3  # K km/s
        co_mom0_flat[np.isnan(co_mom0)] = np.nan

        # flat noise map (mom0)
        co_emom0_flat = np.sqrt(np.nansum((co_cube_noise*flat_mask)**2, axis=0)) * np.abs(v_ch) * 1e-3  # K km/s
        co_emom0_flat[np.isnan(co_mom0)] = np.nan

        # save maps
        fits.writeto(flat_maps_dir + f'{galaxy}_{array}_{line}_flat_narrow_broad_mom0.fits', co_mom0_flat, hdr_co_mom0, overwrite=True)
        fits.writeto(flat_maps_dir + f'{galaxy}_{array}_{line}_flat_narrow_broad_emom0.fits', co_emom0_flat, hdr_co_mom0, overwrite=True)

        ########################################
        # WIDE BROAD MASK
        ########################################

        # mask for velocity integration
        flat_mask = np.zeros_like(co_cube)  # initiate cube
        mask_vrange = np.where(np.abs(vaxis) <= v_width_wide/2)[0]  # index within fixed velocity window
        flat_mask[mask_vrange,:,:] = 1  # set velocity window to one
        flat_mask[co_broadmask == 1] = 1  # extend mask to shuffled strict mask to capture signal outside of the fixed window
        flat_mask = flat_mask.astype(np.uint8)
        
        # save mask (cube) to fits file
        fits.writeto(shuffled_cubes_dir + f'{galaxy}_{array}_{line}_shuffled_wide_broadmask.fits', flat_mask, hdr_co_broadmask, overwrite=True)

        ########################################
        # WIDE BROAD FLAT MAP
        ########################################
        
        # flat map (mom0)
        co_mom0_flat = np.nansum(co_cube*flat_mask, axis=0) * np.abs(v_ch) * 1e-3  # K km/s
        co_mom0_flat[np.isnan(co_mom0)] = np.nan

        # flat noise map (mom0)
        co_emom0_flat = np.sqrt(np.nansum((co_cube_noise*flat_mask)**2, axis=0)) * np.abs(v_ch) * 1e-3  # K km/s
        co_emom0_flat[np.isnan(co_mom0)] = np.nan

        # save maps
        fits.writeto(flat_maps_dir + f'{galaxy}_{array}_{line}_flat_wide_broad_mom0.fits', co_mom0_flat, hdr_co_mom0, overwrite=True)
        fits.writeto(flat_maps_dir + f'{galaxy}_{array}_{line}_flat_wide_broad_emom0.fits', co_emom0_flat, hdr_co_mom0, overwrite=True)

        print(f'[INFO] {galaxy.upper()}: Step 2.2: Broad flat maps; successfully finished.')

        ######################################################################################################
        # (2.3) STRICT FLAT MAPS
        ######################################################################################################

        ########################################
        # NARROW STRICT MASK
        ########################################
        
        # mask for velocity integration
        flat_mask = np.zeros_like(co_cube)  # initiate cube
        mask_vrange = np.where(np.abs(vaxis) <= v_width_narrow/2)[0]  # index within fixed velocity window
        flat_mask[mask_vrange,:,:] = 1  # set velocity window to one
        flat_mask[co_strictmask == 1] = 1  # extend mask to shuffled strict mask to capture signal outside of the fixed window
        flat_mask = flat_mask.astype(np.uint8)
        
        # save mask (cube) to fits file
        fits.writeto(shuffled_cubes_dir + f'{galaxy}_{array}_{line}_shuffled_narrow_strictmask.fits', flat_mask, hdr_co_strictmask, overwrite=True)

        ########################################
        # NARROW STRICT FLAT MAP
        ########################################
        
        # flat map (mom0)
        co_mom0_flat = np.nansum(co_cube*flat_mask, axis=0) * np.abs(v_ch) * 1e-3  # K km/s
        co_mom0_flat[np.isnan(co_mom0)] = np.nan

        # flat noise map (mom0)
        co_emom0_flat = np.sqrt(np.nansum((co_cube_noise*flat_mask)**2, axis=0)) * np.abs(v_ch) * 1e-3  # K km/s
        co_emom0_flat[np.isnan(co_mom0)] = np.nan

        # save maps
        fits.writeto(flat_maps_dir + f'{galaxy}_{array}_{line}_flat_narrow_strict_mom0.fits', co_mom0_flat, hdr_co_mom0, overwrite=True)
        fits.writeto(flat_maps_dir + f'{galaxy}_{array}_{line}_flat_narrow_strict_emom0.fits', co_emom0_flat, hdr_co_mom0, overwrite=True)

        ########################################
        # WIDE STRICT MASK
        ########################################

        # mask for velocity integration
        flat_mask = np.zeros_like(co_cube)  # initiate cube
        mask_vrange = np.where(np.abs(vaxis) <= v_width_wide/2)[0]  # index within fixed velocity window
        flat_mask[mask_vrange,:,:] = 1  # set velocity window to one
        flat_mask[co_strictmask == 1] = 1  # extend mask to shuffled strict mask to capture signal outside of the fixed window
        flat_mask = flat_mask.astype(np.uint8)
        
        ########################################
        # WIDE STRICT FLAT MAP
        ########################################
        
        # flat map (mom0)
        co_mom0_flat = np.nansum(co_cube*flat_mask, axis=0) * np.abs(v_ch) * 1e-3  # K km/s
        co_mom0_flat[np.isnan(co_mom0)] = np.nan

        # flat noise map (mom0)
        co_emom0_flat = np.sqrt(np.nansum((co_cube_noise*flat_mask)**2, axis=0)) * np.abs(v_ch) * 1e-3  # K km/s
        co_emom0_flat[np.isnan(co_mom0)] = np.nan

        # save maps
        fits.writeto(flat_maps_dir + f'{galaxy}_{array}_{line}_flat_wide_strict_mom0.fits', co_mom0_flat, hdr_co_mom0, overwrite=True)
        fits.writeto(flat_maps_dir + f'{galaxy}_{array}_{line}_flat_wide_strict_emom0.fits', co_emom0_flat, hdr_co_mom0, overwrite=True)

        print(f'[INFO] {galaxy.upper()}: Step 2.3: Strict flat maps; successfully finished.')


        ######################################################################################################
        # (2.4) RE-SHUFFLED MASKS
        ######################################################################################################]

        # load velocity field
        fname_vfield = vel_maps_dir + f'{galaxy}_combined_velocity_field.fits'
        vfield, hdr_vfield = fits.getdata(fname_vfield, header=True)

        # produce re-shuffled masks for two mask versions
        mask_list = ['shuffled_wide_broadmask', 'shuffled_narrow_strictmask']
        mask_list_new = ['wide_broadmask', 'narrow_strictmask']

        # loop over masks
        for mask, mask_new in zip(mask_list, mask_list_new):

            # load shuffled mask
            fname_mask = shuffled_cubes_dir + f'{galaxy}_{array}_{line}_{mask}.fits'
            data_mask, hdr_mask = fits.getdata(fname_mask, header=True)

            # load (unshuffled) CO cube
            fname_cube = co_cubes_dir + f'{galaxy}_{array}_{line}.fits'
            data_cube, hdr_cube = fits.getdata(fname_cube, header=True)
            
            ########################################
            # RE-SHUFFLING OF VELOCITY-INTEGRATION MASK
            ########################################
            
            # get velocity axis
            n_ch = hdr_cube['NAXIS3'] # number of channels
            v_ch = hdr_cube['CDELT3'] # channel width (m/s, with sign)
            v_ref = hdr_cube['CRVAL3'] # reference channel (m/s)
            ch_ref = hdr_cube['CRPIX3'] # reference channel
            v_ch0 = v_ref - (ch_ref-1)*v_ch # velocity of first channel
            vaxis = np.linspace(v_ch0, v_ch0+(n_ch-1)*v_ch, n_ch) # velocity axis
            if hdr_mask['CUNIT3'] == 'm/s':
                vaxis *= 1e-3  # convert to km/s
            else:
                print('[WARNING] Velocity axis units not specified, assuming km/s')

            # edit header to match original velocity axis
            hdr_mask['CRPIX3'] = hdr_cube['CRPIX3']
            hdr_mask['CRVAL3'] = hdr_cube['CRVAL3']
            hdr_mask['CDELT3'] = hdr_cube['CDELT3']
            
            # shuffle mask
            mask_reshuffled = shuffle_cube(data_mask, np.flip(vaxis), vfield)

            # save reshuffled mask
            fname_mask_reshuffled = shuffled_cubes_dir + f'{galaxy}_{array}_{line}_{mask_new}.fits'
            fits.writeto(fname_mask_reshuffled, mask_reshuffled, hdr_mask, overwrite=True)

        print(f'[INFO] {galaxy.upper()}: Step 2.4: Reshuffled masks; successfully finished.')

        ######################################################################################################

        if diagnostics:

            ######################################################################################################
            # PLOT PV DIAGRAMS AND FLAT MAPS
            ######################################################################################################

            # get position angle and centre coordinates of galaxy
            if galaxy == 'circinus':
                phangs_table_galaxy = phangs_table[phangs_table['name'] == 'eso097-013']  # select row in table
            else:
                phangs_table_galaxy = phangs_table[phangs_table['name'] == galaxy]  # select row in table
            ctr_ra = float(phangs_table_galaxy['orient_ra'].iloc[0])  # get galaxy centre right ascension
            ctr_dec = float(phangs_table_galaxy['orient_dec'].iloc[0])  # get galaxy centre declination
            posang = float(phangs_table_galaxy['orient_posang'].iloc[0])  # get position angle
            
            # filenames
            co_cube_fname = co_cubes_dir + f'{galaxy}_{array}_{line}.fits'
            co_cube_shuffled_fname = shuffled_cubes_dir + f'{galaxy}_{array}_{line}_shuffled.fits'
            co_cube_shuffled_noise_fname = shuffled_cubes_dir + f'{galaxy}_{array}_{line}_shuffled_noise.fits'
            co_broadmask_shuffled_fname = shuffled_cubes_dir + f'{galaxy}_{array}_{line}_shuffled_wide_broadmask.fits'
            co_strictmask_shuffled_fname = shuffled_cubes_dir+ f'{galaxy}_{array}_{line}_shuffled_narrow_strictmask.fits'
            
            # create lists to be iterated
            fname_list = [co_cube_fname, co_cube_shuffled_fname, co_cube_shuffled_fname, co_cube_shuffled_fname]
            labels = ['original cube', 'shuffled cube', 'shuffled cube', 'shuffled cube']
            labels_mask = ['none', '50 km/s mask', 'wide+broad mask', 'narrow+strict mask']
            
            # labels
            box = dict(boxstyle='round', facecolor='k', linewidth=0.5, alpha=0.8, pad=0.2)  # set box

            # gridspec inside gridspec
            fig = plt.figure(figsize=(16, 7))
            gs0 = gridspec.GridSpec(2, 1, figure=fig, hspace=0.4)
            # upper grid
            gs00 = gs0[0].subgridspec(1, 4, wspace=0.2)
            # lower grid
            gs01 = gs0[1].subgridspec(1, 7, wspace=0.3)

            n = 0
            
            # loop over cubes and get PV diagram data
            for fname, label, label_mask in zip(fname_list, labels, labels_mask):
            
                # get pv-diagram data for CO cube
                pv_data, position, velocity = get_pv_data(fits_cube = fname, ctr_ra=ctr_ra, ctr_dec=ctr_dec, posang=posang)  
                
                # create subplot
                ax = fig.add_subplot(gs00[n])
                
                # plot PV diagram
                vmax = np.nanmax(pv_data)
                im = ax.imshow(pv_data, # plot pv data as 2d array
                            cmap='gist_heat',  # colormap
                            norm=SymLogNorm(linthresh=np.nanpercentile(pv_data, 90), vmin=0, vmax=vmax), # color stretch
                            extent=[position[0]*3600, position[-1]*3600, velocity[0], velocity[-1]],  # set axis values according to position, velocity values
                            aspect='auto'  # stretch pixels according to figure size
                            )
                
                # colorbar
                vmin = 0.01
                cbar_ticks = np.logspace(np.log10(vmin), np.log10(vmax), 5)
                cbar_ticks = [0] + list(cbar_ticks[:-1]) + [vmax]
                cax = ax.inset_axes([0, 1, 1, 0.05], transform=ax.transAxes, zorder=5)
                cbar = plt.colorbar(im, cax=cax, orientation='horizontal', ticks=cbar_ticks)
                cax.xaxis.set_ticks_position('top')
                cbar.ax.set_title('$T_\mathrm{CO(2-1)}$ [K]', pad=10, size=12)
            
                formatter = ScalarFormatter()
                formatter.set_scientific(False)
                cax.xaxis.set_major_formatter(formatter)
                cbar.ax.tick_params(which='minor', top=False, bottom=False)
                cbar.ax.tick_params(which='minor', labeltop=False, labelbottom=False)    
                cbar.ax.set_xticklabels(['{:.2f}'.format(x) for x in cbar_ticks])
                cbar.ax.tick_params(labelsize=8)  

                # axis labels
                ax.set_xlabel('Position along major axis [arcsecond]', size=12)
                if n == 0:
                    ax.set_ylabel('Velocity [km/s]', size=12)
                ax.tick_params(axis='both', which='both', colors='white', labelcolor='k')
                ax.tick_params(labelsize=8)  
                
                if n == 1:
                    # show mask contours
                    ax.axhline(-50/2, ls='solid', c='white')
                    ax.axhline(50/2, ls='solid', c='white')
                elif n == 2:
                    # show mask contours
                    pv_data_mask, position_mask, velocity_mask = get_pv_data(fits_cube=co_broadmask_shuffled_fname, ctr_ra=ctr_ra, ctr_dec=ctr_dec, posang=posang)  
                    ax.contour(pv_data_mask, levels=[0.01], colors='white', extent=[position_mask[0]*3600, position_mask[-1]*3600, velocity_mask[0], velocity_mask[-1]])
                elif n == 3:
                    # show mask contours
                    pv_data_mask, position_mask, velocity_mask = get_pv_data(fits_cube=co_strictmask_shuffled_fname, ctr_ra=ctr_ra, ctr_dec=ctr_dec, posang=posang)  
                    ax.contour(pv_data_mask, levels=[0.01], colors='white', extent=[position_mask[0]*3600, position_mask[-1]*3600, velocity_mask[0], velocity_mask[-1]])

                ax.set_xlim(position[0]*3600, position[-1]*3600)
                    
                # labels
                ax.text(0.97, 0.97, label, ha='right', va='top', transform=ax.transAxes, c='white', bbox=box, zorder=10, size=8)
                if n > 0:
                    ax.text(0.97, 0.03, label_mask, ha='right', va='bottom', transform=ax.transAxes, c='white', bbox=box, zorder=10, size=8)

                n += 1

            
            ########################################
            # PHANGS-ALMA broad mom0 map (DR4)
            ########################################
            
            # load mom0 map
            fname_co_mom0 = co_maps_dir + f'{galaxy}_{array}_{line}_broad_mom0.fits'
            co_mom0, hdr_co_mom0 = fits.getdata(fname_co_mom0, header=True)
            wcs = WCS(hdr_co_mom0)
            
            # combined velocity field
            ax = fig.add_subplot(gs01[0], projection=wcs)
            linthresh = max(np.nanpercentile(co_mom0, 90), vmax/100)
            vmax = np.nanmax(co_mom0)
            im = ax.imshow(co_mom0, 
                        transform=ax.get_transform(wcs), cmap='gist_heat', 
                        norm=SymLogNorm(linthresh=linthresh, vmin=0, vmax=vmax), # color stretch
                        )
            vmin = 3
            cbar_ticks = np.logspace(np.log10(vmin), np.log10(vmax), 4)
            cbar_ticks = [0] + list(cbar_ticks[:-1]) + [vmax]
            cax = ax.inset_axes([0, 1, 1, 0.05], transform=ax.transAxes, zorder=5)
            cbar = plt.colorbar(im, cax=cax, orientation='horizontal', ticks=cbar_ticks)
            cax.xaxis.set_ticks_position('top')
            cbar.ax.set_title(r'$W_\mathrm{CO(2-1)}$ [$\rm{K}\,\rm{km}\,\rm{s}^{-1}$]', pad=10, size=12)

            formatter = ScalarFormatter()
            formatter.set_scientific(False)
            cax.xaxis.set_major_formatter(formatter)
            cbar.ax.tick_params(which='minor', top=False, bottom=False)
            cbar.ax.tick_params(which='minor', labeltop=False, labelbottom=False)    
            cbar.ax.set_xticklabels(['{:.0f}'.format(x) for x in cbar_ticks])
            cbar.ax.tick_params(labelsize=8)   

            ax.set_ylabel('Dec. [J2000]', labelpad=-0.1, size=12)
            ax.set_xlabel('R.A. [J2000]', size=12)
            ax.set_facecolor('grey')
            ax.tick_params(axis='both', which='both', colors='white', labelcolor='k', labelsize=6)

            # labels
            ax.text(0.97, 0.97, r'broad mom0 (DR4)', ha='right', va='top', transform=ax.transAxes, c='white', bbox=box, zorder=10, size=6)


            ########################################
            # FLAT NAIVE MAP (50 km/s)
            ########################################

            # load flat map
            fname_co_mom0_flat_naive = flat_maps_dir + f'{galaxy}_{array}_{line}_flat_50kms_mom0.fits'
            co_mom0_flat_naive = fits.getdata(fname_co_mom0_flat_naive, header=False)
            
            # combined velocity field
            ax = fig.add_subplot(gs01[1], projection=wcs)
            vmax = np.nanmax(co_mom0_flat_naive)
            linthresh = max(np.nanpercentile(co_mom0_flat_naive, 90), vmax/100)
            im = ax.imshow(co_mom0_flat_naive, 
                        transform=ax.get_transform(wcs), cmap='gist_heat', 
                        norm=SymLogNorm(linthresh=linthresh, vmin=0, vmax=vmax), # color stretch
                        )
            vmin = 3
            cbar_ticks = np.logspace(np.log10(vmin), np.log10(vmax), 4)
            cbar_ticks = [0] + list(cbar_ticks[:-1]) + [vmax]
            cax = ax.inset_axes([0, 1, 1, 0.05], transform=ax.transAxes, zorder=5)
            cbar = plt.colorbar(im, cax=cax, orientation='horizontal', ticks=cbar_ticks)
            cax.xaxis.set_ticks_position('top')
            cbar.ax.set_title(r'$W_\mathrm{CO(2-1)}$ [$\rm{K}\,\rm{km}\,\rm{s}^{-1}$]', pad=10, size=12)
            cbar.ax.tick_params(labelsize=8)  

            formatter = ScalarFormatter()
            formatter.set_scientific(False)
            cax.xaxis.set_major_formatter(formatter)
            cbar.ax.tick_params(which='minor', top=False, bottom=False)
            cbar.ax.tick_params(which='minor', labeltop=False, labelbottom=False)    
            cbar.ax.set_xticklabels(['{:.0f}'.format(x) for x in cbar_ticks])

            # ax.set_ylabel('Dec. [J2000]', labelpad=-0.1, size=12)
            ax.set_ylabel(' ')
            ax.set_xlabel('R.A. [J2000]', size=12)
            ax.set_facecolor('grey')
            ax.tick_params(axis='both', which='both', colors='white', labelcolor='k', labelsize=6)

            # labels
            ax.text(0.97, 0.97, r'flat 50 km/s mom0 (DR5)', ha='right', va='top', transform=ax.transAxes, c='white', bbox=box, zorder=10, size=6)

            
            ########################################
            # FLAT NAIVE NOISE MAP (50 km/s)
            ########################################
            
            # load flat moise map
            fname_co_emom0_flat_naive = flat_maps_dir + f'{galaxy}_{array}_{line}_flat_50kms_emom0.fits'
            co_emom0_flat_naive = fits.getdata(fname_co_emom0_flat_naive, header=False)
                
            # combined velocity field
            ax = fig.add_subplot(gs01[2], projection=wcs)
            im = ax.imshow(co_emom0_flat_naive, 
                        transform=ax.get_transform(wcs), cmap='gist_rainbow_r', 
                        vmin=np.nanpercentile(co_emom0_flat_naive, 10), vmax=np.nanpercentile(co_emom0_flat_naive, 80),
                        )

            cax = ax.inset_axes([0, 1, 1, 0.05], transform=ax.transAxes, zorder=5)
            cbar = plt.colorbar(im, cax=cax, orientation='horizontal')
            cax.xaxis.set_ticks_position('top')
            cbar.ax.set_title(r'$\sigma_{W_\mathrm{CO(2-1)}}$ [$\rm{K}\,\rm{km}\,\rm{s}^{-1}$]', pad=10, size=12)
            cbar.ax.tick_params(labelsize=8)  

            # ax.set_ylabel('Dec. [J2000]', labelpad=-0.1, size=12)
            ax.set_ylabel(' ')
            ax.set_xlabel('R.A. [J2000]', size=12)
            ax.set_facecolor('grey')
            ax.tick_params(axis='both', which='both', colors='white', labelcolor='k', labelsize=6)

            # labels
            ax.text(0.97, 0.97, r'flat 50 km/s emom0 (DR5)', ha='right', va='top', transform=ax.transAxes, c='white', bbox=box, zorder=10, size=6)

            
            ########################################
            # FLAT BROAD MAP
            ########################################

            # load flat map
            fname_co_mom0_flat_broad = flat_maps_dir + f'{galaxy}_{array}_{line}_flat_wide_broad_mom0.fits'
            co_mom0_flat_broad = fits.getdata(fname_co_mom0_flat_broad, header=False)

            # combined velocity field
            ax = fig.add_subplot(gs01[3], projection=wcs)
            vmax = np.nanmax(co_mom0_flat_broad)
            linthresh = max(np.nanpercentile(co_mom0_flat_broad, 90), vmax/100)
            im = ax.imshow(co_mom0_flat_broad, 
                        transform=ax.get_transform(wcs), cmap='gist_heat', 
                        norm=SymLogNorm(linthresh=linthresh, vmin=0, vmax=vmax), # color stretch
                        )
            vmin = 3
            cbar_ticks = np.logspace(np.log10(vmin), np.log10(vmax), 4)
            cbar_ticks = [0] + list(cbar_ticks[:-1]) + [vmax]
            cax = ax.inset_axes([0, 1, 1, 0.05], transform=ax.transAxes, zorder=5)
            cbar = plt.colorbar(im, cax=cax, orientation='horizontal', ticks=cbar_ticks)
            cax.xaxis.set_ticks_position('top')
            cbar.ax.set_title(r'$W_\mathrm{CO(2-1)}$ [$\rm{K}\,\rm{km}\,\rm{s}^{-1}$]', pad=10, size=12)
            cbar.ax.tick_params(labelsize=8)  

            formatter = ScalarFormatter()
            formatter.set_scientific(False)
            cax.xaxis.set_major_formatter(formatter)
            cbar.ax.tick_params(which='minor', top=False, bottom=False)
            cbar.ax.tick_params(which='minor', labeltop=False, labelbottom=False)    
            cbar.ax.set_xticklabels(['{:.0f}'.format(x) for x in cbar_ticks])

            # ax.set_ylabel('Dec. [J2000]', labelpad=-0.1, size=12)
            ax.set_ylabel(' ')
            ax.set_xlabel('R.A. [J2000]', size=12)
            ax.set_facecolor('grey')
            ax.tick_params(axis='both', which='both', colors='white', labelcolor='k', labelsize=6)

            # labels
            ax.text(0.97, 0.97, r'wide+broad mom0 (DR5)', ha='right', va='top', transform=ax.transAxes, c='white', bbox=box, zorder=10, size=6)

            
            ########################################
            # FLAT BROAD NOISE MAP
            ########################################
            
            # load flat moise map
            fname_co_emom0_flat_broad = flat_maps_dir + f'{galaxy}_{array}_{line}_flat_wide_broad_emom0.fits'
            co_emom0_flat_broad = fits.getdata(fname_co_emom0_flat_broad, header=False)
                
            # combined velocity field
            ax = fig.add_subplot(gs01[4], projection=wcs)
            im = ax.imshow(co_emom0_flat_broad, 
                        transform=ax.get_transform(wcs), cmap='gist_rainbow_r', 
                        vmin=np.nanpercentile(co_emom0_flat_broad, 10), vmax=np.nanpercentile(co_emom0_flat_broad, 80),
                        )

            cax = ax.inset_axes([0, 1, 1, 0.05], transform=ax.transAxes, zorder=5)
            cbar = plt.colorbar(im, cax=cax, orientation='horizontal')
            cax.xaxis.set_ticks_position('top')
            cbar.ax.set_title(r'$\sigma_{W_\mathrm{CO(2-1)}}$ [$\rm{K}\,\rm{km}\,\rm{s}^{-1}$]', pad=10, size=12)
            cbar.ax.tick_params(labelsize=8)  

            # ax.set_ylabel('Dec. [J2000]', labelpad=-0.1, size=12)
            ax.set_ylabel(' ')
            ax.set_xlabel('R.A. [J2000]', size=12)
            ax.set_facecolor('grey')
            ax.tick_params(axis='both', which='both', colors='white', labelcolor='k', labelsize=6)

            # labels
            ax.text(0.97, 0.97, r'wide+broad emom0 (DR5)', ha='right', va='top', transform=ax.transAxes, c='white', bbox=box, zorder=10, size=6)


            ########################################
            # FLAT STRICT MAP
            ########################################

            # load flat map
            fname_co_mom0_flat_strict = flat_maps_dir + f'{galaxy}_{array}_{line}_flat_narrow_strict_mom0.fits'
            co_mom0_flat_strict = fits.getdata(fname_co_mom0_flat_strict, header=False)

            # combined velocity field
            ax = fig.add_subplot(gs01[5], projection=wcs)
            vmax = np.nanmax(co_mom0_flat_strict)
            linthresh = max(np.nanpercentile(co_mom0_flat_strict, 90), vmax/100)
            im = ax.imshow(co_mom0_flat_strict, 
                        transform=ax.get_transform(wcs), cmap='gist_heat', 
                        norm=SymLogNorm(linthresh=linthresh, vmin=0, vmax=vmax), # color stretch
                        )
            vmin = 3
            cbar_ticks = np.logspace(np.log10(vmin), np.log10(vmax), 4)
            cbar_ticks = [0] + list(cbar_ticks[:-1]) + [vmax]
            cax = ax.inset_axes([0, 1, 1, 0.05], transform=ax.transAxes, zorder=5)
            cbar = plt.colorbar(im, cax=cax, orientation='horizontal', ticks=cbar_ticks)
            cax.xaxis.set_ticks_position('top')
            cbar.ax.set_title(r'$W_\mathrm{CO(2-1)}$ [$\rm{K}\,\rm{km}\,\rm{s}^{-1}$]', pad=10, size=12)
            cbar.ax.tick_params(labelsize=8)  

            formatter = ScalarFormatter()
            formatter.set_scientific(False)
            cax.xaxis.set_major_formatter(formatter)
            cbar.ax.tick_params(which='minor', top=False, bottom=False)
            cbar.ax.tick_params(which='minor', labeltop=False, labelbottom=False)    
            cbar.ax.set_xticklabels(['{:.0f}'.format(x) for x in cbar_ticks])

            # ax.set_ylabel('Dec. [J2000]', labelpad=-0.1, size=12)
            ax.set_ylabel(' ')
            ax.set_xlabel('R.A. [J2000]', size=12)
            ax.set_facecolor('grey')
            ax.tick_params(axis='both', which='both', colors='white', labelcolor='k', labelsize=6)

            # labels
            ax.text(0.97, 0.97, r'narrow+strict mom0 (DR5)', ha='right', va='top', transform=ax.transAxes, c='white', bbox=box, zorder=10, size=6)

            
            ########################################
            # FLAT BROAD NOISE MAP
            ########################################
            
            # load flat moise map
            fname_co_emom0_flat_strict = flat_maps_dir + f'{galaxy}_{array}_{line}_flat_narrow_strict_emom0.fits'
            co_emom0_flat_strict = fits.getdata(fname_co_emom0_flat_strict, header=False)
                
            # combined velocity field
            ax = fig.add_subplot(gs01[6], projection=wcs)
            im = ax.imshow(co_emom0_flat_strict, 
                        transform=ax.get_transform(wcs), cmap='gist_rainbow_r', 
                        vmin=np.nanpercentile(co_emom0_flat_strict, 10), vmax=np.nanpercentile(co_emom0_flat_strict, 80),
                        )

            cax = ax.inset_axes([0, 1, 1, 0.05], transform=ax.transAxes, zorder=5)
            cbar = plt.colorbar(im, cax=cax, orientation='horizontal')
            cax.xaxis.set_ticks_position('top')
            cbar.ax.set_title(r'$\sigma_{W_\mathrm{CO(2-1)}}$ [$\rm{K}\,\rm{km}\,\rm{s}^{-1}$]', pad=10, size=12)
            cbar.ax.tick_params(labelsize=8)  

            # axis settings
            ax.set_ylabel(' ')
            ax.set_xlabel('R.A. [J2000]', size=12)
            ax.set_facecolor('grey')
            ax.tick_params(axis='both', which='both', colors='white', labelcolor='k', labelsize=6)

            # labels
            ax.text(0.97, 0.97, r'narrow+strict emom0 (DR5)', ha='right', va='top', transform=ax.transAxes, c='white', bbox=box, zorder=10, size=6)

            # title
            fig.suptitle(galaxy.upper(), size=30, y=1.05)
            
            ########################################
            # save figure
            plotname = f'{galaxy}_pv_diagrams_and_flat_maps'
            plt.savefig(plot_dir + plotname + '.png', facecolor='white', transparent=False, bbox_inches='tight') 
            plt.close()

    print('[INFO] Step 2:\tSuccessfully finished.')
    print('[INFO] Pipeline successfully finished.')