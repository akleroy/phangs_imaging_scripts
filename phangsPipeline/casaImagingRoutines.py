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
import casaStuff as casa

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
        print "File not found."
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
    
    os.system('cp -r '+input_root+'.image '+output_root+'.image')
    os.system('cp -r '+input_root+'.model '+output_root+'.model')
    os.system('cp -r '+input_root+'.mask '+output_root+'.mask')
    os.system('cp -r '+input_root+'.pb '+output_root+'.pb')
    os.system('cp -r '+input_root+'.psf '+output_root+'.psf')
    os.system('cp -r '+input_root+'.residual '+output_root+'.residual')
    os.system('cp -r '+input_root+'.psf '+output_root+'.weight')
    os.system('cp -r '+input_root+'.residual '+output_root+'.sumwt')

def replace_imaging_with_copy(
    to_root=None,
    from_root=None):
    """
    Replace a cube with a copy.
    """

    wipe_cube(to_root)

    os.system('cp -r '+from_root+'.image '+to_root+'.image')
    os.system('cp -r '+from_root+'.model '+to_root+'.model')
    os.system('cp -r '+from_root+'.mask '+to_root+'.mask')
    os.system('cp -r '+from_root+'.pb '+to_root+'.pb')
    os.system('cp -r '+from_root+'.psf '+to_root+'.psf')
    os.system('cp -r '+from_root+'.residual '+to_root+'.residual')
    os.system('cp -r '+from_root+'.psf '+to_root+'.weight')
    os.system('cp -r '+from_root+'.residual '+to_root+'.sumwt')

def export_imaging_to_fits(
    image_root=None,
    bitpix=-32):
    """
    Export the various products associated with a CASA cube to FITS.
    """

    exportfits(imagename=image_root+'.image',
               fitsimage=image_root+'.fits',
               velocity=True, overwrite=True, dropstokes=True, 
               dropdeg=True, bitpix=bitpix)

    exportfits(imagename=image_root+'.model',
               fitsimage=image_root+'_model.fits',
               velocity=True, overwrite=True, dropstokes=True, 
               dropdeg=True, bitpix=bitpix)

    exportfits(imagename=image_root+'.residual',
               fitsimage=image_root+'_residual.fits',
               velocity=True, overwrite=True, dropstokes=True, 
               dropdeg=True, bitpix=bitpix)

    exportfits(imagename=image_root+'.mask',
               fitsimage=image_root+'_mask.fits',
               velocity=True, overwrite=True, dropstokes=True, 
               dropdeg=True, bitpix=bitpix)
    
    exportfits(imagename=image_root+'.pb',
               fitsimage=image_root+'_pb.fits',
               velocity=True, overwrite=True, dropstokes=True, 
               dropdeg=True, bitpix=bitpix)

    return


#endregion
