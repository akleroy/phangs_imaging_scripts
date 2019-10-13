# NB: CASA doesn't always include the pwd in the module search path. I
# had to modify my init.py file to get this to import. See the README.

import os
import numpy as np
import pyfits
import scipy.ndimage as ndimage
import glob

# Other PHANGS scripts
from phangsPipelinePython import *
import line_list

# Analysis utilities
import analysisUtils as au

# CASA imports
from taskinit import *

# Import specific CASA tasks
from concat import concat
from exportfits import exportfits
from feather import feather
from flagdata import flagdata
from imhead import imhead
from immath import immath
from impbcor import impbcor
from importfits import importfits
from imrebin import imrebin
from imregrid import imregrid
from imsmooth import imsmooth
from imstat import imstat
from imsubimage import imsubimage
from makemask import makemask
from mstransform import mstransform
from split import split
from statwt import statwt
from tclean import tclean
from uvcontsub import uvcontsub
from visstat import visstat

# Physical constants
sol_kms = 2.9979246e5

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# DIRECTORY STRUCTURE AND FILE MANAGEMENT
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

def rebuild_directories(outroot_dir=None):
    """
    Reset and rebuild the directory structure.
    """
    
    if outroot_dir is None:
        print("Specify a root directory.")
        return

    check_string = raw_input("Reseting directory structure in "+\
                                 outroot_dir+". Type 'Yes' to confirm.")
    if check_string != "Yes":
        print("Aborting")
        return

    if os.is_dir(outroot_dir) == False:
        print("Got make "+outroot_dir+" manually. Then I will make the subdirectories.")

    for subdir in ['raw/','process/','products/','feather/','delivery/']:
        os.system('rm -rf '+outroot_dir+'/'+subdir)       
        os.system('mkdir '+outroot_dir+'/'+subdir)

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# SET UP THE CUBES
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

def phangs_stage_cubes(
    gal=None, array=None, product=None,
    root_dir=None, 
    overwrite=False,
    ):
    """
    Stage cubes for further processing in CASA. This imports them back
    into CASA as .image files.
    """

    if gal is None or array is None or product is None or \
            root_dir is None:
        print("Missing required input.")
        return
    
    input_dir = dir_for_gal(gal)
    in_cube_name = input_dir+gal+'_'+array+'_'+product+'.fits'
    in_pb_name = input_dir+gal+'_'+array+'_'+product+'_pb.fits'

    out_dir = root_dir+'raw/'
    out_cube_name = out_dir+gal+'_'+array+'_'+product+'.image'
    out_pb_name = out_dir+gal+'_'+array+'_'+product+'.pb'
    
    print("... importing data for "+in_cube_name)

    if os.path.isfile(in_cube_name):
        importfits(fitsimage=in_cube_name, imagename=out_cube_name,
                   zeroblanks=True, overwrite=overwrite)
    else:
        print("File not found "+in_cube_name)

    if os.path.isfile(in_pb_name):
        importfits(fitsimage=in_pb_name, imagename=out_pb_name,
                   zeroblanks=True, overwrite=overwrite)
    else:
        print("Directory not found "+in_pb_name)

def phangs_stage_singledish(
    gal=None, product=None, root_dir=None, 
    overwrite=False):
    """
    Copy the single dish data for further processing
    """

    if gal is None or product is None or \
            root_dir is None:
        print("Missing required input.")
        return
    
    sdk = read_singledish_key()
    if (gal in sdk.keys()) == False:
        print(gal+" not found in single dish key.")
        return
    
    this_key = sdk[gal]
    if (product in this_key.keys()) == False:
        print(product+" not found in single dish key for "+gal)
        return
    
    sdfile_in = this_key[product]
    
    sdfile_out = root_dir+'raw/'+gal+'_tp_'+product+'.image'    

    print("... importing single dish data for "+sdfile_in)

    importfits(fitsimage=sdfile_in, imagename=sdfile_out+'.temp',
               zeroblanks=True, overwrite=overwrite)

    if overwrite:
        os.system('rm -rf '+sdfile_out)
    imsubimage(imagename=sdfile_out+'.temp', outfile=sdfile_out,
               dropdeg=True)

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# BASIC IMAGE PROCESSING STEPS
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

def phangs_primary_beam_correct(
    gal=None, array=None, product=None, root_dir=None, 
    cutoff=0.25, overwrite=False):
    """
    Construct primary-beam corrected images using PHANGS naming conventions.
    """

    input_dir = root_dir+'raw/'
    input_cube_name = input_dir+gal+'_'+array+'_'+product+'.image'
    input_pb_name = input_dir+gal+'_'+array+'_'+product+'.pb'
    output_dir = root_dir+'process/'
    output_cube_name = output_dir+gal+'_'+array+'_'+product+'_pbcorr.image'

    print("")
    print("... producing a primary beam corrected image for "+input_cube_name)
    print("")
    
    primary_beam_correct(
        infile=input_cube_name, 
        pbfile=input_pb_name,
        outfile=output_cube_name,
        cutoff=cutoff, overwrite=overwrite)

def primary_beam_correct(
    infile=None, pbfile=None, outfile=None, 
    cutoff=0.25, overwrite=False):
    """
    Construct a primary-beam corrected image.
    """

    if infile is None or pbfile is None or outfile is None:
        print("Missing required input.")
        return

    if os.path.isdir(infile) == False:
        print("Input file missing - "+infile)
        return

    if os.path.isdir(pbfile) == False:
        print("Primary beam file missing - "+pbfile)
        return

    if overwrite:
        os.system('rm -rf '+outfile)

    impbcor(imagename=infile, pbimage=pbfile, outfile=outfile, cutoff=cutoff)

def phangs_convolve_to_round_beam(
    gal=None, array=None, product=None, root_dir=None, 
    force_beam=None, overwrite=False):
    """
    Construct primary-beam corrected images using PHANGS naming
    conventions. Runs on both primary beam corrected cube and flat
    cube, forcing the same beam for both.
    """

    if gal is None or array is None or product is None or \
            root_dir is None:
        print("Missing required input.")
        return
    
    input_dir = root_dir+'process/'
    input_cube_name = input_dir+gal+'_'+array+'_'+product+'_pbcorr.image'
    output_dir = root_dir+'process/'
    output_cube_name = output_dir+gal+'_'+array+'_'+product+'_pbcorr_round.image'

    print("")
    print("... convolving to a round beam for "+input_cube_name)
    print("")

    round_beam = \
        convolve_to_round_beam(
        infile=input_cube_name,
        outfile=output_cube_name,
        force_beam=force_beam,
        overwrite=overwrite)

    print("")
    print("... found beam of "+str(round_beam)+" arcsec. Forcing flat cube to this.")
    print("")

    input_dir = root_dir+'raw/'
    input_cube_name = input_dir+gal+'_'+array+'_'+product+'.image'
    output_dir = root_dir+'process/'
    output_cube_name = output_dir+gal+'_'+array+'_'+product+'_flat_round.image'    

    convolve_to_round_beam(
        infile=input_cube_name,
        outfile=output_cube_name,
        force_beam=round_beam,
        overwrite=overwrite)
    
def convolve_to_round_beam(
    infile=None, outfile=None, force_beam=None, overwrite=False):
    """
    Convolve supplied image to have a round beam.
    """
    
    if infile is None or outfile is None:
        print("Missing required input.")
        return

    if os.path.isdir(infile) == False:
        print("Input file missing - "+infile)
        return    

    if force_beam is None:
        hdr = imhead(infile)
        if (hdr['axisunits'][0] != 'rad'):
            print("ERROR: Based on CASA experience. I expected units of radians.")
            print("I did not find this. Returning. Adjust code or investigate file "+infile)
            return

        pixel_as = abs(hdr['incr'][0]/np.pi*180.0*3600.)

        if (hdr['restoringbeam']['major']['unit'] != 'arcsec'):
            print("ERROR: Based on CASA experience. I expected units of arcseconds for the beam.")
            print("I did not find this. Returning. Adjust code or investigate file "+infile)
            return    
        bmaj = hdr['restoringbeam']['major']['value']    
        target_bmaj = np.sqrt((bmaj)**2+(2.0*pixel_as)**2)
    else:
        target_bmaj = force_beam

    imsmooth(imagename=infile,
             outfile=outfile,
             targetres=True,
             major=str(target_bmaj)+'arcsec',
             minor=str(target_bmaj)+'arcsec',
             pa='0.0deg',
             overwrite=overwrite
             )

    return target_bmaj

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# ROUTINES FOR FEATHERING THE DATA
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

def prep_for_feather(
    gal=None, array=None, product=None, root_dir=None, 
    overwrite=False):
    """
    Prepare the single dish data for feathering
    """
    
    if gal is None or array is None or product is None or \
            root_dir is None:
        print("Missing required input.")
        return    

    sdfile_in = root_dir+'raw/'+gal+'_tp_'+product+'.image'
    interf_in = root_dir+'process/'+gal+'_'+array+'_'+product+'_flat_round.image'
    pbfile_name = root_dir+'raw/'+gal+'_'+array+'_'+product+'.pb'    

    if (os.path.isdir(sdfile_in) == False):
        print("Single dish file not found: "+sdfile_in)
        return

    if (os.path.isdir(interf_in) == False):
        print("Interferometric file not found: "+interf_in)
        return

    if (os.path.isdir(pbfile_name) == False):
        print("Primary beam file not found: "+pbfile_name)
        return

    # Align the relevant TP data to the product.
    sdfile_out = root_dir+'process/'+gal+'_tp_'+product+'_align_'+array+'.image'
    imregrid(imagename=sdfile_in,
             template=interf_in,
             output=sdfile_out,
             asvelocity=True,
             axes=[-1],
             interpolation='cubic',
             overwrite=overwrite)

    # Taper the TP data by the primary beam.
    taperfile_out = root_dir+'process/'+gal+'_tp_'+product+'_taper_'+array+'.image'
    if overwrite:
        os.system('rm -rf '+taperfile_out)

    impbcor(imagename=sdfile_out, 
            pbimage=pbfile_name, 
            outfile=taperfile_out, 
            mode='multiply',
            stokes='I')

    return

def phangs_feather_data(
    gal=None, array=None, product=None, root_dir=None, 
    cutoff=-1,overwrite=False):
    """
    Feather the interferometric and total power data.
    """

    if gal is None or array is None or product is None or \
            root_dir is None:
        print("Missing required input.")
        return    

    sdfile_in = root_dir+'process/'+gal+'_tp_'+product+'_taper_'+array+'.image'
    interf_in = root_dir+'process/'+gal+'_'+array+'_'+product+'_flat_round.image'
    pbfile_name = root_dir+'raw/'+gal+'_'+array+'_'+product+'.pb' 

    if (os.path.isdir(sdfile_in) == False):
        print("Single dish file not found: "+sdfile_in)
        return
        
    if (os.path.isdir(interf_in) == False):
        print("Interferometric file not found: "+interf_in)
        return

    if (os.path.isdir(pbfile_name) == False):
        print("Primary beam file not found: "+pbfile_name)
        return

    # Feather the inteferometric and "flat" TP data.
    outfile_name = root_dir+'process/'+gal+'_'+array+'+tp_'+product+ \
        '_flat_round.image'

    if overwrite:        
        os.system('rm -rf '+outfile_name)
    os.system('rm -rf '+outfile_name+'.temp')
    feather(imagename=outfile_name+'.temp',
            highres=interf_in,
            lowres=sdfile_in,
            sdfactor=1.0,
            lowpassfiltersd=False)
    imsubimage(imagename=outfile_name+'.temp', outfile=outfile_name,
               dropdeg=True)
    os.system('rm -rf '+outfile_name+'.temp')
    infile_name = outfile_name

    # Primary beam correct the feathered data.
    outfile_name = root_dir+'process/'+gal+'_'+array+'+tp_'+product+ \
        '_pbcorr_round.image'
    
    if overwrite:        
        os.system('rm -rf '+outfile_name)

    print(infile_name)
    print(pbfile_name)
    impbcor(imagename=infile_name,
            pbimage=pbfile_name, 
            outfile=outfile_name, 
            mode='divide', cutoff=cutoff)

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# CLEANUP
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

def convert_jytok(
    infile=None, outfile=None, overwrite=False, inplace=False):
    """
    Convert a cube from Jy/beam to K.
    """

    c = 2.99792458e10
    h = 6.6260755e-27
    kb = 1.380658e-16

    if infile is None or (outfile is None and inplace==False):
        print("Missing required input.")
        return
    
    if os.path.isdir(infile) == False:
        print("Input file not found: "+infile)
        return
    
    if inplace == False:
        if overwrite:
            os.system('rm -rf '+outfile)
        
        if os.path.isdir(outfile):
            print("Output file already present: "+outfile)
            return

        os.system('cp -r '+infile+' '+outfile)
        target_file = outfile
    else:
        target_file = infile

    hdr = imhead(target_file, mode='list')
    unit = hdr['bunit']
    if unit != 'Jy/beam':
        print("Unit is not Jy/beam. Returning.")
        return
    restfreq_hz = hdr['restfreq'][0]
    bmaj_unit = hdr['beammajor']['unit']
    if bmaj_unit != 'arcsec':
        print("Beam unit is not arcsec, which I expected. Returning.")
        print("Unit instead is "+bmaj_unit)
        return    
    bmaj_as = hdr['beammajor']['value']
    bmin_as = hdr['beamminor']['value']
    bmaj_sr = bmaj_as/3600.*np.pi/180.
    bmin_sr = bmin_as/3600.*np.pi/180.
    beam_in_sr = np.pi*(bmaj_sr/2.0*bmin_sr/2.0)/np.log(2)
    jtok = c**2 / beam_in_sr / 1e23 / (2*kb*restfreq_hz**2)

    myia = au.createCasaTool(iatool)
    myia.open(target_file)
    vals = myia.getchunk()
    vals *= jtok
    myia.putchunk(vals)
    myia.setbrightnessunit("K")
    myia.close()

    imhead(target_file, mode='put', hdkey='JTOK', hdvalue=jtok)

    return

def trim_cube(    
    infile=None, outfile=None, overwrite=False, inplace=False, min_pixperbeam=3):
    """
    Trim and rebin a cube to smaller size.
    """
    
    if infile is None or outfile is None:
        print("TRIM_CUBE: Missing required input.")
        return
    
    if os.path.isdir(infile) == False:
        print("TRIM_CUBE: Input file not found: "+infile)
        return

    # First, rebin if needed
    hdr = imhead(infile)
    if (hdr['axisunits'][0] != 'rad'):
        print("ERROR: Based on CASA experience. I expected units of radians.")
        print("I did not find this. Returning. Adjust code or investigate file "+infile)
        return

    pixel_as = abs(hdr['incr'][0]/np.pi*180.0*3600.)

    if (hdr['restoringbeam']['major']['unit'] != 'arcsec'):
        print("ERROR: Based on CASA experience. I expected units of arcseconds for the beam.")
        print("I did not find this. Returning. Adjust code or investigate file "+infile)
        return    
    bmaj = hdr['restoringbeam']['major']['value']    
    
    pix_per_beam = bmaj*1.0 / pixel_as*1.0
    
    if pix_per_beam > 6:
        imrebin(
            imagename=infile,
            outfile=outfile+'.temp',
            factor=[2,2,1],
            crop=True,
            dropdeg=True,
            overwrite=overwrite,
            )
    else:
        os.system('cp -r '+infile+' '+outfile+'.temp')

    # Figure out the extent of the image inside the cube
    myia = au.createCasaTool(iatool)
    myia.open(outfile+'.temp')
    mask = myia.getchunk(getmask=True)    
    myia.close()

    this_shape = mask.shape

    mask_spec_x = np.sum(np.sum(mask*1.0,axis=2),axis=1) > 0
    pad = 0
    xmin = np.max([0,np.min(np.where(mask_spec_x))-pad])
    xmax = np.min([np.max(np.where(mask_spec_x))+pad,mask.shape[0]-1])

    mask_spec_y = np.sum(np.sum(mask*1.0,axis=2),axis=0) > 0
    ymin = np.max([0,np.min(np.where(mask_spec_y))-pad])
    ymax = np.min([np.max(np.where(mask_spec_y))+pad,mask.shape[1]-1])

    mask_spec_z = np.sum(np.sum(mask*1.0,axis=0),axis=0) > 0
    zmin = np.max([0,np.min(np.where(mask_spec_z))-pad])
    zmax = np.min([np.max(np.where(mask_spec_z))+pad,mask.shape[2]-1])
    
    box_string = ''+str(xmin)+','+str(ymin)+','+str(xmax)+','+str(ymax)
    chan_string = ''+str(zmin)+'~'+str(zmax)

    print("... ... ... box selection: "+box_string)
    print("... ... ... channel selection: "+chan_string)

    if overwrite:
        os.system('rm -rf '+outfile)
        imsubimage(
        imagename=outfile+'.temp',
        outfile=outfile,
        box=box_string,
        chans=chan_string,
        )
    
    os.system('rm -rf '+outfile+'.temp')
    
def phangs_cleanup_cubes(
    gal=None, array=None, product=None, root_dir=None, 
    overwrite=False, min_pixeperbeam=3):
    """
    Clean up cubes.
    """

    if gal is None or array is None or product is None or \
            root_dir is None:
        print("Missing required input.")
        return

    for this_ext in ['flat', 'pbcorr']:

        root = root_dir+'process/'+gal+'_'+array+'_'+product+'_'+this_ext
        infile = root+'_round.image'
        outfile = root+'_round_k.image'
        outfile_fits = root+'_round_k.fits'
    
        # Trim the cube to a smaller size and rebin as needed

        trim_cube(infile=infile, outfile=outfile, 
                  overwrite=overwrite, inplace=False,
                  min_pixperbeam=min_pixeperbeam)

        # Convert to Kelvin

        convert_jytok(infile=outfile, inplace=True)

        # Export to FITS
    
        exportfits(imagename=outfile, fitsimage=outfile_fits,
                   velocity=True, overwrite=True, dropstokes=True, 
                   dropdeg=True, bitpix=-32)
    
        # Clean up headers

        pass

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# LINEAR MOSAICKING ROUTINES
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

def common_res_for_mosaic(
    gal=None, array=None, product=None, root_dir=None, 
    overwrite=False, target_res=None):
    """
    Convolve multi-part cubes to a common res for mosaicking.
    """
    
    if gal is None or array is None or product is None or \
            root_dir is None:
        print("Missing required input.")
        return    
    
    # Look up parts
    this_mosaic_key = mosaic_key()
    if (gal in this_mosaic_key.keys()) == False:
        print("Galaxy "+gal+" not in mosaic key.")
        return
    parts = this_mosaic_key[gal]

    # Figure out target resolution if it is not supplied

    if target_res is None:
        print("Calculating target resolution ... ")
        for part in parts:
            print("Checking "+part)
            pass
        
    print("Convolving to target resolution: "+str(target_res))

    for part in parts:
        print("Convolving "+part)
        pass

    return target_res

def align_for_mosaic(
    gal=None, array=None, product=None, root_dir=None, 
    overwrite=False, target_res=None):
    """
    Convolve multi-part cubes to a common res for mosaicking.
    """
    
    if gal is None or array is None or product is None or \
            root_dir is None:
        print("Missing required input.")
        return    
    
    # Look up parts
    this_mosaic_key = mosaic_key()
    if (gal in this_mosaic_key.keys()) == False:
        print("Galaxy "+gal+" not in mosaic key.")
        return
    parts = this_mosaic_key[gal]

def phangs_mosaic_data(
    gal=None, array=None, product=None, root_dir=None, 
    overwrite=False):
    """
    Linearly mosaic multipart cubes.
    """

    if gal is None or array is None or product is None or \
            root_dir is None:
        print("Missing required input.")
        return    
    
    # Look up parts
    this_mosaic_key = mosaic_key()
    if (gal in this_mosaic_key.keys()) == False:
        print("Galaxy "+gal+" not in mosaic key.")
        return
    parts = this_mosaic_key[gal]

    pass
    
