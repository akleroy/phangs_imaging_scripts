tested_versions = ['4.6.0','4.7.0','4.7.1']
this_version = (casa['build']['version']).split('-')[0]
if this_version not in tested_versions:
    print "The script hasn't been verified for this version of CASA."
    print "This version of CASA is "+this_version
    print "Tested versions are "+str(tested_versions)

import glob
import numpy as np
import scipy

# Hate this ... but directories
execfile('../scripts/auExtensions.py')

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# TIMER
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

import time
total_start_time = time.time()

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# INPUTS
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

print ""
print "....................................."
print "alignMask: begins."

# ......................................
# Input/output files
# ......................................

try:
    cube_root
except NameError:
    print "alignMask: Please define a root output name via cube_root=XXX."
    abort = True

try:
    mask_in
except NameError:
    print "alignMask: Please define a mask input file name via mask_in=XXX."
    abort = True

try:
    mask_root
except NameError:
    print "alignMask: Please define a mask output file name via mask_root=XXX."
    abort = True

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# IMPORT TO FITS
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

os.system('rm -rf '+mask_root+'_temp.image')
importfits(fitsimage=mask_in, imagename=mask_root+'_temp.image', overwrite=True)

os.system('rm -rf '+mask_root+'_temp_align.image')
imregrid(imagename=mask_root+'_temp.image', 
         template=cube_root+'.image', 
         output=mask_root+'_temp_align.image', 
         asvelocity=True,
         interpolation='nearest',         
         replicate=False,
         overwrite=True)

im_head = imhead(cube_root+'.image')
mask_head = imhead(mask_root+'_temp_align.image')

os.system('rm -rf '+mask_root+'.mask')
os.system('cp -r '+cube_root+'.image '+mask_root+'.mask')

if (im_head['axisnames'][3] == 'Frequency') and \
        (im_head['ndim'] == 4):    
    ia.open(mask_root+'_temp_align.image')
    mask = ia.getchunk(dropdeg=True)
    ia.close()

    ia.open(mask_root+'.mask')
    data = ia.getchunk(dropdeg=False)
    data[:,:,0,:] = mask
    ia.putchunk(data)
    ia.close()
elif (im_head['axisnames'][2] == 'Frequency') and \
        (im_head['ndim'] == 4):    
    ia.open(mask_root+'_temp_align.image')
    mask = ia.getchunk(dropdeg=True)
    ia.close()

    ia.open(mask_root+'.mask')
    data = ia.getchunk(dropdeg=False)
    data[:,:,:,0] = mask
    ia.putchunk(data)
    ia.close()
else:
    print "ALERT! Did not find a case."

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# CLEAN UP
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

print "alignMask: removing temporary files."

os.system('rm -rf '+mask_root+'_temp.image')
os.system('rm -rf '+mask_root+'_temp_align.image')

print "alignMask: ends."
print "....................................."
