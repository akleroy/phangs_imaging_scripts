# Script to process cubes produced by clean into FITS files useful for
# future analysis.

tested_versions = ['4.6.0','4.7.0']
this_version = (casa['build']['version']).split('-')[0]
if this_version not in tested_versions:
    print "The script hasn't been verified for this version of CASA."
    print "This version of CASA is "+this_version
    print "Tested versions are "+str(tested_versions)

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# CHECK INPUTS 
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

abort = False
print "---------- statCleanCube.py ----------"

# ......................................
# Input/output files
# ......................................

try:
    cube_root
except NameError:
    print "Please define a root output name for the cube via cube_root=XXX."
    abort = True

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# STAT THE MODEL AND RESIDUAL
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
    
print "... running statistics on model and residual cubes."

imstat_model = imstat(cube_root+'.model')

imstat_residual = imstat(cube_root+'.residual')


