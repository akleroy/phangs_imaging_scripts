#!/usr/bin/env python
# 
# Run this script OUTSIDE CASA in an environment that has astropy,
# spectral-cube, scipy, and numpy.
# 

##############################################################################
# Load routines, change directory, initialize handlers
##############################################################################

import os, sys
import importlib

# Pipeline directory. Set this to the location on your system

pipedir = '/data/tycho/0/leroy.42/reduction/alma/phangs_imaging_scripts/'

# Location of the master key. Set this to the master key that points
# to all of the keys for your project.

key_file = '/data/tycho/0/leroy.42/reduction/alma/phangs_imaging_scripts/phangs-alma_keys/master_key.txt'

# Change directory to the pipeline directory.

os.chdir(pipedir)

# Import the logger and initialize the logging. You can change the
# level of message that you want to see by changing "level" here or
# save to a logfile with the keyword.

from phangsPipeline import phangsLogger as pl

# reloads are useful for debugging but can be commented out
importlib.reload(pl)
pl.setup_logger(level='DEBUG', logfile=None)

# Imports

# sys.path.insert(1, )
from phangsPipeline import handlerKeys as kh
from phangsPipeline import handlerDerived as der

# reloads are useful for debugging but can be commented out
importlib.reload(kh)
importlib.reload(der)

# Initialize the various handler objects. First initialize the
# KeyHandler, which reads the master key and the files linked in the
# master key. Then feed this keyHandler, which has all the project
# data, into the other handlers (here DerivedHandler), which run the
# actual pipeline using the project definitions from the KeyHandler.

this_kh = kh.KeyHandler(master_key=key_file)
this_der = der.DerivedHandler(key_handler=this_kh)

# Make missing directories

this_kh.make_missing_directories(imaging=True, derived=True, postprocess=True, release=True)

##############################################################################
# Set up what we do this run
##############################################################################

# Set the configs (arrays), spectral products (lines), and targets to
# consider.

# Set the targets. Called with only () it will use all targets. The
# only= , just= , start= , stop= criteria allow one to build a smaller
# list. 

# Set the configs. Set both interf_configs and feather_configs just to
# determine which cubes will be processed. The only effect in this
# derive product calculation is to determin which cubes get fed into
# the calculation.

# Set the line products. Similarly, this just determines which cubes
# are fed in. Right now there's no derived product pipeline focused on
# continuum maps.

# ASIDE: In PHANGS-ALMA we ran a cheap parallelization by running
# several scripts with different start and stop values in parallel. If
# you are running a big batch of jobs you might consider scripting
# something similar.

this_der.set_targets()
# this_der.set_targets(only=['ngc1809'])

this_der.set_interf_configs(only=['12m+7m'])
this_der.set_feather_configs(only=['12m+7m+tp'])

this_der.set_line_products(only=['co21'])
this_der.set_no_cont_products(True)

# Use boolean flags to set the steps to be performed when the pipeline
# is called. See descriptions below (but only edit here).

do_convolve = True
do_noise = True
do_strictmask = True
do_broadmask = True
do_moments = True
do_secondary = True

##############################################################################
# Step through derived product creation
##############################################################################

# Run the calculations requested by the user. The steps are annotated
# here, but in general, do not change anything below this line. Just
# use the flags above to steer the calculation.

# Convolve the post-processed data products to the various angular and
# physical resolutions specified in the keys.

if do_convolve:
    this_der.loop_derive_products(do_convolve=True, do_noise=False,
                                  do_strictmask=False, do_broadmask=False,
                                  do_moments=False, do_secondary=False)

# Estimate the noise from the signal-free regions of the data to
# produce a three-dimensional noise model for each cube.

if do_noise:
    this_der.loop_derive_products(do_convolve=False, do_noise=True,
                                  do_strictmask=False, do_broadmask=False,
                                  do_moments=False, do_secondary=False)

# Construct "strict masks" for each cube at each resolution.

if do_strictmask:
    this_der.loop_derive_products(do_convolve=False, do_noise=False,
                                  do_strictmask=True, do_broadmask=False,
                                  do_moments=False, do_secondary=False)

# Combine the strict masks across all linked resolutions to form
# "broad masks" that have high completeness.

if do_broadmask:
    this_der.loop_derive_products(do_convolve=False, do_noise=False,
                                  do_strictmask=False, do_broadmask=True,
                                  do_moments=False, do_secondary=False)

# Apply the masks and use the cubes and noise models to produce moment
# maps with associated uncertainty.

if do_moments:
    this_der.loop_derive_products(do_convolve=False, do_noise=False,
                                  do_strictmask=False, do_broadmask=False,
                                  do_moments=True, do_secondary=False)

# Run a second round of moment calculations. This enables claculation
# of moments that depend on other, earlier moment map calculations

if do_secondary:
    this_der.loop_derive_products(do_convolve=False, do_noise=False,
                                  do_strictmask=False, do_broadmask=False,
                                  do_moments=False, do_secondary=True)
