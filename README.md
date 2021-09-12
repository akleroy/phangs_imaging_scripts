## README for PHANGS Pipeline Version 2.0

***If you are looking for Version 1.0 of the pipeline, you can access it by changing branches to "version1.0".***

***These are the programs to run the PHANGS pipeline. Configuration
files for a large set of PHANGS projects, including the full
PHANGS-ALMA CO survey, exist in a separate repository. If you are need
access or more examples, please request access as needed.***

**Contact:** Email leroy.42@osu.edu or and please feel free to open
issues on the github issues page.

This is "version 2" of the PHANGS post-processing and science-ready
data product pipeline. These programs use CASA, astropy, and
affiliated packages (analysisutils, spectral-cube, reproject) to
process data from the calibrated visibility to science-ready maps.

### Notes from Hao He
1. Hao has added [running_steps.md](running_steps.md) to describe the general steps to produce continuum image using PHANGS version 2 scripts.
2. Hao modified the scripts [handlerVis.py](phangsPipeline/handlerVis.py) and [casaVisRoutines.py](phangsPipeline/casaVisRoutines.py) so you can bin the continuum measurement set to a certain channel width by specifying the number of bins.
    1. To do that, you need to set `collpase_cont=True` and `cont_width= [num of bins]` in `this_uvh.loop_stage_uvdata()` function in [run_pipeline_casa.py](run_pipeline_casa.py). 
3. Previously the continuum processing does not produce the appropriate mask. This is because the input image name is wrong for the command to generate the mask. Hao modified scripts  [handlerVis.py](phangsPipeline/handlerVis.py) and [casaVisRoutines.py](phangsPipeline/casaVisRoutines.py) and now it can properly generate the mask for image cleaning.  

### EXECUTIVE SUMMARY

**Procedure** A full pipeline run has four stages:

1. **Staging** Stage and process uv-data. This stage includes
continuum subtraction, line extraction, and regridding.

2. **Imaging** Image and deconvolve the uv-data. This runs in several
stages: dirty imaging, clean mask alignment, multi-scale
deconvolution, re-masking, and single convolution.

3. **Post-Process** Process deconvolved data into science-ready data
cubes. This stage includes merging with the total power and
mosaicking.

4. **Derive Produts** Convolution, noise estimation, masking, and
calculation of science-ready data products.

**Architecture** The pipeline is organized and run by a series of
"handler" objects. These handlers organize the list of targets, array
configurations, spectral products, and derived moments and execute
loops.

The routines to process individual data sets are in individual
modules, grouped by theme (e.g., casaImagingRoutines or
scNoiseRoutines). These routines do not know about the larger
infrastructure of arrays, targets, etc.. They generally take an input
file, output file, and various keyword arguments.

A project is defined by a series of text key files in a
"key_directory". These define the measurement set inputs,
configurations, spectral line products, moments, and derived
products. 

**User Control** For the most part the user's job is to *define the
key files* and to run some scripts.
