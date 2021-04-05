## README for PHANGS-ALMA Pipeline Version 2.0

***If you are looking for Version 1.0 of the pipeline, you can access it by changing branches to "version1.0". Note that this will mostly be for historical reasons. We suggest using Version 2.0 moving forward.***

***These are the programs to run the PHANGS-ALMA pipeline. Configuration files for a large set of PHANGS projects, including the live version of the files for the PHANGS-ALMA CO survey, exist in a separate repository. We include a frozen set of files that can be used to reduce PHANGS-ALMA as examples here. If you are need access to those other repositories or need examples, please request access as needed.***

**Contact:** For issues, the preferred method is to open an issue on the github issues page. If you have specific other topics to discuss you can email the PHANGS-ALMA data reduction group at mailto:adr@phangs.groups.io . If you want to directly contact a person you can reach out to Adam Leroy, Erik Rosolowsky, or Daizhong Liu via email. But issues are better.

This is "version 2" of the PHANGS post-processing and science-ready data product pipeline. These programs use CASA, astropy, and affiliated packages (analysisutils, spectral-cube, reproject) to process data from the calibrated visibility to science-ready maps. The procedures and background for key parts of the pipeline are discussed in the Astrophysical Journal Supplements Paper "PHANGS-ALMA Data Processing and Pipeline" by Leroy, Hughes, Liu, Pety, Rosolowsky, Saito, Schinnerer, Usero, Faesi, Herrera et al.. Please consult that paper for more background and details.

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

