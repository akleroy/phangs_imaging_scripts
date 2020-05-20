## README for PHANGS Pipeline Version 2.0

*** If you are looking for Version 1.0 of the pipeline, you can access
    it by changing branches to "version1.0".***

**Contact:** Email leroy.42@osu.edu or and please feel free to open
issues on the github issues page.

This is "version 2" of the PHANGS post-processing and science-ready
data product pipeline. These programs use CASA, astropy, and
affiliated packages (analysisutils, spectral-cube, reproject) to
process data from the calibrated visibility to science-ready maps.

### EXECUTIVE SUMMARY

***Procedure*** Pipeline processing has four stages:

1. ***Staging*** Stage and process uv-data.

2. ***Imaging*** Image and deconvolve the uv-data.

3. ***Post-Process*** Process deconvolved data into science-ready data cubes.

4. ***Derive Produts*** Convolution, noise estimation, masking, and calculation of science-ready data products.

The pipeline is organized and run by a series of "handler"
objects. These handlers organize the list of targets, array
configurations, spectral products, and derived moments and execute loops.

