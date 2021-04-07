## README for PHANGS-ALMA Pipeline Version 2.0

### PREFACE

**Contents:** This is "version 2" of the PHANGS post-processing and science-ready data product pipeline. These programs use CASA, astropy, and affiliated packages (analysisutils, spectral-cube, reproject) to process data from the calibrated visibility to science-ready maps. The procedures and background for key parts of the pipeline are discussed in the Astrophysical Journal Supplements Paper "PHANGS-ALMA Data Processing and Pipeline" by Leroy, Hughes, Liu, Pety, Rosolowsky, Saito, Schinnerer, Usero, Faesi, Herrera et al.. Please consult that paper for more background and details.

**Pipeline and Configuration Files:** These are the programs to run the PHANGS-ALMA pipeline. Configuration files for a large set of PHANGS projects, including the live version of the files for the PHANGS-ALMA CO survey, exist in a separate repository. We include a frozen set of files that can be used to reduce PHANGS-ALMA as examples here. If you are need access to those other repositories or need examples, please request access as needed.

**Contact:** For issues, the preferred method is to open an issue on the github issues page. If you have specific other topics to discuss you can email the PHANGS-ALMA data reduction group at adr@phangs.groups.io . If you want to directly contact a person you can reach out to Adam Leroy, Erik Rosolowsky, or Daizhong Liu via email. But issues are better.

**Earlier Versions:** If you are looking for Version 1.0 of the pipeline, you can access it by changing branches to "version1.0". Note that this will mostly be for historical reasons. We suggest using Version 2.0 moving forward.

**Total Power:** The total power scripts and pipeline described by Herrera et al. are here: https://github.com/PhangsTeam/TP_ALMA_data_reduction/ this pipeline assumes that you have already run these to produce cubes.

### TWO WAYS TO USE THE PIPELINE

There are two ways that this pipeline might be useful. First, it provides an end-to-end path to process calibrated ALMA data (or VLA data) of the sort produced by scriptForPI into cubes and maps. That end-to-end approach is described in "Workflow for most users." Second, the `phangsPipeline` directory contains a number of modules for use inside and outside CASA that should have general utility. These are written without requiring any broader awareness of the pipeline infrastructure and should just be generally useful. These are files named `casaSOMENAME.py` and `scSOMEOTHERNAME.py` and, to a lesser extent, `utilsYETANOTHERNAME.py`.

### WORKFLOW FOR MOST USERS

If you just want to *use* the pipeline then you will need to do three things:

( 0. Run scriptForPI to apply the observatory-provided calibration to your data. The pipeline picks up from there, it does not replace the outstanding ALMA observatory calibration and flagging pipeline. )

1. Make configuration files ("key files") that describe your project. Usually you can copy and modify an existing project to get a good start. We provide PHANGS-ALMA as an example.

2. Put together a small script to run the pipeline. Well, really put together two small scripts: one to runt he CASA stuff and another to run the pure python stuff. In theory these could be combined or generalized, but we usually just write a few small programs.

3. Run these scripts in order. The CASA stuff runs inside a CASA shell - the pipeline seems to work up through CASA 5.7 and has been heavily used in 5.4 and 5.6, In theory it should be workable in CASA 6.1+ but this isn't for sure yet. The pure python stuff expects a distribution with numpy, astropy, spectral-cube, and scipy and python 3.6+ or so.

**The Easiest Way** This release includes the full PHANGS-ALMA set of keys and the scripts we use to run the pipeline for PHANGS-ALMA. These are *heavily documented* - copy them to make your own script and configuration and follow the documented in those scripts to get started. To be specific:

The PHANGS-ALMA imaging keys with heavy documentation are in:

The script to run the CASA part of the pipeline is:

The script to create derived procuts is:

### A FEW DETAILS ON PROCEDURE

1. **Staging (in CASA)** Stage and process uv-data. This stage includes continuum subtraction, line extraction, and spectral regridding.

2. **Imaging (in CASA)** Image and deconvolve the uv-data. This runs in several stages: dirty imaging, clean mask alignment, multi-scale deconvolution, re-masking, and single convolution.

3. **Post-Process (in CASA)** Process deconvolved data into science-ready data cubes. This stage includes merging with the total power and mosaicking.

4. **Derive Produts (in python)** Convolution, noise estimation, masking, and calculation of science-ready data products.

The simples way to run these is to write two small scripts and do the following:

1. Initialize CASA
2. Run a script that initializes a `keyHandler` object pointed at your key directory (see below). Then use this keyHandler to initialize handler objects for uv data, imaging, and postprocessing. Optionally restrict thoe objects of interest for each handler to a subset of targets, array configurations, or lines.
3. Inside the same script, run the main loop commands for each handler object.

Then exit CASA and

4. Initialize a python environment with scipy, numpy, astropy, and spectral-cube installed.
5. Run a script that initializes a `keyHandler` again pointed at your key directory, then use this keyHandler to initialize a derived poduct handler.
6. Run the main loop for the derived product handler.

we provide a pair of example scripts that do this for PHANGS-ALMA in the repository: 

### CONTENTS OF THE PIPELINE IN MORE DETAIL

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

