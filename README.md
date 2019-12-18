## README for PHANGS Data Reduction

Email leroy.42@osu.edu or otherwise get in touch with questions or
suggestions.

### EXECUTIVE SUMMARY

Download the data. Reduce them using the standard pipeline. This
pipeline will collect them, extract the lines and continuum, combine
measurement sets, image the various products, combine the
interferometer and total power data, and produce products for
distribution.

It runs in five steps right now:

1) Download and calibrate the data using the observatory scripts.

2) Stage the imaging using `stage_imaging.py`

3) Image the data using `image_data.py`

4) Create final cubes using one of two methods:

   (a) `IDL build_cubes.pro`

   (b) `process_cubes.py`

   These steps convolve to a round beam, primary beam correct,
feather, combine multi-part galaxies, and then export final FITS
cubes. Eventually the two steps will be equivalent and we'll add a
third version operating purely using spectral cube and astropy. For
now the IDL pipeline is most complete.

5) Create higher level data products using the IDL scripts
`build_products_7m.pro`, `build_products_12m.pro`. These steps build
masks, create fixed resolution cubes, and build moment maps and
similar products.

6) compile cubes into a release using the IDL build_release.

Parts 1-3 and 4b run in CASA. Parts 2-4 depend on `phangsPipeline.py`
and use the `analysisUtils` available from [here](ftp://ftp.cv.nrao.edu/pub/casaguides/analysis_scripts.tar).

Part 4 has been shifted to prefer CASA.

Parts 5 and 6 currently run in IDL.

### SETUP

This file and the rest of the scripts should sit in a directory called 

```
imaging/scripts/
```

under whatever parent directory you choose to use. For the parent directory is

```
/data/tycho/0/leroy.42/reduction/alma/PHANGS/
```

and I'll just refer to the top level thing as `PHANGS/`

within PHANGS/ and a couple other directories also checked by the
scripts I have the following directories holding calibrated data:

```
2013.1.00803.S
2013.1.01161.S
2015.1.00782.S
2015.1.00925.S
2015.1.00956.S
2017.1.00392.S
2017.1.00766.S
2017.1.00886.L
2018.1.01321.S
2018.1.01651.S
2018.A.00062.S
Cycle1_650
imaging
```

These hold the uv data and are referenced in the imaging directory by `ms_file_key.txt`.

You need `ms_file_key.txt` to point to the calibrated uv data in order
for the first part of the scripts to work.

Note that in `stage_imaging.py` you can specify a root directory to
search. `ms_file_key.txt` only needs to work relative to the roots that
you supply in `stage_imaging.py`.

You can hack the scripts to get around this, but the easiest thing is
to just follow this approach. The option to specify the root directory
in stage_imaging.py makes it pretty easy to do this.

### SETTING UP THE PIPELINE

The actual imaging part of the pipeline is now mostly concentrated in
the phangsPipeline.py module. You want this imported as `pp` in able
to be able to run the scripts. I do this by adding the following lines
to my `~/.casa/init.py`:

```
sys.path.append("/home/maury/leroy.42/casapy/analysis_scripts/")
sys.path.append("/data/tycho/0/leroy.42/reduction/alma/PHANGS/imaging/scripts/")
import analysisUtils as au
import phangsPipeline as pp
```

That should get you both the analysis utilities, which you need, and
the phangs pipeline. Change the directories as appropriate, of course.

As long as those two things import successfully, you should be able to
run the pipeline. I use CASA 5.4.0 for imaging and staging as of this
writing but in general just try to use the latest CASA. I'm using
5.6.1 for `process_cubes`.

Note that there are two other big pipeline
modules. `phangsPipelinePython.py` has file accessing parts of the
pipeline that don't depend on CASA, and so can be imported
elsewhere. `phangsCubePipeline.py` contains routines for dealing with
images and post processing.

### REDUCTION

After untarring, I move things to the appropriate project directory
under `PHANGS/` and then I go to the scripts directory and run the
script for PI. Some tips:

- You may need multiple versions of CASA installed. For the large
  program I've needed both 4.7.2 and 5.1 and for later versions I've
  needed 5.4.0 and 5.6.1. Be sure you have the pipeline version.

- When you rerun remember that you will need to delete the calibrated
  directory.

- If you are adding new data, remember that you will need run the
  scripts (usually just the calibration application in `scriptForPI.py`)
  then look in `../calibrated/` and add the new measurement sets to
  the `ms_file_key.txt`. After this, the pipeline is ready to ingest
  the data into imaging.

### COMBINING CASES WITH MULTIPLE MEASUREMENT SETS

When multiple measurement sets (.ms files in `calibrated/` ) show up for
a single target, we can either reflect these in the `ms_file_key.txt` as,
e.g., `7m_1` `7m_2` and so on OR we can concatenate them in to a single
data set.

I did this concatenation PRIOR to the large program (but do not do it
for later projects). To do this, I had been copying the script
`combineCalibrated.py` from the `PHANGS/imaging/scripts/` directory to
the `scripts/` directory of whatever project I'm working on. Then I just
run that script in CASA. You should run this whenever you see
`calibrated_final.ms` in the `ms_file_key` but not in the `calibrated/`
directory.

With the start of the large program, we have changed conventions. Now
because various file naming conventions have changed. Moving forward,
we will just start going with `7m_1`, `7m_2`, `12m_1`, and so on. We
put an entry in the `ms_file_key.txt` for each observation.

### GETTING THE DATA INTO SHAPE FOR IMAGING

After you have the data calibrated, you need to:

- Make sure the `ms_file_key.txt` points at the calibrated data. This
  file allows us to map the calibrated data into the files we will use
  for imaging.

- Make sure that you have `dir_key.txt` in place. This provides a
  mapping of directories for all of our non-standard cases. These are
  almost always split mosaics, where we image galaxies in multiple
  parts. Anything without an entry is taken to be
  one-galaxy-one-science-goal.

- Make sure that you have an up-to-date `mosaic_definitions.txt`. This
  file includes the phase centers and velocities to be used in
  imaging and line extraction.

Now you should be read to run the script `stage_imaging.py`

This is a script that you are intended to edit according to your
current needs (aside: this means you should only check it in to github
when you are actually editing the functionality).

When you run stage_imaging it will try to do three things:

1) Make the directory for imaging if it doesn't already exist.

2) copy the calibrated data in to the directory where you are going to image.

3) extract lines from all of the individual measurement sets. Then
concatenate them into a single measurement set for each spectral line
on a shared velocity grid to be used for imaging. Also make "channel
0" measurement sets.

4) Extract continuum measurements sets from each individual
measurement set and then concatenate them into a single measurement
set for each spectral line.

You can turn these steps on and off at the top of stage_imaging by
flipping the boolean flags on and off.

There you can also specify which array you want to stage (7m, 12m or
both). And you can specify which galaxies to skip or a list of which
galaxies to consider.

At the end of this, you should have files that look like this:

```
gal_7m_co21.ms
gal_7m_co21_chan0.ms
gal_7m_c18o21.ms
gal_7m_c18o21_chan0.ms
gal_7m_cont.ms
```

As well as a bunch of intermediate files. Those files above will be
used for imaging. Note that some galaxies will lack the c18o21.

Note that the way we do this it that the velocity gridding is done at
this stage. This might change in the future. But for right now, you
would need to restage the data at a higher velocity resolution.

### IMAGING

After running stage_imaging you are ready to begin imaging. This is
done using the script `image_data.py`.

Similar to stage_data, the idea here is that you edit the top of the
script to select which array, galaxy, and data product you would like
to image.

The imaging itself is current set up in these stages:

1) Dirty map creation.

2) Multiscale clean down to S/N ~ 4 with a broad mask or no mask.

3) Single scale clean of the remaining flux down to low S/N or
convergence in flux. Uses a narrower mask.

4) Export to FITS files.

These steps are carried out inside the phangsPipeline, which in
principle (with a bit of work) can be deployed in a variety of more
general ways.

### CLEAN MASKS

The multiscale clean involves a broad clean mask. So far, I have been
creating these as an output of the last step of the pipeline (data
product creation) based on the feathered cubes with some heavy
processing. Then these are placed back in the clean mask directory and
imaging gets rerun. This means that imaging is an iterative
process. In practice this mostly just takes one iteration: image
without a clean mask, build a clean mask, then image again.

### POST-PROCESSING

After imaging, you have one of two options to process the cubes into a
final form: either use the IDL script `build_cubes.pro` or use the
(under development) CASA script `process_cubes.py`. Both apply primary
beam corrections, convolve the data to have a round beam, feather the
data with the single dish data, and export to a final cube.

`process_cubes.py` works just like `stage_imaging.py` and
`image_data.py`. The IDL side of the pipeline is better vetted, and
works similarly. You just need to call build_cubes with the various
switches flipped.

I'll add more on this and on product creation as the pipeline
approaches public release.


