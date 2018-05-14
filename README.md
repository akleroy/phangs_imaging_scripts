README for PHANGS Data Reduction
--------------------------------

EXECUTIVE SUMMARY

Download the data. Reduce them using the standard pipeline. This
pipeline will collect them, extract the lines and continuum, combine
measurement sets, image the various products, combine the
interferometer and total power data, and produce products for
distribution.

It runs in four steps right now:

1) Download and calibrate the data using the observatory scripts.
2) Stage the imaging using stage_imaging.py
3) Image the data using image_data.py
4) Create the higher level data products using build_release_v0p6.pro

Parts 1-3 run in CASA. Parts 2-3 depend on phangsPipeline.py and use
the analysisUtils. Part 4 runs in IDL with a stop in CASA to feather.

Email me: leroy.42@osu.edu or otherwise get in touch with questions or
suggestions.

SETUP

This file and the rest of the scripts should sit in a directory called 

imaging/scripts/

under whatever parent directory you choose to use. For the parent directory is

/data/tycho/0/leroy.42/reduction/alma/PHANGS/

and I'll just refer to the top level thing as PHANGS/

within PHANGS/ I have the following directories:

2013.1.00803.S  
2013.1.01161.S  
2015.1.00925.S  
2015.1.00956.S 
2017.1.00886.L
Cycle1_650
imaging
raw

For me, raw/ contains some subdirectories that hold 12m/ 7m/ tp/ data
and I keep the tarballs here until I'm sure that I'm done with
them. This isn't used by the scripts at all. You can do this however
you want.

The rest of the directories DO matter. These hold the uv data and are
referenced in the imaging directory by

ms_file_key.txt

You need ms_file_key.txt to point to the calibrated uv data in order
for the first part of the scripts to work.

You could override the filename in read_ms_key(yourfilename) or you
can just mimic the directory structure above. I suggest the latter.

I keep a list of all the names of data sets and the associated galaxy
and project in SCIENCE_GOAL_KEY.txt - I suggest that you reference
this in setting things up.

SETTING UP THE PIPELINE

The actual imaging part of the pipeline is now mostly concentrated in
the phangsPipeline.py module. You want this imported as "pp" in able
to be able to run the scripts. I do this by adding the following lines
to my ~/.casa/init.py:

sys.path.append("/home/maury/leroy.42/casapy/analysis_scripts/")
sys.path.append("/data/tycho/0/leroy.42/reduction/alma/PHANGS/imaging/scripts/")
import analysisUtils as au
import phangsPipeline as pp

That should get you both the analysis utilities, which you need, and
the phangs pipeline. Change the directories as appropriate, of course.

As long as those two things import successfully, you should be able to
run the pipeline. I use CASA 5.1.1 as of this writing but in general
just try to use the latest CASA.

REDUCTION

After untarring, I move things to the appropriate project directory
under PHANGS/ and then I go to the scripts directory and run the
script for PI. Some tips:

- you may need multiple versions of CASA installed. Already for the
  large program I've needed both 4.7.2 and 5.1. Be sure you have the
  pipeline version.

- when you rerun you will need to delete the calibrated directory.

COMBINING CASES WITH MULTIPLE MEASUREMENT SETS

When multiple measurement sets (.ms files in calibrated/ ) show up for
a single target, we can either reflect these in the ms_file_key as,
e.g., 7m_1 7m_2 and so on OR we can concatenate them in to a single
data set. 

I did this concatenation PRIOR to the large program. To do this, I
copy the script

combineCalibrated.py

from the PHANGS/imaging/scripts/ directory to the scripts/ directory
of whatever project I'm working on. Then I just run that script in
CASA. You should run this whenever you see calibrated_final.ms in the
ms_file_key but not in the calibrated/ directory.

With the start of the large program, we have changed conventions. Now
because various file naming conventions have changed. Moving forward,
we will just start going with 7m_1 7m_2 and so on, with an entry in
the ms_file_key.txt for each observation.

GETTING THE DATA INTO SHAPE FOR IMAGING

After you have the data calibrated, you need to:

- make sure the ms_file_key.txt points at the calibrated data. This
  file allows us to map the calibrated data into the files we will use
  for imaging.

- make sure that you have dir_key.txt in place. This provides a
  mapping of directories for all of our non-standard cases, almost
  always split mosaics. Anything without an entry is taken to be
  one-galaxy-one-science-goal.

- make sure that you have an up-to-date mosaic_definitions.txt . This
  file includes the phase centers and velocities to be used in
  imaging and line extraction.

Now you should be read to run the script

stage_imaging.py

This is a script that you are intended to edit according to your
current needs (aside: this means you should only check it in to github
when you are actually editing the functionality).

When you run stage imaging it will try to do three things:

1) Make the directory for imaging if it doesn't already exist.

2) copy the calibrated data in to the directory where you are going to image.

3) extract lines from all of the individual measurement sets. Then
concatenate them into a single measurement set for each spectral line
on a shared velocity grid to be used for imaging. Also make "channel
0" measurement sets.

4) Extract continuum measurements sets from each individual
measurement set and then concatenate them into a single measurement
set for each spectral line.

You can turn these steps on and off at the top of stage_imaging.

There you can also specify which array you want to stage (7m, 12m or
both). And you can specify which galaxies to skip or a list of which
galaxies to consider.

At the end of this, you should have files that look like this:

gal_7m_co21.ms
gal_7m_co21_chan0.ms
gal_7m_c18o21.ms
gal_7m_c18o21_chan0.ms
gal_7m_cont.ms

As well as a bunch of intermediate files. Those files above will be
used for imaging. Note that some galaxies will lack the c18o21.

IMAGING

After running stage_imaging you are ready to begin imaging. This is
done using the script

image_data.py

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

CLEAN MASKS

The multiscale clean involves a broad clean mask. So far, I have been
creating these as an output of the last step of the pipeline (data
product creation) based on the feathered cubes with some heavy
processing. Then these are placed back in the clean mask directory and
imaging gets rerun. This means that in principle imaging is an
iterative process. In practice this mostly just takes one iteration.

I'll aim to add some testplot creation to show residuals and the
current mask and check that in.

DATA PRODUCT CREATION

Right now data product creation occurs using an IDL pipeline. The
latest version of this as of this writing is build_release_v0p6. But
I'll revise this after revising the imaging part of the pipeline. So
expect a new version soon, but this currently should run more or less
fine once a few changed naming conventions are resolved.

