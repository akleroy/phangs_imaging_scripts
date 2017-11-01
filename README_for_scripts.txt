README for PHANGS Data Reduction
--------------------------------

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

I have been doing the latter (concatenation) for the cases PRIOR to
the large program. To do this, I copy the script

combineCalibrated.py

from the PHANGS/imaging/scripts/ directory to the scripts/ directory
of whatever project I'm working on. Then I just run that script in
CASA. You should run this whenever you see calibrated_final.ms in the
ms_file_key but not in the calibrated/ directory.

Note that we are going to change conventions for the large program
because various file naming conventions have changed. Moving forward,
we will just start going with 7m_1 7m_2 and so on.

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

When you run stage imaging it will try to do four things:

1) copy the calibrated data in to the directory where you are going to image.

2) run any custom scripts for that individual galaxy. This isn't
implemented right now, but will eventually allow handling of things
like post-calibration flagging or (e.g., for 5128) continuum
subtraction or similar stuff.

3) extract lines from all of the individual measurement sets and then
concatenate them into a single measurement set for each spectral line.

4) extract continuum measurements sets from each individual
measurement set and then concatenate them into a single measurement
set for each spectral line.

You can turn these steps on and off at the top of stage_imaging.

There you can also specify which array you want to stage (7m, 12m or
both). And you can specify which galaxies to skip or a list of which
galaxies to consider.

At the end of this, you should have files that look like this:

gal_co21.ms
gal_co21_chan0.ms
gal_c18o21.ms
gal_c18o21_chan0.ms
gal_cont.ms

As well as a bunch of intermediate files. Those files above will be
used for imaging. Note that some galaxies will lack the c18o21.

IMAGING

Right now this still uses the script infrarestructure
(image_all_data.py etc.) - I'm working on updating this to 
