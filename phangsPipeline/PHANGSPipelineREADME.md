## README for the PHANGS Pipeline

This README specifically describes the PHANGS CASA and python
pipeline. 

Contact: Email leroy.42@osu.edu or otherwise get in touch with
questions or suggestions.

### REQUIRED SETUP

You need the analysisUtils to use the ppieline. I import these
automatically by adding the following lines to my `~/.casa/init.py`:

```
sys.path.append("/home/maury/leroy.42/casapy/analysis_scripts/")
import analysisUtils as au
```

### OVERVIEW OF DATA AND DIRECTORY MODEL

The pipeline assumes the following structure:

working directory/
├── scripts/
│   ├── stage_imaging.py
│   ├── image_data.py
│   ├── process_cubes.py
│   ├── phangsPipeline/
│   └── data_files/
├── target1/
├── target2/
├── target3/
├── ...
├── targetn/
├── release/
│   └── vX/
│   	├── raw/
│   	├── process/
│   	├── products/
│   	└── delivery/
└── singledish*/

* Optional - the single data can in theory live anywhere.

This is not exhaustive, but gives the idea. Key points:

* All operations assume a relative structure with scripts/ parallel to
  the directories for individual targets or groups of targets.

* The pipeline code, which the user does not have to modify, lives in
  a subdirectory of scripts/ called phangsPipeline/ .

* The data files that define targets, uv data, single dish data,
  etc. for individual projects live in a subdirectory of scripts/
  called called data_files/ . A separate set of data_files can be used
  to create a separate application of the pipeline for different
  projects, e.g., the PHANGS-ALMA CO 2-1 surveys use one set of files,
  while the dense gas mapping uses another set of files.

* Each source ("target") has a working directory in parallel to the
  scripts/ directory. The pipeline changes to that directory and works
  there when working on that target. Sometimes multiple targets are
  grouped together in a single directory. For example, this happens if
  a single galaxy gets observed across multiple independent
  mosaics. This is controlled by a data file called dir_key.txt and
  entirely up to the user.

* The user manipulates stage_imaging, image_data, etc. in the scripts/
  to run the staging, imaging, etc.. This is the main user interface
  to the pipeline.

* The process_cubes and following scripts (some still IDL) copy the
  imaged data to the release/vX/ (where X is the release name)
  directory and then post-process the imaging in these directories.

* A delivery for distribution is built into the release/vX/delivery/
  directory.

### COMMENTARY ON THE ASSUMED DATA MODEL

The pipeline is currently built to operate on *one line at a time*,
because it assumes the ability to map between velocity and
frequency. It can operate on many different lines in serial, e.g.,
CO2-1 then 13CO2-1 then C18O2-1, with the lines defined in 'line_list'
with frequencies taken from splatalogue.

The pipeline makes some assumptions that are worth spelling out:

* Right now a single measurement set is assumed to contain a single
  target. The properties of each target are defined in a data
  file. 

  This could be modified by adding a field selection to the
  "ms_file_key", but we do not currently have any data that require
  it.

* Right now the pipeline distinguishes between the 12m array and the
  7m array, but doesn't distinguish beyond that. Individual
  measurement sets get the label, e.g., 7m_1, 7m_2, etc.. These are
  grouped and processed together. 

  Starting in Cycle 8, we will beging to combine multiple 12m
  configurations and we should revise the formalism a little. It will
  make sense to assign each measurement set to a user defined array
  (e.g., "12m_ext") and separate the labeling from the array code.

### PIPELINE CONTENTS

The pipeline is organized as a package called phangsPipeline. Within
that there are the following modules.

* keyHandler : this handles the data files that describe the location
  of uv data, properties of targets, and more or less all the other
  user-defined input to the pipeline. A KeyHandler can be intiated
  pointing at different files to run the pipeline for different
  projects.


