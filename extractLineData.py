# This script extracts a line measurement set from calibrated
# visibility data sets. The idea is to yield a single, simple
# visibility set that can be fed to our imaging scripts.

# USAGE: Run *once* with "do_copy=True" and "do_concat=True" for each
# data set. Then run with "do_extract=True" once *for each line.*

# do_copy:
# - The user supplies a list of measurement sets. 
# - The script copies these to a unified naming scheme.

# do_concat:
# - Data are optionally smoothed in time.
# - All data are concatenated into a single measurement set.
# - Intermediate products are deleted.

# do_extract:
# - Data are regridded to the desired velocity grid.
# - A channel zero data set is created.

# TO DO:
#
# - For maximum generality, this script should take the central
# velocity and a width and do a lot of the supporting calculations,
# saving the user a headache.
#
# - Would also then make sense to start with known lines. Later.
#

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# CHECKS AND DEFAULTS
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

tested_versions = ['4.6.0','4.7.0']
this_version = casa['build']['version']
if this_version not in tested_versions:
    print "The script hasn't been verified for this version of CASA."
    print "This version of CASA is "+this_version
    print "Tested versions are "+str(tested_versions)
else:
    print "The script has been verified for this version of CASA."

sol_kms = 2.99e5

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# TIMER
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

import time
start_time = time.time()

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# CONTROL FLOW
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

try:
    do_copy
except NameError:    
    do_copy = True

try:
    do_concat
except NameError:    
    do_concat = True

try:
    do_extract
except NameError:    
    do_extract = True

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# INPUTS
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

print "---------- extractLineData.py ----------"

try:
    calibrated_files
except NameError:
    print "Please define a calibrated_files file dictionary."
    print "Format is calibrated_files[tag]:ms_to_copy."
    print "(I will turn off all steps for this round)."
    do_copy = False
    do_concat = False
    do_extract = False

try:
    fields_to_use
except NameError:
    print "Please define a fields dictionary with the same keys as calibrate_files."
    print "Format is calibrated_files[tag]:'fields_to_use'."
    print "Can be an empty string for all fields."
    print "(I will turn off all steps for this round)."
    do_copy = False
    do_concat = False
    do_extract = False

try:
    timebin
except NameError:
    print "No timebin specified. I will leave this blank."
    timebin = '0s'

try:
    do_statwt
except NameError:
    print "I *will* use statwt to reset the weights set do_statwt=False to turn this off."
    do_statwt = True
 
try:
    out_root
except NameError:
    print "Please define a root for the output file names. Use out_root="
    print "Our current practice is to set this to a source identifier."
    print "(I will turn off all steps for this round)."
    do_copy = False
    do_concat = False
    do_extract = False

try:
    tag
except NameError:
    print "Please define a tag for the output file names. Use tag="
    print "Our current practice is to set this to a line identifier." 
    print "(I will turn off all steps for this round)."
    do_copy = False
    do_concat = False
    do_extract = False

# .............................................................
# Define the line to extract and averaging to do 
# .............................................................

# Eventually we can imagine automating aspects of this, including
# using analysisUtils to figure out the ideal channel averaging
# calls. For now, let's get a working path and then automate.

try:
    linetag
except NameError:
    print "I will use a default linetag of 'co21'."
    linetag = 'co21'

try:
    restfreq_ghz
except NameError:
    print "I will assume that we target CO 2-1 and set restfreq_ghz=230.538."
    restfreq_ghz = 230.53800

try:
    chan_dv_kms
except NameError:
    print "I will assume a desired channel width of 2.5kms."
    chan_dv_kms = 2.5

try:
    source_vel_kms
except NameError:
    print "I will assume a desired starting velocity of 0km/s."
    source_vel_kms = '0km/s'

try:
    vwidth_kms
except NameError:
    print "I will assume a desire velocity width of 500km/s."
    vwidth_kms = 500

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# COPY DATA FROM ITS ORIGINAL LOCATION
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

if do_copy:

    print "--------------------------------------------------------"
    print "extractLineData: (1) Copying data from original location"
    print "--------------------------------------------------------"

    for this_tag in calibrated_files.keys():
        in_file = calibrated_files[this_tag]
        out_file = out_root+'_'+this_tag+'.ms'
        os.system('rm -rf '+out_file)
        os.system('rm -rf '+out_file+'.flagversions')
        os.system('scp -r '+in_file+' '+out_file)

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# TIME AVERAGE AND CONCATENATE
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

if do_concat:

    print "--------------------------------------------------------"
    print "extractLineData: (2) Processing line into single MS"
    print "--------------------------------------------------------"
      
    print "... time smoothing and then concatenating all of the data."

    files_to_concat = []

    for this_tag in calibrated_files.keys():
        this_infile = out_root+'_'+this_tag+'.ms'
        this_outfile = out_root+'_'+this_tag+'_temp_timebin.ms'
        
        os.system('rm -rf '+this_outfile)
        os.system('rm -rf '+this_outfile+'.flagversions')

        split(vis=this_infile
              , field=fields_to_use[this_tag]
              , datacolumn='DATA'
              , timebin=timebin
              , outputvis=this_outfile)        

        files_to_concat.append(this_outfile)

    concat_out_file = out_root+'_'+tag+'_concat.ms'
    os.system('rm -rf '+concat_out_file)
    os.system('rm -rf '+concat_out_file+'.flagversions')
    concat(vis=files_to_concat,
           concatvis=concat_out_file)

    # ......................................
    # Delete intermediate files
    # ......................................

    print "... cleaning up copied files."

    for this_tag in calibrated_files.keys():
        copied_file = out_root+'_'+this_tag+'.ms'
        os.system('rm -rf '+copied_file)
        os.system('rm -rf '+copied_file+'.flagversions')

    print "... cleaning up time-smoothed individual files."

    for file_name in files_to_concat:
        os.system('rm -rf '+file_name)
        os.system('rm -rf '+file_name+'.flagversions')    

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# EXTRACT THE LINE OF INTEREST
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
  
if do_extract:

    print "--------------------------------------------------------"
    print "extractLineData: (3) Extract the line of interest"
    print "--------------------------------------------------------"

    print "... extracting line data using MSTRANSFORM."

    # .........................................
    # Calculate and recast as strings
    # .........................................

    # Crude, but probably fine for extragalactic work. Figure out what
    # SPWs contain the line in question.

    target_freq_ghz = restfreq_ghz*(1.-source_vel_kms/sol_kms)

    spw_list = au.getScienceSpwsForFrequency(concat_out_file,
                                             target_freq_ghz*1e9)    
    spw_list_string = ''    
    first = True
    for spw in spw_list:
        if not first:
            spw_list_string += ','
        first = False
        spw_list_string += str(spw)

    # Figure out the starting velocity and number of channels.
    start_vel_kms = (source_vel_kms - vwidth_kms/2.0)
    nchan = int(floor(vwidth_kms / chan_dv_kms)+1)

    # I'm hardcoding the precision with extragalactic ALMA work in mind.
    restfreq_string = "{:10.6f}".format(restfreq_ghz)+'GHz'
    chan_dv_string =  "{:5.2f}".format(chan_dv_kms)+'km/s'
    start_vel_string =  "{:6.1f}".format(start_vel_kms)+'km/s'

    # Figure out how much averaging is needed to reach
    chan_width_hz = au.getChanWidths(concat_out_file, spw_list_string)
    target_width_hz = chan_dv_kms/sol_kms*restfreq_ghz*1e9
    rebin_factor = target_width_hz / chan_width_hz
    
    if rebin_factor < 2:
        chanbin = 1
    else:
        chanbin = int(floor(rebin_factor/2.))

    # .........................................
    # Concert to strings and summarize
    # .........................................

    print "LINE TAG: ", linetag
    print "REST FREQUENCY: ", restfreq_string
    print "SPECTRAL WINDOWS: ", spw_list
    print "TARGET CHANNEL WIDTH: ", chan_dv_string
    print "START VELOCITY: ", start_vel_string
    print "NUMBER OF CHANNELS: ", str(nchan)
    print "CHANNELS TO BIN TOGETHER FIRST: ", chanbin
    
    # .........................................
    # Regrid
    # .........................................

    concat_out_file = out_root+'_'+tag+'_concat.ms'

    print "... carrying out channel averaging first."

    final_out_file = out_root+'_'+linetag+'_regrid.ms'

    os.system('rm -rf '+final_out_file+'.temp')
    os.system('rm -rf '+final_out_file+'.temp.flagversions')

    if chanbin > 1:
        chanaverage = True
    else:
        chanaverage = False

    mstransform(vis=concat_out_file,
                spw=spw_list_string,
                outputvis=final_out_file+'.temp',
                datacolumn='DATA',
                chanaverage=chanaverage,
                chanbin=chanbin,
                hanning=True,
                # Does this matter? It shouldn't, but AS thinks it still might
                interpolation='cubic',
                )

    print "... now regridding to the desired velocity grid."

    os.system('rm -rf '+final_out_file)
    os.system('rm -rf '+final_out_file+'.flagversions')

    mstransform(vis=final_out_file+'.temp',
                outputvis=final_out_file,
                datacolumn='DATA',
                combinespws=True,
                regridms=True,
                mode='velocity',
                interpolation='cubic',
                start=start_vel_string,
                nchan=nchan,
                width=chan_dv_string,
                restfreq=restfreq_string,
                outframe='lsrk',
                veltype='radio',
            )

    os.system('rm -rf '+final_out_file+'.temp')
    os.system('rm -rf '+final_out_file+'.temp.flagversions')

    # .........................................
    # Use "statwt" to derive empirical weights
    # .........................................

    if do_statwt:

        print "... deriving empirical weights using STATWT."

        statwt(vis=final_out_file,
               datacolumn='DATA')

    # ......................................
    # Make "Channel 0" CO measurement sets
    # ......................................

    print '... making a "channel zero" data set using SPLIT.'

    chan0_vis = out_root+'_'+linetag+'_chan0.ms'
    os.system('rm -rf '+chan0_vis)
    os.system('rm -rf '+chan0_vis+'.flagversions')
    split(vis=final_out_file
          , datacolumn='DATA'
          , width=4000
          , outputvis=chan0_vis)

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# PRINT OUR TIME BENCHMARK
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

stop_time = time.time()
elapsed_time = (stop_time - start_time)/60.
print "This run took "+str(elapsed_time)+" minutes"
