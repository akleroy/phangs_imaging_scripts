# This script extracts a line measurement set from calibrated
# visibility data sets. The idea is to yield a single, simple
# visibility set that can be fed to our imaging scripts.

# USAGE: Run *once* with "do_copy=True" and "do_split=True" for each
# data set. Then run with "do_extract=True" and "do_combine=True" once
# *for each line.*
#

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# CHECKS AND DEFAULTS
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

tested_versions = ['4.6.0','4.7.0','4.7.1']
this_version = (casa['build']['version']).split('-')[0]
if this_version not in tested_versions:
    print "The script hasn't been verified for this version of CASA."
    print "This version of CASA is "+this_version
    print "Tested versions are "+str(tested_versions)

sol_kms = 2.99e5

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# TIMER
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

import time
extract_start_time = time.time()

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# CONTROL FLOW
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

try:
    do_copy
except NameError:    
    do_copy = True

try:
    do_split
except NameError:    
    do_split = True

try:
    do_extract
except NameError:    
    do_extract = True

try:
    do_combine
except NameError:    
    do_combine = True

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# INPUTS
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

print "****************************************"
print "********** extractLineData.py **********"
print "****************************************"

abort = False

try:
    calibrated_files
except NameError:
    print "Please define a calibrated_files file dictionary."
    print "Format is calibrated_files[tag]:ms_to_copy."
    print "(I will turn off all steps for this round)."
    abort = True

try:
    field_names
except NameError:
    print "Considering all fields."
    field_names = ''

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
    abort = True

try:
    tag
except NameError:
    print "Please define a tag for the output file names. Use tag="
    print "Our current practice is to set this to a line identifier." 
    print "(I will turn off all steps for this round)."
    abort = True

if abort:
    do_copy = False
    do_split = False
    do_extract = False
    do_combine = False

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

try:
    use_contsub
except NameError:
    print "I will NOT look for a continuum subtracted version."
    use_contsub = False

contsub_string = ''
if use_contsub:
    print "I will use the continuum subtracted version of the data."
    contsub_string = '.contsub'

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# COPY DATA FROM ITS ORIGINAL LOCATION
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

if do_copy:

    print "--------------------------------------------------------"
    print "extractLineData: (1) Copying data from original location"
    print "--------------------------------------------------------"

    for this_tag in calibrated_files.keys():
        in_file = calibrated_files[this_tag]
        out_file = out_root+'_'+this_tag+'_copied.ms'
        os.system('rm -rf '+out_file)
        os.system('rm -rf '+out_file+'.flagversions')
        command = 'cp -r '+in_file+' '+out_file
        print command
        os.system(command)

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# TIME AVERAGE AND SPLIT OUT ONLY THE SOURCE
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

if do_split:

    print "--------------------------------------------------------"
    print "extractLineData: (2) Splitting out source"
    print "--------------------------------------------------------"
      
    print "... time smoothing and extracting source data for all files."

    for this_tag in calibrated_files.keys():
        this_infile = out_root+'_'+this_tag+'_copied.ms'
        this_outfile = out_root+'_'+this_tag+'.ms'
        
        os.system('rm -rf '+this_outfile)
        os.system('rm -rf '+this_outfile+'.flagversions')

        print "... splitting file "+this_infile

        split(vis=this_infile
              , intent ='OBSERVE_TARGET#ON_SOURCE'
              , field = field_names
              , datacolumn='DATA'
              , timebin=timebin
              , outputvis=this_outfile)        

        if do_statwt:
            print "... deriving emprical weights using STATWT"
            statwt(vis=this_outfile,
                   datacolumn='DATA')
            
    # ......................................
    # Delete intermediate files
    # ......................................

    print "... cleaning up copied files."

    for this_tag in calibrated_files.keys():
        copied_file = out_root+'_'+this_tag+'_copied.ms'
        os.system('rm -rf '+copied_file)
        os.system('rm -rf '+copied_file+'.flagversions')

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# EXTRACT THE LINE OF INTEREST FROM EACH FILE
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
  
if do_extract:

    print "--------------------------------------------------------"
    print "extractLineData: (3) Extract the line of interest"
    print "--------------------------------------------------------"

    for this_tag in calibrated_files.keys():

        this_infile = out_root+'_'+this_tag+'.ms'+contsub_string
        this_outfile = out_root+'_'+this_tag+'_'+linetag+'.ms'

        print "... extracting line data for "+this_infile

        # .........................................
        # Calculate and recast as strings
        # .........................................

        # Crude, but probably fine for extragalactic work. Figure out what
        # SPWs contain the line in question.

        target_freq_ghz = restfreq_ghz*(1.-source_vel_kms/sol_kms)
        
        spw_list = au.getScienceSpwsForFrequency(this_infile,
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
        chan_width_hz = au.getChanWidths(this_infile, spw_list_string)
        target_width_hz = chan_dv_kms/sol_kms*restfreq_ghz*1e9
        rebin_factor = min(target_width_hz / chan_width_hz)
    
        if rebin_factor < 2:
            chanbin = 1
        else:
            chanbin = int(floor(rebin_factor/2.))

        # .........................................
        # Concert to strings and summarize
        # .........................................

        print "FILE:", this_infile
        print "LINE TAG: ", linetag
        print "REST FREQUENCY: ", restfreq_string
        print "SPECTRAL WINDOWS: ", spw_list
        print "TARGET CHANNEL WIDTH: ", chan_dv_string
        print "SUPPLIED SOURCE VELOCITY: ", str(source_vel_kms)
        print "DESIRED VELOCITY WIDTH: ", str(vwidth_kms)
        print "START VELOCITY: ", start_vel_string
        print "NUMBER OF CHANNELS: ", str(nchan)
        print "CHANNELS TO BIN TOGETHER FIRST: ", chanbin

        # .........................................
        # Regrid
        # .........................................

        print "... carrying out channel averaging and hanning smooth first."

        os.system('rm -rf '+this_outfile+'.temp')
        os.system('rm -rf '+this_outfile+'.temp.flagversions')

        if chanbin > 1:
            chanaverage = True
        else:
            chanaverage = False

        mstransform(vis=this_infile,
                    outputvis=this_outfile+'.temp',
                    spw=spw_list_string,
                    datacolumn='DATA',
                    chanaverage=chanaverage,
                    chanbin=chanbin,
                    hanning=False,
                    # Does this matter? It shouldn't, but AS thinks it still might
                    interpolation='cubic',
                    )

        print "... now regridding to the desired velocity grid."

        os.system('rm -rf '+this_outfile)
        os.system('rm -rf '+this_outfile+'.flagversions')

        mstransform(vis=this_outfile+'.temp',
                    outputvis=this_outfile,
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

        os.system('rm -rf '+this_outfile+'.temp')
        os.system('rm -rf '+this_outfile+'.temp.flagversions')

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# CONCATENATE, STATWT, AND PRODUCE CHANNEL 0
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

if do_combine:

    print "--------------------------------------------------------"
    print "extractLineData: (4) Combining all data for the line"
    print "--------------------------------------------------------"

    files_to_concat = []
    for this_tag in calibrated_files.keys():
        this_infile = out_root+'_'+this_tag+'_'+linetag+'.ms'
        files_to_concat.append(this_infile)

    final_out_file = out_root+'_'+tag+'_'+linetag+'.ms'
    os.system('rm -rf '+final_out_file)
    os.system('rm -rf '+final_out_file+'.flagversions')
    concat(vis=files_to_concat,
           concatvis=final_out_file)

    # .........................................
    # A final statwt
    # .........................................

    if do_statwt:
        print "... deriving empirical weights using STATWT."
        statwt(vis=final_out_file,
               datacolumn='DATA')

    # ......................................
    # Make "Channel 0" CO measurement sets
    # ......................................

    print '... making a "channel zero" data set using SPLIT.'

    chan0_vis = out_root+'_'+tag+'_'+linetag+'_chan0.ms'
    os.system('rm -rf '+chan0_vis)
    os.system('rm -rf '+chan0_vis+'.flagversions')
    split(vis=final_out_file
          , datacolumn='DATA'
          , spw=''
          , outputvis=chan0_vis
          , width=10000)

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# PRINT OUR TIME BENCHMARK
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

extract_stop_time = time.time()
extract_elapsed_time = (extract_stop_time - extract_start_time)/60.

print "****************************************"
print "extractLineData took "+"{:6.2f}".format(extract_elapsed_time)+" minutes"
print "****************************************"
