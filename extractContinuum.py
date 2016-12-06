# This script extracts a data set ready for continuum imaging from
# calibrated visibility data sets.
#
# USAGE:
#
# TO DO:
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
execfile('../scripts/line_list.py')
import numpy as np

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# TIMER
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

import time
extract_start_time = time.time()

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# CONTROL FLOW
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

try:
    do_recopy
except NameError:
    do_recopy = True

try:
    do_flag
except NameError:    
    do_flag = True

try:
    do_average
except NameError:    
    do_average = True

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# INPUTS
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

print "*****************************************"
print "********** extractContinuum.py **********"
print "*****************************************"

abort = False

try:
    calibrated_files
except NameError:
    print "Please define a calibrated_files file dictionary."
    print "Format is calibrated_files[tag]:ms_to_copy."
    print "(I will turn off all steps for this round)."
    abort = True

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

# .............................................................
# Trigger abort
# .............................................................
    
if abort:
    do_recopy = False
    do_flag = False
    do_average = False

# .............................................................
# Define the lines to flag
# .............................................................

try:
    lines_to_flag
except NameError:
    print "I will default to flagging only the CO and 13CO lines."
    lines_to_flag = lines_co.append(lines_13co)

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# MAKE A CONTINUUM COPY OF THE DATA
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

if do_recopy:

    print "--------------------------------------------------------"
    print "extractContinuum: (1) Making a continuum copy"
    print "--------------------------------------------------------"

    for this_tag in calibrated_files.keys():
        in_file = out_root+'_'+this_tag+'.ms'
        out_file = out_root+'_'+this_tag+'_cont_copy.ms'
        os.system('rm -rf '+out_file)
        os.system('rm -rf '+out_file+'.flagversions')
        command = 'cp -r '+in_file+' '+out_file
        print command
        os.system(command)

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# FLAG LINES
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

if do_flag:

    print "--------------------------------------------------------"
    print "extractContinuum: (2) Flagging spectral lines"
    print "--------------------------------------------------------"

    for this_tag in calibrated_files.keys():
        this_infile = out_root+'_'+this_tag+'_cont_copy.ms'
        
        print "... flagging lines in file "+this_infile

        vm = au.ValueMapping(this_infile)

        spw_flagging_string = ''
        for line in lines_to_flag:
            rest_linefreq_ghz = line_list[line]

            shifted_linefreq_hz = rest_linefreq_ghz*(1.-source_vel_kms/sol_kms)
            hi_linefreq_hz = rest_linefreq_ghz*(1.-(source_vel_kms-vwidth_kms)/sol_kms)
            lo_linefreq_hz = rest_linefreq_ghz*(1.-(source_vel_kms+vwidth_kms)/sol_kms)

            spw_list = au.getScienceSpwsForFrequency(this_infile,
                                                     shifted_linefreq_ghz*1e9)
            if spw_list == []:
                continue
            
            for this_spw in spw_list:
                freq_ra = vm.spwInfo[this_spw]['chanFreqs']
                chan_ra = np.arange(len(freq_ra))
                to_flag = (freq_ra >= low_linefreq_hz)*(freq_ra <= hi_linefreq_hz)
                to_flag[np.argmin(np.abs(freq_ra - shifted_linefreq_hz))]
                low_chan = np.min(chan_ra[to_flag])
                hi_chan = np.min(chan_ra[to_flag])                
                this_spw_string = str(this_spw)+':'+str(low_chan)+'~'+str(hi_chan)
                spw_flagging_string += this_spw_string+','

        print "... proposed flagging "+spw_flagging_string

        flagdata(vis=this_infile,
                 spw=spw_string,
                 )

        if do_statwt:
            print "... deriving emprical weights using STATWT"
            statwt(vis=this_infile,
                   datacolumn='DATA')
            
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# EXTRACT THE LINE OF INTEREST FROM EACH FILE
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
  
if do_extract:

    print "--------------------------------------------------------"
    print "extractContinuum: (3) Average and collapse"
    print "--------------------------------------------------------"

    for this_tag in calibrated_files.keys():

        this_infile = out_root+'_'+this_tag+'.ms'
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
                    hanning=True,
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
          , outputvis=chan0_vis)

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# PRINT OUR TIME BENCHMARK
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

extract_stop_time = time.time()
extract_elapsed_time = (extract_stop_time - extract_start_time)/60.

print "****************************************"
print "extractLineData took "+"{:6.2f}".format(extract_elapsed_time)+" minutes"
print "****************************************"
