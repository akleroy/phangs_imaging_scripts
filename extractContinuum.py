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
    print "I will default to flagging only the CO, 13CO, and C18O lines."
    lines_to_flag = lines_co+lines_13co+lines_c18o

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
        print "... flagging lines "+str(lines_to_flag)

        vm = au.ValueMapping(this_infile)

        spw_flagging_string = ''
        first = True
        for line in lines_to_flag:
            rest_linefreq_ghz = line_list[line]

            shifted_linefreq_hz = rest_linefreq_ghz*(1.-source_vel_kms/sol_kms)*1e9
            hi_linefreq_hz = rest_linefreq_ghz*(1.-(source_vel_kms-vwidth_kms/2.0)/sol_kms)*1e9
            lo_linefreq_hz = rest_linefreq_ghz*(1.-(source_vel_kms+vwidth_kms/2.0)/sol_kms)*1e9

            spw_list = au.getScienceSpwsForFrequency(this_infile,
                                                     shifted_linefreq_hz)
            if spw_list == []:
                continue

            print "Found overlap for "+line
            for this_spw in spw_list:
                freq_ra = vm.spwInfo[this_spw]['chanFreqs']
                chan_ra = np.arange(len(freq_ra))
                to_flag = (freq_ra >= lo_linefreq_hz)*(freq_ra <= hi_linefreq_hz)
                to_flag[np.argmin(np.abs(freq_ra - shifted_linefreq_hz))]
                low_chan = np.min(chan_ra[to_flag])
                hi_chan = np.max(chan_ra[to_flag])                
                this_spw_string = str(this_spw)+':'+str(low_chan)+'~'+str(hi_chan)
                if first:
                    spw_flagging_string += this_spw_string
                    first = False
                else:
                    spw_flagging_string += ','+this_spw_string

        print "... proposed flagging "+spw_flagging_string

        if spw_flagging_string != '':
            flagdata(vis=this_infile,
                     spw=spw_flagging_string,
                     )

        if do_statwt:
            print "... deriving emprical weights using STATWT"
            statwt(vis=this_infile,
                   datacolumn='DATA')
            
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# EXTRACT THE LINE OF INTEREST FROM EACH FILE
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
  
if do_average:

    print "--------------------------------------------------------"
    print "extractContinuum: (3) Average and collapse"
    print "--------------------------------------------------------"

    files_to_concat = []
    for this_tag in calibrated_files.keys():
        this_infile = out_root+'_'+this_tag+'_cont_copy.ms'
        this_outfile = out_root+'_'+this_tag+'_cont.ms'

        print "... averaging data for "+this_infile

        os.system('rm -rf '+this_outfile)
        os.system('rm -rf '+this_outfile+'.flagversions')
        split(vis=this_infile,
              outputvis=this_outfile,
              width=10000,
              datacolumn='DATA',
              keepflags=False)
        files_to_concat.append(this_outfile)
    
    final_out_file = out_root+'_'+tag+'_cont.ms'
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

    # .........................................
    # Clean up
    # .........................................

    for this_tag in calibrated_files.keys():
        this_infile = out_root+'_'+this_tag+'_cont_copy.ms'
        os.system('rm -rf '+this_infile)
        os.system('rm -rf '+this_infile+'.flagversions')

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# PRINT OUR TIME BENCHMARK
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

extract_stop_time = time.time()
extract_elapsed_time = (extract_stop_time - extract_start_time)/60.

print "****************************************"
print "extractContinuum took "+"{:6.2f}".format(extract_elapsed_time)+" minutes"
print "****************************************"
