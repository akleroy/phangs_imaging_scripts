# This script subtracts the continuum from a data set, attempting to
# exclude line channels from the fit.
#
# USAGE:
#
# TO DO:
#

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# CHECKS AND DEFAULTS
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

tested_versions = ['5.0.0']
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
    do_contsub
except NameError:
    do_contsub = True

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# INPUTS
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

print "*****************************************"
print "********** subtractContinuum.py *********"
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
    do_contsub

# .............................................................
# Define the lines to flag
# .............................................................

try:
    lines_to_flag
except NameError:
    print "I will default to flagging only the CO, 13CO, and C18O lines."
    lines_to_flag = lines_co+lines_13co+lines_c18o


# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# FLAG LINES
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

if do_contsub:

    print "--------------------------------------------------------"
    print "extractContinuum: (2) Flagging spectral lines"
    print "--------------------------------------------------------"

    for this_tag in calibrated_files.keys():
        this_infile = out_root+'_'+this_tag+'.ms'
        
        print "... continuum subtracting "+this_infile
        print "... avoiding lines "+str(lines_to_flag)

        vm = au.ValueMapping(this_infile)

        spw_flagging_string = ''
        first = True
        for spw in vm.spwInfo.keys():
            this_spw_string = str(spw)+':0'
            if first:
                spw_flagging_string += this_spw_string
                first = False
            else:
                spw_flagging_string += ','+this_spw_string            

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

        print "... proposed channels to avoid "+spw_flagging_string

        os.system('rm -rf '+this_infile+'.contsub')
        uvcontsub(vis=this_infile,
                  fitspw=spw_flagging_string,
                  excludechans=True)

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# PRINT OUR TIME BENCHMARK
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

extract_stop_time = time.time()
extract_elapsed_time = (extract_stop_time - extract_start_time)/60.

print "****************************************"
print "extractContinuum took "+"{:6.2f}".format(extract_elapsed_time)+" minutes"
print "****************************************"
