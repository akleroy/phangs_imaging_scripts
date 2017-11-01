# At some point they stopped shipping "scriptForImagingPrep.py" - this
# simple script does the aggolmeration. For the Large Program this is
# mostly vestigal, we will just label each delivered MS as its own
# thing in the MS key. But for the earlier cases with many executions,
# this script gets copied into the scripts/ directory of the
# individual reduction and used to build "calibrated_final.ms"

import glob
in_dir = '../calibrated/'
vislist=glob.glob(in_dir+'*.ms.split.cal')

concatvis=in_dir+'calibrated.ms'
os.system('rm -rf '+concatvis)
os.system('rm -rf '+concatvis+'.flagversions')
concat(vis=vislist,
       concatvis=concatvis)

sourcevis=in_dir+'calibrated_source.ms'
rmtables(sourcevis)
os.system('rm -rf ' + sourcevis + '.flagversions')
split(vis=concatvis,
      intent='*TARGET*',
      outputvis=sourcevis,
      datacolumn='data')

os.system('mv -i ' + sourcevis + ' ' + in_dir+'calibrated_final.ms')
