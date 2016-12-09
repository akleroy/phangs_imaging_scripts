import numpy as np

# Modification of analysis utilities version to deal with some lack of
# files in our processed visibilities.

def estimateSynthesizedBeam(vis, field='', useL80method=True, verbose=False,
                            frequencyGhz=None):
    """
    Modified from Todd Hunter's to allow the representative frequency
    to be set by hand because the natural path appears to lose the
    summary XML file that the script wants.
    """
    if useL80method:
        result = au.getBaselineStats(vis, field=field, percentile=80, verbose=verbose)
    else:
        result = au.getBaselineStats(vis, field=field, verbose=verbose)
    if (result == None): return

    if frequencyGhz == None:
        vm = au.ValueMapping(vis)
        nspw = 0
        sumFreq = 0.0
        for spw in vm.spwInfo.keys():
            nspw += 1
            sumFreq += vm.spwInfo[spw]['meanFreq']/1e9        
        frequencyGHz = sumFreq/(nspw*1.0)
    print "Representative frequency for beam estimate: ", frequencyGHz

    if useL80method:
        L80 = result[0] # meters
        arcsec = 3600*np.degrees(1)*0.574*au.c_mks/(frequencyGHz*1e9*L80)
    else:
        maxBaseline = result[2] # meters
        arcsec = au.printBaselineAngularScale(maxBaseline, frequencyGHz)

    return arcsec
