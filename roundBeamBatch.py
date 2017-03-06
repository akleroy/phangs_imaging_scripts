# This is temporary until we get the pure python one working.

gals = {'ngc0628':'arcsec' \
            , 'ngc1672':'arcsec' \
            , 'ngc3351':'arcsec' \
            , 'ngc3627':'arcsec' \
            , 'ngc4254':'arcsec' \
            , 'ngc4303':'arcsec' \
            , 'ngc4321':'arcsec' \
            , 'ngc4535':'arcsec' \
            , 'ngc5068':'arcsec' \
            , 'ngc6744north':'arcsec' \
            , 'ngc6744south':'arcsec' \
            }

for this_gal in gals.keys():
    cube_root='../release/v0p4/process/'+this_gal+'_co21'
    target_beam = gals[this_gal]
    execfile('roundBeamCASA.py')
