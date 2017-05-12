# Script to feather for PHANGS

gal_list = [
    'ic5332',
    'ngc0628',
    'ngc1087',
    'ngc1300',
    'ngc1385',
    'ngc1433',
    'ngc1512',
    'ngc1566',
    'ngc1672',
    'ngc2835',
    'ngc3351',
    'ngc3627',
    'ngc4254',
    'ngc4303',
    'ngc4321',
    'ngc4535',
    'ngc5068',
    'ngc6744north',
    'ngc6744south'
    ]

for gal in gal_list:
    print "Importing "+gal

    importfits(fitsimage=gal+'_co21_pbcorr_round.fits',
               imagename=gal+'_co21_pbcorr_round.image',
               zeroblanks=True, overwrite=True)
    importfits(fitsimage=gal+'_tp.fits',
               imagename=gal+'_tp.image',
               zeroblanks=True,overwrite=True)

for gal in gal_list:
    print "Feathering "+gal
    
    feather(imagename=gal+'_co21_feathered.image',
            highres=gal+'_co21_pbcorr_round.image',
            lowres=gal+'_tp.image')

for gal in gal_list:
    print "Exporting "+gal

    exportfits(imagename=gal+'_co21_feathered.image',
               fitsimage=gal+'_co21_feathered.fits',
               velocity=True, dropdeg=True, dropstokes=True, 
               overwrite=True, bitpix=16)
