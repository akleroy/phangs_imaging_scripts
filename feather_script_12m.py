# Script to feather the 12m data for PHANGS

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
    'ngc3627north',
    'ngc3627south',
    'ngc4254north',
    'ngc4254south',
    'ngc4303',
    'ngc4321north',
    'ngc4321south',
    'ngc4535',
    'ngc5068north',
    'ngc5068south',
    'ngc6744north',
    'ngc6744south'
    ]

for gal in gal_list:
    print "Importing "+gal
    importfits(fitsimage=gal+'_co21_12m+7m_flat_round.fits',
               imagename=gal+'_co21_12m+7m_flat_round.image',
               zeroblanks=True, overwrite=True)
    importfits(fitsimage=gal+'_tp_tapered_12m+7m.fits',
               imagename=gal+'_tp_tapered_12m+7m.image',
               zeroblanks=True,overwrite=True)

for gal in gal_list:
    print "Feathering "+gal    
    feather(imagename=gal+'_co21_12m+7m_feathered.image',
            highres=gal+'_co21_12m+7mflat_round.image',
            lowres=gal+'_tp_tapered_12m+7m.image')

for gal in gal_list:
    print "Exporting "+gal
    exportfits(imagename=gal+'_co21_12m+7m_feathered.image',
               fitsimage=gal+'_co21_12m+7m_feathered.fits',
               velocity=True, dropdeg=True, dropstokes=True, 
               overwrite=True, bitpix=16)
