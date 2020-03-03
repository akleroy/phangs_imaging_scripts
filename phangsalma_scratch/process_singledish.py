# run in working_dir = '../singledish/data/'

os.system('rm -rf ALMA_TP.M74.v0p2.image')
importfits(fitsimage='ALMA_TP.M74.v0p2.blanked.fits',
           imagename='ALMA_TP.M74.v0p2.image')

imreframe(imagename='ALMA_TP.M74.v0p2.image',
          restfreq='230.53800GHz')

os.system('rm -rf ALMA_TP.M74.v0p2_rebin.image')
imrebin(imagename='ALMA_TP.M74.v0p2.image',
        outfile='ALMA_TP.M74.v0p2_rebin.image',
        factor=[1,1,7,1])
        
exportfits(imagename='ALMA_TP.M74.v0p2_rebin.image',
           fitsimage='ALMA_TP.M74.v0p2_rebin.fits',
           velocity=True, dropdeg=True, overwrite=True)

# import/export to strip degenerate axes

os.system('rm -rf ALMA_TP.NGC_1087.CO21.v0p2.gildas.image')
importfits(fitsimage='ALMA_TP.NGC_1087.CO21.v0p2.gildas.blanked.fits',
           imagename='ALMA_TP.NGC_1087.CO21.v0p2.gildas.image')

exportfits(imagename='ALMA_TP.NGC_1087.CO21.v0p2.gildas.image',
           fitsimage='ALMA_TP.NGC_1087.CO21.v0p2.casa.fits',
           velocity=True, dropdeg=True, overwrite=True)

# import/export to strip degenerate axes, delete "BLANK" keyword

importfits(fitsimage='ALMA_TP.NGC_1300.CO21.v0p2.image.VLSRK.gildas.blanked.fits',
           imagename='ALMA_TP.NGC_1300.CO21.v0p2.image.VLSRK.gildas.image' 
           , overwrite=True)

exportfits(imagename='ALMA_TP.NGC_1300.CO21.v0p2.image.VLSRK.gildas.image',
           fitsimage='ALMA_TP.NGC_1300.CO21.v0p2.image.VLSRK.casa.fits',
           velocity=True, dropdeg=True, overwrite=True)

# import/export to strip degenerate axes, delete "BLANK" keyword

importfits(fitsimage='ALMA_TP.NGC_1385.CO21.v0p2.image.VLSRK.gildas.blanked.fits',
           imagename='ALMA_TP.NGC_1385.CO21.v0p2.image.VLSRK.gildas.image' 
           , overwrite=True)

exportfits(imagename='ALMA_TP.NGC_1385.CO21.v0p2.image.VLSRK.gildas.image',
           fitsimage='ALMA_TP.NGC_1385.CO21.v0p2.image.VLSRK.casa.fits',
           velocity=True, dropdeg=True, overwrite=True)

# import/export to strip degenerate axes, delete "BLANK" keyword

importfits(fitsimage='ALMA_TP.NGC_1433.CO21.v0p2.gildas.NuclearFix.blanked.fits',
           imagename='ALMA_TP.NGC_1433.CO21.v0p2.gildas.NuclearFix.image' 
           , overwrite=True)

exportfits(imagename='ALMA_TP.NGC_1433.CO21.v0p2.gildas.NuclearFix.image',
           fitsimage='ALMA_TP.NGC_1433.CO21.v0p2.casa.fits',
           velocity=True, dropdeg=True, overwrite=True)

# import/export to strip degenerate axes, delete "BLANK" keyword

importfits(fitsimage='ALMA_TP.NGC_1566.CO21.v0p2.gildas.blanked.fits',
           imagename='ALMA_TP.NGC_1566.CO21.v0p2.gildas.image',
           overwrite=True)

exportfits(imagename='ALMA_TP.NGC_1566.CO21.v0p2.gildas.image',
           fitsimage='ALMA_TP.NGC_1566.CO21.v0p2.casa.fits',
           velocity=True, dropdeg=True, overwrite=True)

# redo export

importfits(fitsimage='ALMA_TP.NGC1365.CO21.v0p2.image.fits',
           imagename='ALMA_TP.NGC1365.CO21.v0p2.image',
           overwrite=True)

exportfits(imagename='ALMA_TP.NGC1365.CO21.v0p2.image',
           fitsimage='ALMA_TP.NGC1365.CO21.v0p2.image.casa.fits',
           velocity=True, dropdeg=True, overwrite=True)
