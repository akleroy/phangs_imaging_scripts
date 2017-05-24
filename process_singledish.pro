tol = 1d-6

cube = readfits('ALMA_TP.M74.v0p2.fits', hdr)
blank = where(abs(cube - sxpar(hdr,'BLANK')) lt tol)
cube[blank] = !values.f_nan
writefits, 'ALMA_TP.M74.v0p2.blanked.fits', cube, hdr

cube = readfits('ALMA_TP.NGC_1087.CO21.v0p2.gildas.fits', hdr)
blank = where(abs(cube - sxpar(hdr,'BLANK')) lt tol)
cube[blank] = !values.f_nan
writefits, 'ALMA_TP.NGC_1087.CO21.v0p2.gildas.blanked.fits', cube, hdr

cube = readfits('ALMA_TP.NGC_1300.CO21.v0p2.image.VLSRK.gildas.fits', hdr)
blank = where(abs(cube - sxpar(hdr,'BLANK')) lt tol)
cube[blank] = !values.f_nan
writefits, 'ALMA_TP.NGC_1300.CO21.v0p2.image.VLSRK.gildas.blanked.fits', cube, hdr

cube = readfits('ALMA_TP.NGC_1385.CO21.v0p2.image.VLSRK.gildas.fits', hdr)
blank = where(abs(cube - sxpar(hdr,'BLANK')) lt tol)
cube[blank] = !values.f_nan
writefits, 'ALMA_TP.NGC_1385.CO21.v0p2.image.VLSRK.gildas.blanked.fits', cube, hdr

cube = readfits('ALMA_TP.NGC_1433.CO21.v0p2.gildas.NuclearFix.fits', hdr)
blank = where(abs(cube - sxpar(hdr,'BLANK')) lt tol)
cube[blank] = !values.f_nan
writefits, 'ALMA_TP.NGC_1433.CO21.v0p2.gildas.NuclearFix.blanked.fits', cube, hdr

cube = readfits('ALMA_TP.NGC_1566.CO21.v0p2.gildas.fits', hdr)
blank = where(abs(cube - sxpar(hdr,'BLANK')) lt tol)
cube[blank] = !values.f_nan
writefits, 'ALMA_TP.NGC_1566.CO21.v0p2.gildas.blanked.fits', cube, hdr

