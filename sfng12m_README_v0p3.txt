
---

This is the 12m-only SFNG survey of nearby galaxies. The survey
targeted the CO 2-1 line, the C18O 2-1 line, and the 1mm continuum. We
produced images at approximately the native resolution of the data for
a Brigg's robust parameter of 0.5, convolving asymmetric beams to the
nearest reasonable round beam (e.g., 1.30x0.95 goes to 1.4x1.4). We
also constructed a tapered version of each image using a 2.5x2.5
arcsecond uv taper and then smoothing to a final resolution of 3x3
arcseconds.

This is release to the team v0p2. Feedback is welcome - to the whole
list for discussion or to leroy.42@osu.edu and schruba@mpe.mpg.de for
directed feedback.

Included data products are:

- CO21 and C18O cubes, "Chan 0" for continuum, CO21, and C18O

- _cube or _taper : data cube, primary beam corrected, in K

- _taper_mask : a mask at the tapered resolution, generally this
  includes some faint emission or empty space in the high resolution
  cube but does a good job of getting all the flux.

_ _taper_mask_hires : the tapered mask but on the high resolution
  astrometry. Useful for using the low resolution emission to mask the
  high resolution data.

- _bright_mask : a mask containing the bright emission at the native
  resolution. It will miss some flux compared to the tapered mask but
  contains only secure detections, mostly. Used for making the moment
  1.

- mom0 and emom0 : the first moment and statistical uncertainty on the
  first moment in K km/s.

- mom1 and emom1 : the first moment (velocity field) and statistical
  uncertainty on it in km/s,

- tpeak and tpeak_12p5kms : the peak intensity along the line of sight
  in K at the native 2.5km/s resolution and the peak intensity along
  the line of sight after smoothing the cube with a 5-channel boxcar
  kernel (i.e., 12.5km/s wide).

- pb : the coverage of the mosaic. This is the response of the
  telescope to intensity on the sky and should be inversely
  proportional to the noise in the primary beam corrected map.

The next releases are anticipated to:

- adjust velocity ranges slightly to encompass the galaxy with a more
  or less uniform amount of line free channels per side.

- deal with an issue where NGC 3627 is currently cleaned without a
  clean mask because CASA "test clean" in 4.6.0 diverges for that
  galaxy with a mask.

- provide feathered data products for those galaxies with IRAM 30-m
  HERA data (NGC 3351, 3627, 4254, 4303, 4321, 4535).

- include the delivered 7m data.

Future developments are anticipated to include checking the scripts in
to a git repository and, eventually, imaging following the imaging
subteam's (Andreas, Annie, Jerome, Kaz) approved path.

Suggestions welcome.

History

v0p1:

- first release, CO 2-1 only

- round beam cubes, tapered cubes, moment maps, masks

v0p2:

- added NGC 0628 processed in same way

- added C18O and continuum + channel 0

- fixed a mislabeling so that "round" is now NOT pb corrected (i.e., it is flat noise)
