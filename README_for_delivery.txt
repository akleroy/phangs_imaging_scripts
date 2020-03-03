PHANGS Data Delivery Readme

This is a delivery of CO 2-1 data for PHANGS.

Files are named via:

galaxy_line_array_datatype

* The currenty delivery includes all PHANGS Cycle 2 and 3 targets and
  two archival targets imaged in a very similar way (NGC 1365 and NGC
  5128) also using CO 2-1. A single map is included for each galaxy.

* The current delivery includes only co21 data. Continuum and some
  C18O and CS54 can be constructed from the data, but tend to be low
  signal to noise. Processing these will come in future releases.

* The current observations all use all three parts of ALMA. The file
  names and headers refer to the arrays using the following code:

  12m : the ALMA main array of 12m dishes

  7m : the interferometric part of the ACA of 7m dishes

  tp : the single dish ("total power") part of the ACA

  12m+7m : the ALMA main array and 7m interferometric array imaged together

  12m+7m+tp : all interferometric data combined with the tp

  7m+tp : the 7m data combined with the tp

  In this release combination of interferometric and total power data
  is done via feathering. See the survey paper and the git repository
  of data reduction scripts for more details.

  All galaxies have been observed with the 7m and tp. As of this
  release a subset of galaxies have not been observed with the 12m
  data. As of the current release, one galaxy (the archival NGC 5128)
  has pathological problems with the single dish data.

  The header should reflect the combination of data used in the cube
  in the ARRAY keyword.

* All data are delivered in units of Kelvin, which is Rayleigh Jeans
  brightness temperature. The integrated intensity in the moment 0
  maps is in units of K*km/s and the moment 1 is in units of
  km/s. Units should be reflected in the BUNIT keyword.

* All cubes have convolved to have round beams and a single common
  beam (i.e., not frequency dependent). The intial images have
  frequency dependent beams and elliptical beams that reflect the u-v
  coverage of the data. The beam is reflected in the BMAJ and BMIN
  keywords.

* Two versions of the cubes are delivered: "pbcorr" and "flat." The
  pbcorr cubes account for the response of the primary beam. These are
  the correct intensity images of the sky, as best we know them. The
  noise varies across these images because there were less
  observations near the edge of the map.

  The flat versions are not corrected for the effect of the primary
  beam. These are useful for identification of signal because the
  noise is more nearly flat across the map..

* The release includes noise estimates position by position for each
  cube. These are appended with the _noise extension. These are
  constructed with a moving box and using some outlier rejection. The
  box on the sky is large enough that some care should be used near
  the edge of the mosaic, where the primary beam response changes
  quickly. There should be a few improvements in future relases.

* Cubes are provided at multiple resolution. With no extension, the
  resolution is the minimum round beam achievable from the basic
  imaging. Other resolutions are labeled in pc assuming the distance
  in the header.

* At each resolution, a mask is provided. This is basic bright signal
  masking, two channels at five sigma extended out to a lower
  significance contour. There are a number of important applications
  that would use different masking. Many of these can be achieved by
  applying a mask from a coarser resolution to the higher resolution
  data. This should be easy to do using the current release, and we
  will add these improved prodcuts in the upcoming releases.

* The mask is applied to the cube to produce a number of two
  dimensional data products:

  mom0, emom0 : the integrated line intensity within the mask

  mom1, emom1 : the intensity weighted mean velocity

  tpeak : the peak intensity inside the entire velocity range covered
  by the mask over each line of sight.

  tpeak_12p5kms : the peak intensity as above, but now using a version
  of the cube that has first been smoothed with a 5-channel tophat
  kernel across each spectrum.


