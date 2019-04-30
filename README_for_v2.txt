PHANGS Data Delivery Readme

This is a delivery of CO 2-1 data for PHANGS.

-------
VERSION
-------

This is the first delivery of version 2, release 2.1.

This includes the pilots, the PHANGS-ALMA Large Program, and several
archival data sets. As of 2.1 not all Large Program data are delivered
and total power data also lag behind other data. In future releases
only data including all arrays will be delivered, but for now we
deliver maps for 12m+7m (i.e., no total power) AND 12m+7m+tp.

This deliver focuses only on data including the 12m array, the next
version will also include the 7m-only maps. Those data are already
available in a preliminary format on the drive.

For reference, these data have all been reduced in the pipeline for
which they were delivered. Then the visibility data are staged in CASA
5.3.0 and the maps have been imaged in 5.1.1 with square maps.

Future versions will override this release, it is NOT intended from my
side that we keep a large set of releases live on the server. Along
the same lines, I expect that these represent better science products
than the v0p6 release.

---------------------
WHAT IS IN HERE
---------------------

These delivery have three directories, two of them populated.

--------------
strict_maps/
--------------

This directory contains moment maps at a succession of
resolutions. These are "strict" in the sense that they have been made
by only signal masking the cube at this resolution. The advantage of
this approach is that the selection function applied to the cube to
build them is very straightforward. These are almost the same as the
maps used in Sun et al. 2018 and then Utomo et al. 2018 and Gallagher
et al. 2018. Use them if you want a very reproducible calculation.

The downside of this approach is that the maps suffer from substantial
incompleteness at higher resolution. The next release will include the
number in the header, but for now see Sun et al. 2018 for examples of
the degree of incompleteness. Here "incompleteness" means that a
significant amount of the flux in the cube (even flux that was
cleaned) does not pass the signal to noise cuts used to build the
maps.

For this release the algorithm is - mask starts at 3 channels above
3.5 sigma and expands in all dimensions out to 2 channels at 2
sigma. Then for moment 0 and equivalent width calculation, the mask is
expanded two channels in each direction along the velocity axis
only. For moment 1 and moment 2 calculation the mask is used as is
with no expansion.

For the tpeak and tpeak12p5kms calculations, the mask used covers all
spatial pixels in all channels where the signal mask has any coverage.

Provided products are:

mom0 - integrated intensity
emom0 - uncertainty in mom0

mom1 - intensity weighted mean velocity
emom1 - uncertainty in mom1

ew - equivalent width (line width measurement)
eew - uncertainty in ew, not currently populated

mom2 - rms velocity dispersion (not yet corrected for sensitivity or resolution)
emom2 - uncertainty in mom2, not currently populated
p2e - peak to edge ratio for mom2 (can be used to correct for sensitivity)

tpeak - peak temperature 
tpeak12p5kms - the peak temperature after smoothing the cube with a 5-channel boxcar along the velocity direction

For the maps, we provide a map at many resolutions. The resolution is
given in the file name. It can be calculated from the BMAJ keyword
along with the distance keyword in the header.

The combination of arrays is indicated by the file name, e.g., 12m+7m or 12m+7m+tp

--------------
cubes/
--------------

This contains the data cubes for the 12m+7m and 12m+7m+tp data that we
currently have in hand.

We provide round-beam cubes, primary beam corrected in Kelvin
units. The array combination is indicated in the file name. For now,
we push only the native resolution cubes to save space on the server.

--------------
hybrid_maps/
--------------

This directory is not currently populated. It will hold the moment
maps constructed using cross-scale techniques. These, e.g., include
most of the flux but at lower signal to noise or have better spatial
coverage compared to the strict maps. But the selection function is
more complex.

These are being made but not yet quite ready to be distributed.
