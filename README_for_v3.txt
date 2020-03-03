PHANGS-ALMA Data Delivery README

This is a delivery of CO 2-1 data for PHANGS-ALMA.

------------------------
IMPORTANT PREFACE
------------------------

This release, release 3.4, is expected to be the last "v3"
release. Then v4 is expected to represent a candidate delivery to the
observatory and public release. Compared to previous v3 versions, the
most important change is that v3.4 fixes a bug introduced sometime
after v0p6 in the stitching of mosaics. If you are using a multi-part
target you should use these data moving forward if at all possible.

The release also includes a additional galaxies and some other
improvements to the data and processing.

------------
VERSION: 3.4
------------

This is the first delivery of our full PHANGS-ALMA data set. It is the
third major delivery to the team, and this is the first version of
this release. 

We label this version 3.4. Subsequent versions of this release using
the same basic processing will be 3.5, 3.6, 3.7, etc. Deliveries that
implement significant processing improvements will be released using
higher base version numbers (e.g. version 4.X).

This is NOT a frozen release. Versions 4.0, 5.0, etc. are intended to
be frozen releases. These will stay stable over many months (or
years). Versions labeled 3.X (3.1, 3.2, etc.) are incremental
updates. These will change as we implement improvements.

This data delivery includes several pilot programs, the PHANGS-ALMA
Large Program, several follow-up programs, and several archival data
sets. If you use these data, please include the acknowledgment,
project codes, and references below.

As of this release, almost all data have been delivered. The handful
of outstanding follow up observations will be incorporated into future
version releases, but team members should not anticipate another dozen
galaxies.

Starting with version 3.4, the data products in a PHANGS-ALMA data
release are associated with a specific version of the PHANGS sample
table. For this release, the associated sample table is:

phangs_sample_table_v1p2.fits

That file is copied into the release. Please note that the format of
the sample table will change with the v4 ALMA release. Then we expect
the format to remain relatively stable after this.

---------------------
ACKNOWLEDGMENTS
---------------------

If you use these data products, please cite the PHANGS-ALMA survey
paper:

A. K. Leroy, E. Schinnerer et al. (in preparation). 

*Temporarily*, while the survey paper remains in preparation, if you
use the "pilot" (2015.1.00925 and 2015.1.00956) sample before
publication of the survey paper please add a citation that "These maps
have been previous presented in J. Sun et al. (2018ApJ...860..172S)."
If you specifically use the NGC 0628 data
(ADS/JAO.ALMA#2012.1.00650.S) before the publication of the survey
paper, please reference that these data have also been presented in
Leroy et al. 2016 (2016ApJ...831...16L) and Kreckel et al. 2018
(2018ApJ...863L..21K)."

Please also add the following acknowledgments:

This paper makes use of the following ALMA data:
ADS/JAO.ALMA#2012.1.00650.S, ADS/JAO.ALMA#2013.1.00803.S,
ADS/JAO.ALMA#2013.1.01161.S, ADS/JAO.ALMA#2015.1.00782.S,
ADS/JAO.ALMA#2015.1.00925.S, ADS/JAO.ALMA#2015.1.00956.S,
ADS/JAO.ALMA#2015.1.00121.S, ADS/JAO.ALMA#2016.1.00386.S,
ADS/JAO.ALMA#2017.1.00392.S, ADS/JAO.ALMA#2017.1.00886.L,
ADS/JAO.ALMA#2017.1.00766.S, ADS/JAO.ALMA#2018.1.01321.S,
ADS/JAO.ALMA#2018.1.01651.S, ADS/JAO.ALMA#2018.A.00062.S

ALMA is a partnership of ESO (representing its member states), NSF
(USA) and NINS (Japan), together with NRC (Canada), MOST and ASIAA
(Taiwan), and KASI (Republic of Korea), in cooperation with the
Republic of Chile. The Joint ALMA Observatory is operated by ESO,
AUI/NRAO and NAOJ. The National Radio Astronomy Observatory is a
facility of the National Science Foundation operated under cooperative
agreement by Associated Universities, Inc.

---------------------
HELP AND FEEDBACK
---------------------

We welcome feedback and suggestions on data products and
documentation. Please direct suggestions to the PHANGS ALMA data
reduction group (adr@phangs.groups.io). For the immediate next few
releases (while we work out remaining bugs), it may be useful to
direct detailed questions directly to Adam Leroy
(leroy.42@osu.edu). Both is probably best.

---------------------
PROCESSING
---------------------

These data have all been calibrated in the pipeline for which they
were delivered. For most, this is a version of CASA 5.X.X but earlier
data sets include some data calibrated in 4.X.X. After calibration,
data were staged for imaging, regridded, and combined using CASA
5.4.0. Imaging proceeded in CASA 5.4.0 for this release.

After imaging, data were exported as FITS files and post
processed. They were convolved to have a round synthesized beam,
converted to brightness temperature (Kelvin) units. Interferometric
and total power data were "feathered" together. Multi-part mosaics
were combined via inverse-noise-squared weighting using the primary
beam response to indicate the sensitivity. Data were then down-sampled
to have the minimum pixel size needed to critically sample the beam
and remove empty space.

To create the products in this release, each cube is convolved to a
succession of physical resolutions. The assumed distance to the target
galaxy is recorded in the FITS header for all convolved data products
using the DIST keyword (see below). At each resolution, the pipeline
estimates the three dimensional noise distribution. "Signal" masks are
created at each resolution that identify regions with detectable
signal. By combining masks created at low resolution and high
sensitivity, high completeness masks are also constructed. These are
used to build "strict" and "broad" data products at different physical
resolutions.

---------------------
WHAT IS IN HERE
---------------------

This delivery includes three directories. The differences are
explained below, then the individual map types are explained.

--------------
strict_maps/
--------------

This directory contains "moment" maps (for a loose definition of
moment) at a succession of physical resolutions. These are "strict" in
the sense that the maps only reflect signal identified in the cube at
the target resolution. 

ADVANTAGES: This approach has two main advantages: (1) The masking
suppresses noise, so any calculation that is unstable in the presence
of noise (moment 2 calculation is the classic example) may want to use
these data. (2) The selection function applied to the cube to build
them is very straightforward. These are almost the same as the maps
used in Sun et al. 2018, Utomo et al. 2018, and Gallagher et al. 2018.

DISADVANTAGES: The downside of this approach is that some (but not
all) of these maps suffer from substantial incompleteness at higher
resolution. The next release will include the number in the header,
but for now see Sun et al. 2018 for examples of the degree of
incompleteness. Here "incompleteness" means that a significant amount
of the flux in the cube (even flux that was cleaned) does not pass the
signal to noise cuts used to build the maps.

Use the strict maps if you want low noise, can live with some
incompleteness, and want a very reproducible calculation.

--------------
broad_maps/
--------------

This directory contains "moment" maps (for a loose definition of
moment) at a succession of physical resolutions. These are "broad" in
the sense that the maps are created using regions of the cube known to
have signal at coarsers resolutions. There is no requirement tha the
signal be detected at the current resolution.

ADVANTAGES: This approach has the major advantage of high completeness
and high covering fraction.

DISADVANTAGES: This approach lowers the signal-to-noise and includes
noise, including potentially noise-only lines of sight, in the final
data products. Calculations that become unstable in the presence of
noise do not work as well in these maps.

Use the broad maps if you want high completeness and can live with
some additional noise.

--------------
cubes/
--------------

This directory contains data cubes. For the moment, we provide the
cubes at their native (round beam) resolution, primary beam corrected,
downsampled, and converted to Kelvin.

N.B.: The set of released data cubes is still TBD. The data volume for
the full set of created products is large (pushing 10TB). Our plan,
not implemented yet, is to add a "support" or "intermediate" directory
that contains a very large data volume.

---------------------
Provided maps
---------------------

Currently we provide the following maps. The file name gives the line
(e.g., co21 indicates 12CO J=2->1 data), the array combination, the
resolution, and the data product. When no resolution is specified, the
map contains native resolution data.

RESOLUTION: The approximate physical resolution is given in the file
name. The exact angular resolution is given by the BMAJ and BMIN
keywords (almost always equal) in the header. We also include the
adopted distance as a keyword in the header. Note that we allow a
tolerance of +/- 10% in convolution to target physical
resolution. That is, a map can have native resolution 85pc or 75pc and
still have the filename labeled 80pc.

INCLUDED ARRAYS: The combination of arrays is indicated by the file
name, e.g., 12m+7m indicates joint 12m array and 7m array imaging, but
NOT total power. 12m+7m+tp indicates data including information from
all arrays.

.....................
Intensity metrics
.....................

mom0 - integrated intensity
emom0 - uncertainty in mom0

EXPLANATION: This is direct integration of the cube along the velocity
axis inside the relevant mask. The uncertainty comes from error
propagation assuming independent velocity channels and using the
empirical noise estimates. The units are K*km/s.

The mom0 is our best tracer of how much molecular gas there is along
the line of sight.

tpeak - peak temperature 
tpeak12p5kms - the peak temperature over any 12.5km/s window

EXPLANATION: Peak intensity (in Kelvin) along the velocity axis for
each line of sight. For tpeak12p5kms, the peak along velocity axis
after smoothing the cube with a 5 channel boxcar along the spectral
axis (but not downsampling). Note that for the tpeak and tpeak12p5kms
calculations, the mask used covers all spatial pixels in all channels
where the signal mask has any coverage. This ensures relatively
uniform noise properties across the map.

The tpeak12p5kms is a good way to see the detailed structure of the
galaxy at high signal to noise.

.....................
Velocity metrics
.....................

mom1 - intensity weighted mean velocity
emom1 - uncertainty in mom1

EXPLANATION: This is the intensity-weighted mean velocity calculated
inside the mask. In the case of the "broad" moment1, we reject pixels
that deviate strongly from the low resolution velocity field as likely
outliers. This has the advantage of producing smooth velocity fields,
but could in principle miss some very weird kinematics. For sightlines
with a single velocity component, the moment1 reflects the average
velocity of the emission. In the case of two or more components, the
moment1 will sit intermediate between the two.

The uncertainty in the moment1 reflects the statistical uncertainty in
the cube, propagated into the map assuming independent velocity
channels.

.....................
Line profile metrics
.....................

ew - equivalent width
eew - uncertainty in ew, not currently populated

EXPLANATION: This is the "equivalent width" line width
diagnostic. This is the line integral (moment 0) divided by the peak
intensity. That is, this is the rectangular width needed to supply the
full line width at the peak intensity. A prefactor is then applied to
recast as the equivalent sigma for the case of a Gaussian line
profile. This is a highly robust statistic in the sense that it
behaves well in the presence of noise or multiple components. It can
easily miss subtleties in the line profile and has some dependence on
spectral resolution.

mom2 - rms velocity dispersion
emom2 - uncertainty in mom2, not currently populated

EXPLANATION: This is the "moment2" line width diagnostic. This is the
intensity-weighted second moment, measuring the rms scatter about the
intensity weighted mean velocity. This metric is highly sensitive to
the inclusion of noise. In the presence of multiple spectral peaks, it
will also be sensitive to the spread between peaks. In that sense,
convergence between this and the ew is a crude diagnostic of
Gaussianity. Because of the sensitivity to noise, it is recommended to
only use moment2 from the 'strict' maps.

Note that the line width maps provided are not yet corrected for the
line spread function (channel width + channel-to-channel
correlation). Nor are they corrected for biases due to finite
sensitivity (this is a particular issue for moment 2).
