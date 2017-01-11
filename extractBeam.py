# Snippet to extract the beam from a header and calculate a few
# related quantities.

# ... this is complicated by the possible presence of per-plane beams
if (header.keys()).count('restoringbeam') == 1:
    # ... the simple case
    beam = str(header['restoringbeam']['major']['value'])+'arcsec'
elif (header.keys()).count('perplanebeams') == 1:
    # ... per plane beams, pick the largest beam
    ppbdict = header['perplanebeams']['beams']
    beam = 0.0
    for plane in ppbdict.keys():
        this_plane = ppbdict[plane]
        for key in this_plane.keys():
            this_major = this_plane[key]["major"]["value"]
            if this_major > beam:
                beam = this_major
    beam = str(beam)+'arcsec'
else:
    print "extractBeam: could not find a beam."
    beam = None

pix_arcsec = abs(header['incr'][0]*180./np.pi*3600.)
beam_arcsec = float((beam.split('arcsec')[0]))
pix_per_beam = beam_arcsec/pix_arcsec
beam_area_pix = (beam_arcsec/pix_arcsec/2.)**2*np.pi/log(2)
