try:
    input_vis
except NameError:
    print "Please define an input visibility via input_vis=XXX."
    abort = True

# Factor by which to oversample the beam
try:
    oversamp
except NameError:
    oversamp = 5

# Cell size implied by baseline distribution
au_cellsize, au_imsize, au_centralField = \
    au.pickCellSize(input_vis, imsize=True, npix=oversamp)
xextent = au_cellsize*au_imsize[0]*1.2
yextent = au_cellsize*au_imsize[1]*1.2

# If a u-v taper is present, then take that as an alternate beam and
# add it in quadrature.
try:
    uvtaper
except NameError:
    uvtaper = None

if uvtaper != None:
    print "I account for a uv taper of "+str(uvtaper)+" arcsec"
    cell_size_taper = uvtaper / oversamp
    au_cellsize = sqrt(au_cellsize**2 + cell_size_taper**2)

# Make the cell size a nice round number
if au_cellsize < 0.1:
    cell_size = au_cellsize
if au_cellsize >= 0.1 and au_cellsize < 0.5:
    cell_size = floor(au_cellsize/0.05)*0.05
if au_cellsize >= 0.5 and au_cellsize < 1.0:
    cell_size = floor(au_cellsize/0.1)*0.1
if au_cellsize >= 1.0 and au_cellsize < 2.0:
    cell_size = floor(au_cellsize/0.25)*0.25
if au_cellsize >= 2.0 and au_cellsize < 5.0:
    cell_size = floor(au_cellsize/0.5)*0.5
if au_cellsize >= 5.0:
    cell_size = floor(au_cellsize/1.0)*0.5

# If we taper, then the cell size may be forced
print "I calculated cell size: ", cell_size
try:
    force_cell_size
except NameError:
    force_cell_size = None

if force_cell_size != None:
    cell_size = force_cell_size
    print "I will force cell size: ", cell_size

    # Now make the image size a nice round number
    need_cells_x = xextent / cell_size
    need_cells_y = yextent / cell_size

    # Valid image sizes are even and multiples of 3, 5, 7
    valid_sizes = []
    for ii in range(10):
        for kk in range(3):
            for jj in range(3):
                valid_sizes.append(2**(ii+1)*5**(jj)*3**(kk))
    valid_sizes.sort()
    valid_sizes = np.array(valid_sizes)

    cells_x = np.min(valid_sizes[valid_sizes > need_cells_x])
    cells_y = np.min(valid_sizes[valid_sizes > need_cells_y])

    image_size = [int(cells_x), int(cells_y)]
    cell_size_string = str(cell_size)+'arcsec'

    print "I propose the following:"
    print "... cell size = "+cell_size_string
    print "... image size = ", image_size
