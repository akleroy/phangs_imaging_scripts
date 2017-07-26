import numpy as np

infile = open('list_of_vis.txt', 'r')

dir_list = []
vis_list = []
while True:
    line  = infile.readline()    
    if len(line) == 0:
        break
    if line[0] == '#':
        continue
    words = line.split()    
    if len(words) < 2:
        continue
    dir_list.append(words[0])
    vis_list.append(words[1])

infile.close()

try:
    oversamp
except NameError:
    oversamp = 5

# Valid image sizes are even and multiples of 3, 5, 7
valid_sizes = []
for ii in range(10):
    for kk in range(3):
        for jj in range(3):
            valid_sizes.append(2**(ii+1)*5**(jj)*3**(kk))
valid_sizes.sort()
valid_sizes = np.array(valid_sizes)

outfile = open('image_params.txt', 'w')
outfile.write('# column 1: visibility file\n')
outfile.write('# column 2: cell size string\n')
outfile.write('# column 3: x pixel imsize\n')
outfile.write('# column 4: y pixel imsize\n')

for ii in range(len(vis_list)):
    
    this_dir = dir_list[ii]
    this_vis = vis_list[ii]
    input_vis = '../' + this_dir + '/' + this_vis

    # Cell size implied by baseline distribution
    au_cellsize, au_imsize, au_centralField = \
        au.pickCellSize(input_vis, imsize=True, npix=oversamp)
    xextent = au_cellsize*au_imsize[0]*1.2
    yextent = au_cellsize*au_imsize[1]*1.2

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

    # Now make the image size a nice round number
    need_cells_x = xextent / cell_size
    need_cells_y = yextent / cell_size

    cells_x = np.min(valid_sizes[valid_sizes > need_cells_x])
    cells_y = np.min(valid_sizes[valid_sizes > need_cells_y])

    image_size = [int(cells_x), int(cells_y)]
    cell_size_string = str(cell_size)+'arcsec'

    # Write to a text file
    line_out = ''
    line_out += this_vis+' '
    line_out += cell_size_string+' '
    line_out += str(image_size[0])+' '
    line_out += str(image_size[1])+' '
    line_out += '\n'

    outfile.write(line_out)
    print ""
    print line_out
    print ""

outfile.close()
