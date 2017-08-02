import os

only = ['ngc5128']
skip = []

infile = open('../scripts/extraction_key.txt', 'r')

dir_list = []
out_root_list = []

while True:
    line  = infile.readline()    
    if len(line) == 0:
        break
    if line[0] == '#':
        continue
    words = line.split()
    if len(words) < 8:
        continue

    this_out_root = words[0]
    
    add_this = True
    if len(only) > 0:
        add_this = False
        for this_root in only:
            if this_out_root == this_root:
                add_this = True

    if len(skip) > 0:
        for this_root in skip:
            if this_out_root == this_root:
                add_this = False

    if add_this:
        out_root_list.append(this_out_root)
        if this_out_root.find('north') != -1:
            dir_list.append(this_out_root[0:this_out_root.find('north')])
        elif this_out_root.find('south') != -1:
            dir_list.append(this_out_root[0:this_out_root.find('south')])
        else:
            dir_list.append(this_out_root)

infile.close()

for ii in range(len(dir_list)):
    os.chdir('../'+dir_list[ii])
    out_root = out_root_list[ii]
    execfile('../scripts/extract_data.py')

os.chdir('../scripts/')
