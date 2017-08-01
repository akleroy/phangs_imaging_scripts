execfile('../scripts/header_ngc4254north.py')
execfile('../scripts/line_list.py')

do_use_pbmask = False
linetag = 'co21'
specmode = 'cube'    
restfreq_ghz = line_list[linetag]
max_loop = 20
pb_limit = 0.75
uvtaper = None    

input_vis_7m = 'ngc4254north_7m_co21.ms'
cube_root_7m = 'ngc4254north_co21_7m'

input_vis_combo = 'ngc4254north_956_co21.ms'
cube_root_combo = 'ngc4254north_co21_12m+7m'

input_vis_12m = 'ngc4254north_12m_co21.ms'
cube_root_12m = 'ngc4254north_co21_12m'

do_image_7m = True
do_image_combo = False
do_image_12m = False

execfile('../scripts/phangsImagingPipeline.py')
