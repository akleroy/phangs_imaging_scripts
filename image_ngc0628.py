execfile('../scripts/header_ngc0628.py')
execfile('../scripts/line_list.py')

do_use_pbmask = False
linetag = 'co21'
specmode = 'cube'    
restfreq_ghz = line_list[linetag]
max_loop = 20
pb_limit = 0.75
uvtaper = None    

input_vis_7m = 'ngc0628_7m_co21.ms'
cube_root_7m = 'ngc0628_co21_7m'

input_vis_combo = 'ngc0628_650_co21.ms'
cube_root_combo = 'ngc0628_co21_12m+7m'

input_vis_12m = 'ngc0628_12m_co21.ms'
cube_root_12m = 'ngc0628_co21_12m'

do_image_7m = True
do_image_combo = False
do_image_12m = False

execfile('../scripts/phangsImagingPipeline.py')
    
