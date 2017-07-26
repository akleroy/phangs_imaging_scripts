execfile('../scripts/header_ngc2835.py')
execfile('../scripts/line_list.py')

do_use_pbmask = False
linetag = 'co21'
specmode = 'cube'    
restfreq_ghz = line_list[linetag]
max_loop = 20
pb_limit = 0.25
uvtaper = None    

input_vis_7m = 'ngc2835_7m_co21.ms'
cube_root_7m = 'ngc2835_co21_7m'

input_vis_combo = 'ngc2835_925_co21.ms'
cube_root_combo = 'ngc2835_co21_12m+7m'

input_vis_12m = 'ngc2835_12m_co21.ms'
cube_root_12m = 'ngc2835_co21_12m'

do_image_7m = True
do_image_combo = True
do_image_12m = True

execfile('../scripts/phangsImagingPipeline.py')
