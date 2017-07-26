execfile('../scripts/header_ic5332.py')
execfile('../scripts/line_list.py')

do_use_pbmask = False
linetag = 'co21'
specmode = 'cube'    
restfreq_ghz = line_list[linetag]
max_loop = 20
pb_limit = 0.25
uvtaper = None    

input_vis_7m = 'ic5332_7m_co21.ms'
cube_root_7m = 'ic5332_co21_7m'

input_vis_combo = 'ic5332_925_co21.ms'
cube_root_combo = 'ic5332_co21_12m+7m'

input_vis_12m = 'ic5332_12m_co21.ms'
cube_root_12m = 'ic5332_co21_12m'

do_image_7m = True
do_image_combo = True
do_image_12m = True

execfile('../scripts/phangsImagingPipeline.py')
