execfile('../scripts/header_ngc1672.py')
execfile('../scripts/line_list.py')

do_use_pbmask = False
linetag = 'co21'
specmode = 'cube'
restfreq_ghz = line_list[linetag]
max_loop = 20
pb_limit = 0.75
uvtaper = None    

input_vis_7m = 'ngc1672_7m_co21.ms'
cube_root_7m = 'ngc1672_co21_7m'

input_vis_combo = 'ngc1672_956_co21.ms'
cube_root_combo = 'ngc1672_co21_12m+7m'

input_vis_12m = 'ngc1672_12m_co21.ms'
cube_root_12m = 'ngc1672_co21_12m'

smallscalebias_7m = 0.6
do_image_7m = True
smallscalebias_combo = 0.8
do_image_combo = False
smallscalebias_12m = 0.8
do_image_12m = False

execfile('../scripts/phangsImagingPipeline.py')
