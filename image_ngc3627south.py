execfile('../scripts/line_list.py')
execfile('../scripts/header_ngc3627south.py')

do_use_pbmask = False
linetag = 'co21'
specmode = 'cube'    
restfreq_ghz = line_list[linetag]
max_loop = 20
pb_limit = 0.75
uvtaper = None    

input_vis_7m = 'ngc3627south_7m_co21.ms'
cube_root_7m = 'ngc3627south_co21_7m'

input_vis_combo = 'ngc3627south_956_co21.ms'
cube_root_combo = 'ngc3627south_co21_12m+7m'

input_vis_12m = 'ngc3627south_12m_co21.ms'
cube_root_12m = 'ngc3627south_co21_12m'

do_image_7m = True
do_image_combo = False
do_image_12m = False

smallscalebias_7m = 0.6
snr_thresh = 4.0

execfile('../scripts/phangsImagingPipeline.py')
