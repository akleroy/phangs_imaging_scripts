execfile('../scripts/header_ngc3627north.py')
execfile('../scripts/line_list.py')

do_use_pbmask = False
linetag = 'co21'
specmode = 'cube'    
restfreq_ghz = line_list[linetag]
max_loop = 20
pb_limit = 0.75
uvtaper = None    

input_vis_7m = 'ngc3627north_7m_co21.ms'
cube_root_7m = 'ngc3627north_co21_7m'

input_vis_combo = 'ngc3627north_956_co21.ms'
cube_root_combo = 'ngc3627north_co21_12m+7m'

input_vis_12m = 'ngc3627north_12m_co21.ms'
cube_root_12m = 'ngc3627north_co21_12m'

# Custom call for the 7m here due to divergence.
smallscalebias_7m = 0.6
snr_thresh = 4.0

do_image_7m = True
do_image_combo = False
do_image_12m = False
#do_multiscale_clean = True

execfile('../scripts/phangsImagingPipeline.py')

#do_multiscale_clean = True
#do_image_7m = False

#snr_thresh = 0.0
#smallscalebias_combo = 0.6
#do_image_combo = True
#smallscalebias_12m = 0.6
#do_image_12m = True

#execfile('../scripts/phangsImagingPipeline.py')
