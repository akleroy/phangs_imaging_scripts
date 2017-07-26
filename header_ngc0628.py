out_root = 'ngc0628'
tag = '650'
phase_center = 'J2000 01h36m41.7s +15d47m01'
source_vel_kms = 657
vwidth_kms = 200

dir_12m = '../../Cycle1_650/2012.1.00650.S/science_goal.uid___A002_X5a9a13_X579/group.uid___A002_X5a9a13_X57a/member.uid___A002_X5a9a13_X57b/calibrated/'

dir_7m = '../../Cycle1_650/2012.1.00650.S/science_goal.uid___A002_X5a9a13_X579/group.uid___A002_X5a9a13_X57a/member.uid___A002_X5a9a13_X57d/calibrated/'

calibrated_files = {'12m_1':dir_12m+'uid___A002_X5b2f01_X3f.ms.split.cal.M74/',
                    '12m_2':dir_12m+'uid___A002_X8081ba_X4018.ms.split.cal.M74/',
                    '12m_3':dir_12m+'uid___A002_X8081ba_X4527.ms.split.cal.M74/',
                    '12m_4':dir_12m+'uid___A002_X80c782_X1911.ms.split.cal.M74/',
                    '12m_5':dir_12m+'uid___A002_X8204db_X611.ms.split.cal.M74/',
                    '12m_6':dir_12m+'uid___A002_X95b353_X471.ms.split.cal.M74/',
                    '12m_7':dir_12m+'uid___A002_X960614_X2d5b.ms.split.cal.M74/',
                    '12m_8':dir_12m+'uid___A002_X966cea_X96c.ms.split.cal.M74/',
                    '7m_1':dir_7m+'uid___A002_X6f1341_X80b.ms.split.cal',
                    '7m_2':dir_7m+'uid___A002_X6f2c6e_Xb0a.ms.split.cal',
                    '7m_3':dir_7m+'uid___A002_X7fc9da_X1f0d.ms.split.cal',
                    '7m_4':dir_7m+'uid___A002_X7fc9da_X4b45.ms.split.cal',
                    '7m_5':dir_7m+'uid___A002_X8081ba_X11fb.ms.split.cal',
                    '7m_6':dir_7m+'uid___A002_X8081ba_X3e04.ms.split.cal',
                    '7m_7':dir_7m+'uid___A002_X8081ba_X44b8.ms.split.cal',
                    '7m_8':dir_7m+'uid___A002_X8081ba_X85e.ms.split.cal',
                    '7m_9':dir_7m+'uid___A002_X8081ba_Xce1.ms.split.cal',
                    '7m_10':dir_7m+'uid___A002_X8204db_X4f.ms.split.cal',
                    }

clean_mask_file = '../clean_masks/ngc0628_co21_clean_mask.fits'
