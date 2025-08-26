# INSIDE CASA
import sys
sys.path.insert(1, os.path.abspath('../'))
from phangsPipeline import utilsFieldSelection

"""
    2015.1.00956.S  uid://A001/X2fb/X27d    12m  DataSet_01
    2015.1.00956.S  uid://A001/X2fb/X27f    7m   DataSet_14
    2015.1.00956.S  uid://A001/X2fb/X283    tp   DataSet_20
    2023.1.01182.S  uid___A001_X3621_X13ef  12m  Dataset_10
"""

utilsFieldSelection.process_ms_list(
    ms_list = [
    '/data1/dzliu/Work_PHANGS/2017_Large_Program/2015.1.00956.S/Level_2_Calib/DataSet_01/calibrated/uid___A002_Xaf4574_X17bc.ms.split.cal', 
    '/data1/dzliu/Work_PHANGS/2017_Large_Program/2015.1.00956.S/Level_2_Calib/DataSet_14/calibrated/uid___A002_Xb62a5b_X8184.ms.split.cal', 
    '/data1/dzliu/Work_PHANGS/2017_Large_Program/2015.1.00956.S/Level_2_Calib/DataSet_14/calibrated/uid___A002_Xb64387_X61ae.ms.split.cal', 
    '/data1/dzliu/Work_PHANGS/2024_Highres_Centers/2023.1.01182.S/Level_2_Calib/DataSet_10/calibrated/uid___A002_X110a1a2_X2ad5b_targets.ms', 
    '/data1/dzliu/Work_PHANGS/2024_Highres_Centers/2023.1.01182.S/Level_2_Calib/DataSet_10/calibrated/uid___A002_X110fd13_X9858_targets.ms', 
    ], 
    lower_left_ra_dec = '10:43:59.3138 +11:41:52.303', 
    upper_right_ra_dec = '10:43:56.6774 +11:42:34.171', 
    verbose = True, 
)

