"""Utilities for angular and physical resolutions.
"""

import re
import numpy as np

regex_psep = re.compile(r'([0-9eE.+-]+)(p)([0-9eE.+-]+)(.*)', re.IGNORECASE)
regex_nounit = re.compile(r'^([0-9eE.+-]+)$', re.IGNORECASE)
regex_arcsec = re.compile(r'^([0-9eE.+-]+)[ \t]*(arcsec|\")$', re.IGNORECASE)
regex_arcmin = re.compile(r'^([0-9eE.+-]+)[ \t]*(arcmin|\')$', re.IGNORECASE)
regex_degree = re.compile(r'^([0-9eE.+-]+)[ \t]*(degree|deg)$', re.IGNORECASE)
regex_pc = re.compile(r'^([0-9eE.+-]+)[ \t]*(parsec|pc)$', re.IGNORECASE)
regex_kpc = re.compile(r'^([0-9eE.+-]+)[ \t]*(kpc)$', re.IGNORECASE)
regex_Mpc = re.compile(r'^([0-9eE.+-]+)[ \t]*(Mpc)$', re.IGNORECASE)

def is_angular_resolution(res, return_value=False):
    """Check if a string is an angular resolution. 
    
    If return_value then return the value in arcsec as well. 
    
    If a string without any unit is given, we assume it is angular resolution.
    """
    out_flag = False
    out_value = None
    if isinstance(res, str):
        res = res.strip()
        if regex_psep.match(res):
            res = regex_psep.sub(r'\1.\3\4', res) # if the input is like '5p00', we convert it to '5.00'
        for regex_obj, mult_factor in list(zip([regex_nounit, regex_arcsec, regex_arcmin, regex_degree], [1.0, 1.0, 60.0, 3600.0])):
            if regex_obj.match(res):
                out_flag = True
                out_value = float(regex_obj.search(res).group(1)) * mult_factor
    else:
        if isinstance(res, float) or isinstance(res, int):
            out_flag = True
            out_value = float(res)
    # 
    if return_value:
        return out_flag, out_value
    else:
        return out_flag


def is_physical_resolution(res, return_value=False):
    """Check if a string is a physical resolution. 
    
    If return_value then return the value in parsec as well. 
    """
    out_flag = False
    out_value = None
    if isinstance(res, str):
        res = res.strip()
        if regex_psep.match(res):
            res = regex_psep.sub(r'\1.\3\4', res) # if the input is like '5p00', we convert it to '5.00'
        for regex_obj, mult_factor in list(zip([regex_pc, regex_kpc, regex_Mpc], [1.0, 1e3, 1e6])):
            if regex_obj.match(res):
                out_flag = True
                out_value = float(regex_obj.search(res).group(1)) * mult_factor
    else:
        if isinstance(res, float) or isinstance(res, int):
            out_flag = True
            out_value = float(res)
    # 
    if return_value:
        return out_flag, out_value
    else:
        return out_flag


def is_distance(distance, return_value=False):
    """Check if a string is a distance ending with parsec units. 
    
    If return_value then return the value in Mega parsec. 
    """
    out_flag = False
    out_value = None
    if isinstance(distance, str):
        for regex_obj, mult_factor in list(zip([regex_pc, regex_kpc, regex_Mpc], [1e-6, 1e-3, 1.0])):
            if regex_obj.match(distance):
                out_flag = True
                out_value = float(regex_obj.search(distance).group(1)) * mult_factor # return value in units of Mpc
    else:
        if isinstance(distance, float) or isinstance(distance, int):
            out_flag = True
            out_value = float(distance)
    # 
    if return_value:
        return out_flag, out_value
    else:
        return out_flag
    

def get_tag_for_angular_resolution(res):
    """Input an angular resolution string or value, output a formatted string tag to be used in filenames.
    """
    check_flag, res_value = is_angular_resolution(res, return_value=True)
    if check_flag:
        sigfigs = 2
        fmt_str = '{:.%df}'%(sigfigs)
        res_str = fmt_str.format(np.round(res_value, decimals=sigfigs)).replace('.','p')
        return res_str
    else:
        raise Exception('The input resolution string "'+str(res)+'" is not an angular resolution!')
        return None


def get_tag_for_physical_resolution(res):
    """Input a physical resolution string or value, output a formatted string tag to be used in filenames.
    """
    check_flag, res_value = is_physical_resolution(res, return_value=True)
    if check_flag:
        sigfigs = 0
        fmt_str = '{:.%df}'%(sigfigs)
        res_str = fmt_str.format(np.round(res_value, decimals=sigfigs)).replace('.','p') + "pc"
        return res_str
    else:
        raise Exception('The input resolution string "'+str(res)+'" is not a physical resolution!')
        return None


def get_tag_for_res(res):
    """Return a tag string to be used in filenames given a resolution string like either '5.0arcsec' 
         or '80pc', or a resolution value in arcsec. 
       
       The returned tag string is always formatted like '5p00' for an angular resolution, or 
         like '80pc' for a physical resolution.
    """
    if is_angular_resolution(res):
        return get_tag_for_angular_resolution(res)
    elif is_physical_resolution(res):
        return get_tag_for_physical_resolution(res)
    else:
        raise Exception('The input resolution string "'+str(res)+'" seems neither an angular nor a physical resolution!')
        return None


def get_angular_resolution_from_physical_resolution(res, distance):
    """Return the angular resolution in arcsec, given a physical resolution and a distance.
    
    The input physical resolution can either be a string or a value in parsec. 
    
    The input distance can also either be a string or a value in Mpc. 
    """
    # 
    res_check_flag, res_value_in_pc = is_physical_resolution(res, return_value=True)
    # 
    distance_check_flag, distance_value_in_Mpc = is_distance(distance, return_value=True)
    # 
    if res_check_flag and distance_check_flag:
        kpc2arcsec = 1e-3/distance_value_in_Mpc/np.pi*180.0*3600.0
        return res_value_in_pc / 1e3 * kpc2arcsec
        # res_value_in_arcsec = res_value_in_pc / 1e3 * kpc2arcsec
        # res_value_in_pc = res_value_in_arcsec * (distance_value_in_Mpc*1e6) * np.pi / 180.0 / 3600.0
    else:
        raise Exception('The input resolution "'+str(res)+'" and distance "'+str(distance)+'" seem incorrect. Please input a physical resolution in parsec and a distance in parsec!')
        return None


def get_angular_resolution_for_res(res, distance=None):
    """Return an angular resolution value in units of arcsec given a resolution string like either '5.0arcsec' 
         or '80pc', or a resolution value in arcsec. 
    """
    res_check_flag, res_value_in_arcsec = is_angular_resolution(res, return_value=True)
    if res_check_flag:
        #print('res "'+str(res)+'" is an angular resolution')
        return res_value_in_arcsec
    else:
        #print('res "'+str(res)+'" is a physical resolution')
        if distance is None:
            raise Exception('Need a distance for calculating angular resolution from the given physical resolution of '+str(res))
        return get_angular_resolution_from_physical_resolution(res, distance)



#def get_angular_resolution_for_config(config):
#    """Return the rough angular resolution range for a given config like '12m+7m', '12m', '7m'.
#    """
#    WIP









