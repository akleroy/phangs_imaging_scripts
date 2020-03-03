"""
Utilities for defining file names.
"""

import os, ast

import logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

##############################################################
# Cube Names
##############################################################

def get_cube_filename(target=None, config=None, product=None,
                      ext=None, casa=False, casaext='.image'
                      ):
    """
    Get the file name for a data cube using the pipeline convention.

    {target}_{config}_{product}_{ext}.{fits or casaext}

    """

    if target is None:
        logging.error("Need a target.")            
        return(None)
    
    if config is None:
        logging.error("Need a configuration.")            
        return(None)
    
    if product is None:
        logging.error("Need a product.")            
        return(None)
    
    if type(target) is not type(''):
        logging.error("Target needs to be a string.", target)
        return(None)
    
    if type(product) is not type(''):
        logging.error("Product needs to be a string.", product)
        return(None)
    
    if type(config) is not type(''):
        logging.error("Config needs to be a string.", config)
        return(None)
    
    filename = target+'_'+config+'_'+product
    if ext is not None:
        if type(ext) is not type(''):
            logging.error("Ext needs to be a string or None.", ext)
            return(None)
        filename += '_'+ext

    if not casa:
        filename += '.fits'
    else:
        filename += casaext

    return(filename)

##############################################################
# Visibility Names After Processing
##############################################################

def get_vis_filename(target=None, config=None, product=None, 
                     ext=None, suffix=None,
                     ):
    """
    Get the file name for a target, config, product visibility
    combination. Optionally include an extension and a suffix. 
    Convention is:
    
    {target}_{config}_{product}_{ext}.ms{.suffix}
    
    """

    if target is None:
        logging.error("Need a target.")            
        return(None)
    
    if config is None:
        logging.error("Need a configuration.")            
        return(None)
    
    if product is None:
        logging.error("Need a product.")            
        return(None)
    
    if type(target) is not type(''):
        logging.error("Target needs to be a string.", target)
        return(None)

    if type(product) is not type(''):
        logging.error("Product needs to be a string.", product)
        return(None)

    if type(config) is not type(''):
        logging.error("Config needs to be a string.", config)
        return(None)
        
    filename = target+'_'+config+'_'+product
    if ext is not None:
        if type(ext) is not type(''):
            logging.error("Ext needs to be a string or None.", ext)
            return(None)
        if len(ext) > 0:
            filename += '_'+ext
        
    filename += '.ms'

    if suffix is not None:
        if type(suffix) is not type(''):
            logging.error("Suffix needs to be a string or None.", suffix)
            return(None)
        filename += '.'+suffix        

    return(filename)

##############################################################
# Visibility Names During Staging
##############################################################

def get_staged_msname(target=None, project=None, array_tag=None,
                      obsnum=None, ext=None, suffix=None,
                      ):
    """
    Get the file name for a staged measurement set (where we have
    not yet extracted the spectral product). Optionally include an
    extension and a suffix.  Convention is:
    
    {target}_{project}_{array_tag}_{ext}.ms{.suffix}
    
    """

    if target is None:
        logging.error("Need a target.")            
        return(None)
    
    if project is None:
        logging.error("Need a project code.")            
        return(None)
    
    if array_tag is None:
        logging.error("Need an array_tag.")            
        return(None)
        
    if type(target) is not type(''):
        logging.error("Target needs to be a string."+str(target))
        return(None)

    if type(project) is not type(''):
        logging.error("Project needs to be a string."+str(project))
        return(None)

    if type(array_tag) is not type(''):
        logging.error("Array tag needs to be a string."+str(config))
        
    if type(obsnum) is not type(''):
        logging.error("Obsnum needs to be a string."+str(obsnum))

    filename = target+'_'+project+'_'+array_tag+'_'+obsnum
    if ext is not None:
        if type(ext) is not type(''):
            logging.error("Ext needs to be a string or None.", ext)
            return(None)
        if len(ext) > 0:
            filename += '_'+ext
        
    filename += '.ms'

    if suffix is not None:
        if type(suffix) is not type(''):
            logging.error("Suffix needs to be a string or None.", suffix)
            return(None)
        filename += '.'+suffix        

    return(filename)
