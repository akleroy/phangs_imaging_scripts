"""
This is a dummy CleanCall class for dry run only, or to be inheritted by casaImagingRoutines.CleanCall.
"""

import numpy as np
import re

import logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

#region class CleanCall

class CleanCall:
    
    def __init__(self, infile_list=[]):

        ############################################
        # Parameters used by the CleanCall object. #
        ############################################

        self.infile_list = infile_list
        self.logfile = None

        ###################################
        # Initialize the clean parameters #
        ###################################
        
        self.clean_params = {}
        self.reset_params()

    #########################
    # Attributes From Files #
    #########################

    # Read parameters from one file

    def reset_params(self):
        """
        Reset the clean call paameters to a dictionary containing only
        the parameters in the supplied files, read in the order that
        they are supplied by the user.
        """

        self.clean_params = {}
        for infile in self.infile_list:
            self.read_one_file(infile)

    def read_one_file(self, fname=None):
        """
        Read the clean call attributes from an input file of the kind
        generated by save_inputs.
        """

        logger.info("Reading a clean input parameter file.")
        logger.info("Reading: "+fname)

        infile = open(fname, 'r')
        
        lines_read = 0
        while True:
            line = infile.readline()
            line = line.replace(" ","")
            line = line.replace("\n","")
            if len(line) == 0:
                break
            if line[0] == '#' or line == '\n':
                continue
            if line.count('=') == 0:
                continue
            exec(re.sub(r'^([a-zA-Z0-9_]+) *= *(.+) *$', r'self.clean_params["\1"] = \2', line) )
            continue

        return()
        
    #####################
    # Attribute setting #
    #####################

    # Include specific routines only where there is logic to be
    # implemented, else just use "set."

    def set_param(self, key, value, nowarning=False):
        """
        Set a clean parameter to the specified value. Set nowarning to
        suppress a check that the field already exists.
        """
        if key in self.clean_params.keys() or nowarning:
            self.clean_params[key] = value
            return()
        else:
            raise Exception('The clean_params do not currently include key  "'+str(key)+'"')
        return()

    def set_restfreq_ghz(self, value=None):
        """
        Set the rest frequency used for cube imaging. Expects a value in GHz.
        """
        self.clean_params['restfreq']=str(value)+'GHz'
        return()
    
    def set_reffreq_ghz(self, value=None):
        """
        Set the refrence frequency used in continuum imaging. Expects
        a value in GHz.
        """
        self.clean_params['reffreq']=str(value)+'GHz'
        return()
    
    def set_multiscale_arcsec(self, scales_in_arcsec=[]):
        """
        Set the scales for deconvolution in acseconds. Requires that a
        cell size already be defined so that these can be translated
        into pixel units.
        """
        cell_in_pix = self.get_cell_in_arcsec()
        if cell_size is None:
            return()

        scales_in_pix = []
        for this_scale_arcsec in scales_in_arcsec:
            this_scale_pix = this_scale_arcsec/cell_in_pix
            scales_in_pix.append()

        scales_in_pix.sort()
        self.clean_params['scales'] = scales_in_pix

    def set_round_uvtaper_arcsec(self, taper=0.0):
        """
        Sets a round outer UV taper of the value provided in units of arcseconds.
        """
        self.clean_params['uvtaper'] = [str(taper)+'arcsec',str(taper)+'arcsec','0deg']
        return()
    
    ####################
    # Attribute access #
    ####################
    
    def get_cell_in_arcsec(self):
        """
        Return the cell size (assumed in arcsec string format) as a float in arcsec.
        """
        if len(self.clean_params['cell']) == 0:
            return(None)
        if self.clean_params['cell'] is None:
            return(None)
        if type(self.clean_params['cell']) == type([]):
            cell_string = self.clean_params['cell'][0]
        else:
            cell_string = self.clean_params['cell']
        # We don't worry about arcmin or deg right now - can adjust if needed
        return(float((cell_string).split('arcsec')))
    
    def get_param(self, key=None):
        """
        Get an attribute by its string name.
        """
        if key is None:
            return(None)
        if key in self.clean_params.keys():
            return(self.clean_params[key])
        else:
            raise Exception('clean_params does not have an attribute named "'+str(key)+'"')
            return(None)

    ##################################
    # Hook up to executions of clean #
    ##################################

    def __str__(self):
        """
        Return the clean call as a string that can be pasted into
        CASA. Uses the default CASA version.
        """
        clean_text = 'tclean('
        kwargs_for_clean = self.kwargs_for_clean()
        first = True
        for this_kwarg in kwargs_for_clean:
            if first:
                first = False
            else:
                clean_text += ','
            string_for_kwarg = str(kwargs_for_clean[this_kwarg])
            if type(kwargs_for_clean[this_kwarg]) == type(''):
                string_for_kwarg = '"'+string_for_kwarg+'"'
            clean_text += str(this_kwarg)+'='+string_for_kwarg
        clean_text += ')'
        return(clean_text)
    
    def kwargs_for_clean(
        self,        
        casaversion='5.6',
        ):
        """
        Output keyword arguments to be fed into (t)clean.
        """

        known_versions = ['5.4','5.6']
        if casaversion not in known_versions:
            logger.error('CASA simple version '+str(casaversion)+' not in known version list.')
            return(None)

        #<TODO># Work on adapting to different CASA versions
        
        return(self.clean_params)
        
#endregion
