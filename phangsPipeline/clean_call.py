"""
This is a dummy CleanCall class for dry run only, or to be inheritted by casaImagingRoutines.CleanCall.
"""

#region class CleanCall

class CleanCall:
    
    def __init__(self):
        self.vis = None
        self.antenna = ""
        self.image_root = None
        self.phase_center = ""
        self.image_size = None
        self.cell_size = None
        self.restfreq_ghz = -1.0
        self.calcres = True
        self.calcpsf = True
        self.specmode = 'cube'
        self.deconvolver = 'hogbom'
        self.threshold = '0.0mJy/beam'
        self.scales_as_pix = [0]
        self.scales_as_angle = None
        self.smallscalebias = 0.9
        self.briggs_weight = 0.5
        self.niter = 0
        self.cycle_niter = 200
        self.minpsffraction = 0.5
        self.pblimit = 0.25
        self.uvtaper = None
        self.restoringbeam = 'common'
        self.usemask = 'user'
        self.mask = ''
        self.pbmask = 0.25
        self.interactive = False
        self.reset = False
        self.logfile = None
        self.clean_mask_file = None
        self.gridder = 'mosaic'
        self.outframe = 'lsrk'
        self.veltype = 'radio'
        self.normtype = 'flatnoise'
        self.weighting = 'briggs'
        self.cyclefactor = 3.0
    
    def __str__(self):
        return 'CleanCall: '+str(self.kwargs())
    
    def set(self, key, value, nowarning=False):
        if hasattr(self, key) or nowarning:
            setattr(self, key, value)
        else:
            raise Exception('CleanCall does not have the variable "'+str(key)+'"')
    
    def get(self, key):
        if hasattr(self, key):
            return getattr(self, key)
        else:
            raise Exception('CleanCall does not have the variable "'+str(key)+'"')
            return None
    
    def kwargs(self):
        casa_clean_param_dict = {} #<TODO># In the future if CASA clean parameters changed, we can just change this code.
        casa_clean_param_list = ['vis', 'imagename', 'phasecenter', 'cell', 'imsize', 'gridder', 'specmode', 'restfreq', 
                                 'outframe', 'veltype', 'calcres', 'calcpsf', 'deconvolver', 'scales', 'smallscalebias', 
                                 'pblimit', 'normtype', 'restoringbeam', 'weighting', 'robust', 'uvtaper', 'niter', 
                                 'threshold', 'cycleniter', 'cyclefactor', 'minpsffraction', 'usemask', 'mask', 'pbmask', 
                                 'interactive'] #<TODO># 'antenna' is not considered for now
        
        restfreq_ghz_str = str(self.restfreq_ghz)+'GHz' if self.restfreq_ghz > 0 else ''
        self.set('restfreq',    restfreq_ghz_str,          nowarning=True)
        self.set('imagename',   self.get('image_root'),    nowarning=True)
        self.set('phasecenter', self.get('phase_center'),  nowarning=True)
        self.set('cell',        self.get('cell_size'),     nowarning=True)
        self.set('imsize',      self.get('image_size'),    nowarning=True)
        self.set('scales',      self.get('scales_as_pix'), nowarning=True)
        self.set('robust',      self.get('briggs_weight'), nowarning=True)
        self.set('cycleniter',  self.get('cycle_niter'),   nowarning=True)
        self.set('pbmask',      self.get('pblimit'),       nowarning=True)
        
        for key in casa_clean_param_list:
            casa_clean_param_dict[key] = self.get(key)
        
        casa_clean_param_dict['uvtaper'] = [str(self.uvtaper)+'arcsec',str(self.uvtaper)+'arcsec','0deg'] if self.uvtaper else ''
        
        return casa_clean_param_dict
        
    def execute(self):
        """
        Execute the clean call by calling CASA tclean(**self.kwargs()).
        """
        pass

#endregion






