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
        self.interactive = False
        self.rest = False
        self.logfile = None
        self.clean_mask_file = None

    def execute(self):
        """
        Execute the clean call.
        """
        pass

#endregion

