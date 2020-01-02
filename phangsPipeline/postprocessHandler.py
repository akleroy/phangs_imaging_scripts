"""
Parts of the PHANGS pipeline to handle post-processing cubes.
"""

import os
import glob
import casaCubeRoutines as ccr
import casaMosaicRoutines as cmr
import casaFeatherRoutines as cfr

class PostProcessHandler:
    """
    Class to handle post-processing of ALMA data. Post-processing here
    begins with the results of imaging and proceeds through reduced,
    science-ready data cubes.
    """

    def __init__(
        self,
        key_handler = None,
        dry_run = False,
        quiet = False,
        dochecks = True
        ):
        
        self._quiet = quiet
        self._dochecks = dochecks
        
        self._targets_list = None
        self._mosaics_list = None
        self._line_products_list = None
        self._cont_products_list = None
        self._interf_configs_list = None
        self._feather_configs_list = None

        self._no_cont = False
        self._no_line = False

        if key_handler is not None:
            self._kh = key_handler

        # Initialize the list variables
        self.set_targets(nobuild=True)
        self.set_mosaic_targets(nobuild=True)
        self.set_line_products(nobuild=True)
        self.set_cont_products(nobuild=True)
        self.set_interf_configs(nobuild=True)
        self.set_feather_configs(nobuild=True)

        self._build_lists()

        self.set_dry_run(dry_run)

        return(None)

#region Control what data gets processed

    def set_targets(
        self, 
        first=None, 
        last=None, 
        skip=[], 
        only=[],
        nobuild=False):
        """
        Set conditions on the list of targets to be considered. By
        default, consider all targets.
        """
        self._targets_first = first
        self._targets_last = last
        self._targets_skip = skip
        self._targets_only = only

        if not nobuild:
            self._build_lists()
        return(None)

    def set_mosaic_targets(
        self, 
        first=None, 
        last=None, 
        skip=[], 
        only=[],
        nobuild=False):
        """
        Set conditions on the list of mosaics to be considered. By
        default, consider all mosaics.
        """
        self._mosaics_first = first
        self._mosaics_last = last
        self._mosaics_skip = skip
        self._mosaics_only = only

        if not nobuild:
            self._build_lists()
        return(None)

    def set_line_products(
        self, 
        skip=[], 
        only=[], 
        nobuild=False,
        ):
        """
        Set conditions on the list of line products to be
        considered. By default, consider all products.
        """
        self._lines_skip = skip
        self._lines_only = only

        if not nobuild:
            self._build_lists()
        return(None)

    def set_cont_products(
        self, 
        skip=[], 
        only=[], 
        nobuild=False,
        ):
        """
        Set conditions on the list of continuum products to be
        considered. By default, consider all products.
        """
        self._cont_skip = skip
        self._cont_only = only

        if not nobuild:
            self._build_lists()
        return(None)

    def set_interf_configs(
        self, 
        skip=[], 
        only=[], 
        nobuild=False,
        ):
        """
        Set conditions on the list of interferometric array
        configurations to be considered. By default, consider all
        configurations.
        """
        self._interf_configs_skip = skip
        self._interf_configs_only = only

        if not nobuild:
            self._build_lists()
        return(None)

    def set_feather_configs(
        self, 
        skip=[], 
        only=[],
        nobuild=False,
        ):
        """
        Set conditions on the list of feathered array
        configurations to be considered. By default, consider all
        configurations.
        """
        self._feather_configs_skip = skip
        self._feather_configs_only = only

        if not nobuild:
            self._build_lists()
        return(None)

    def set_no_line(
        self,
        no_line = False):
        """
        Toggle the program to line products.
        """
        self._no_line = no_line

    def set_no_cont(
        self,
        no_cont = False):
        """
        Toggle the program to skip continuum products.
        """
        self._no_cont = no_cont

    def set_dry_run(
        self,
        dry_run = False):
        """
        Toggle the program using a 'dry run', i.e., not actually executing.
        """
        self._dry_run = dry_run

    def set_key_handler(
        self,
        key_handler = None):
        """
        Set the keyhandler being used by the pipeline.
        """
        self._kh = key_handler
        self._build_lists()

    def _build_lists(
        self
        ):
        """
        Build the target lists.
        """

        if self._kh is None:
            if not self._quiet:
                print("Cannot build lists without a keyHandler.")
            return(None)

        self._targets_list = self._kh.get_targets(            
            only = self._targets_only,
            skip = self._targets_skip,
            first = self._targets_first,
            last = self._targets_last,
            )

        self._mosaics_list = self._kh.get_linear_mosaic_targets(            
            only = self._mosaics_only,
            skip = self._mosaics_skip,
            first = self._mosaics_first,
            last = self._mosaics_last,
            )

        if self._no_line:
            self._line_products_list = []
        else:
            self._line_products_list = self._kh.get_line_products(
                only = self._lines_only,
                skip = self._lines_skip,
                )

        if self._no_cont:
            self._cont_products_list = []
        else:
            self._cont_products_list = self._kh.get_continuum_products(
                only = self._cont_only,
                skip = self._cont_skip,
                )

        self._interf_configs_list = self._kh.get_interf_configs(
            only = self._interf_configs_only,
            skip = self._interf_configs_skip,
            )

        self._feather_configs_list = self._kh.get_feather_configs(
            only = self._feather_configs_only,
            skip = self._feather_configs_skip,
            )

    def _all_products(
        self
        ):
        """
        Get a combined list of line and continuum products.
        """
        
        if self._cont_products_list is None:
            if self._line_products_list is None:
                return([])
            else:
                return(self._line_products_list)

        if self._line_products_list is None:
            if self._cont_products_list is None:
                return([])
            else:
                return(self._cont_products_list)

        if len(self._cont_products_list) is 0:
            if self._line_products_list is None:
                return ([])
            else:
                return(self._line_products_list)

        if len(self._line_products_list) is 0:
            if self._cont_products_list is None:
                return([])
            else:
                return(self._cont_products_list)
        
        return(self._line_products_list + self._cont_products_list)

#endregion

#region Master loop

    def _master_loop(
        self,
        do_stage = False,
        do_pbcorr = False,
        do_round = False,
        do_singledish = False,
        do_feather = False,
        do_cleanup = False,
        do_export = False,
        ):
        """
        The master loop that steps over all targets, products, and configurations.
        """              
        
        if self._targets_list is None or self._interf_configs_list is None:            
            if not self._quiet:
                print("MASTER LOOP: Need a target and interferometer configuration list. Returning.")
            return(None)
    
        for this_target in self._targets_list:

            imaging_dir = self._kh.get_imaging_dir_for_target(this_target)
            
            postprocess_dir = self._kh.get_postprocess_dir_for_target(this_target)

            # Loop over line and continuum products

            for this_product in self._all_products():
                
                if this_product is None:
                    continue
                
                # Loop over interferometric configurations only
                
                full_config_list = []
                config_type_list = []

                for this_config in self._interf_configs_list:
                    full_config_list.append(this_config)
                    config_type_list.append('interf')

                for this_config in self._feather_configs_list:
                    full_config_list.append(this_config)
                    config_type_list.append('feather')

                for ii in range(len(full_config_list)):

                    this_config = full_config_list[ii]
                    config_type = config_type_list[ii]

                    orig_file = self._kh.get_cube_filename(
                        target = this_target,
                        config = this_config,
                        product = this_product,
                        ext = None,
                        casa = True,
                        casaext = '.image')

                    pb_file = self._kh.get_cube_filename(
                        target = this_target,
                        config = this_config,
                        product = this_product,
                        ext = None,
                        casa = True,
                        casaext = '.pb')

                    sd_file = self._kh.get_singledish_filename(
                        target = this_target,
                        product = this_product)

                    pbcorr_file = self._kh.get_cube_filename(
                        target = this_target,
                        config = this_config,
                        product = this_product,
                        ext = 'pbcorr',
                        casa = True,
                        casaext = '.image')

                    round_file = self._kh.get_cube_filename(
                        target = this_target,
                        config = this_config,
                        product = this_product,
                        ext = 'round',
                        casa = True,
                        casaext = '.image')

                    pbcorr_round_file = self._kh.get_cube_filename(
                        target = this_target,
                        config = this_config,
                        product = this_product,
                        ext = 'pbcorr_round',
                        casa = True,
                        casaext = '.image')

                    round_k_file = self._kh.get_cube_filename(
                        target = this_target,
                        config = this_config,
                        product = this_product,
                        ext = 'round_k',
                        casa = True,
                        casaext = '.image')

                    pbcorr_round_k_file = self._kh.get_cube_filename(
                        target = this_target,
                        config = this_config,
                        product = this_product,
                        ext = 'pbcorr_round_k',
                        casa = True,
                        casaext = '.image')

                    # Copy the data from the original location to the
                    # postprocessing directories.

                    if do_stage and config_type == 'interf':

                        for fname in [orig_file, pb_file]:
                            
                            infile = imaging_dir + fname
                            outfile = postprocess_dir + fname

                            if not self._quiet:
                                print()
                                print("Staging via ccr.copy_dropdeg: ")
                                print("... from "+infile)
                                print("... to "+outfile)
                                print()
                            
                            if not self._dry_run:
                                ccr.copy_dropdeg(
                                    infile=infile,
                                    outfile=outfile,
                                    overwrite=True,
                                    quiet=self._quiet)

                    # Apply the primary beam correction to the data.

                    if do_pbcorr and config_type == 'interf':

                        infile = postprocess_dir + orig_file
                        outfile = postprocess_dir + pbcorr_file
                        pbfile = postprocess_dir + pb_file

                        if not self._quiet:
                            print()
                            print("Primary beam correcting via ccr.primary_beam_correct: ")
                            print("... from "+infile)
                            print("... to "+outfile)
                            print("... using "+pbfile)
                            print()

                        if not self._dry_run:
                            ccr.primary_beam_corrrect(
                                infile=infile,
                                outfile=outfile,
                                pbfile=pbfile,
                                overwrite=True,
                                quiet=self._quiet)

                    # Convolve the data to have a round beam.

                    if do_round and config_type == 'interf':
                        
                        for fname in [orig_file, pbcorr_file]:
                            
                            infile = postprocess_dir + fname
                            if fname == orig_file:
                                outfile = postprocess_dir + round_file

                            if fname == pbcorr_file:
                                outfile = postprocess_dir + pbcorr_round_file

                            if not self._quiet:
                                print()
                                print("Convolving to have a round beam via ccr.convolve_to_round_beam: ")
                                print("... from "+infile)
                                print("... to "+outfile)
                                print()
                        
                            if not self._dry_run:
                                ccr.convolve_to_round_beam(
                                    infile=infile,
                                    outfile=outfile,
                                    overwrite=True,
                                    quiet=self._quiet)

                    # Stage the singledish data for feathering

                    if do_singledish and config_type == 'interf':

                        template = postprocess_dir + orig_file
                        infile = sdfile
                        outfile = postprocess_dir + prepped_sd_file

                        if not self._quiet:
                            print()
                            print("Prepping single dish target for feather via cfr.prep_sd_for_feather: ")
                            print("... original file "+infile)
                            print("... prepped file "+outfile)
                            print("... template "+template)
                            print()
                        
                        if not self._dry_run:
                            cfr.prep_sd_for_feather(
                                sdfile_in=infile,
                                sdfile_out=outfile,
                                interf_file=template,
                                doimport=True,
                                checkunits=True,
                                doalign=True,
                                overwrite=True,
                                quiet=self._quiet)
                            
                    # Feather the single dish and interferometer data

                    if do_feather and config_type == 'interf':

                        pass
                                
                    # "Clean up" - meaning reduce volume and change
                    # units.

                    if do_cleanup:

                        pass

                    # Export to FITS and clean up output

                    if do_export:

                        pass

        return()

#endregion

#region Stage and correct data

    def stage_interferometer_data(
        self
        ):
        """
        Copy the interferometer data to start post-processing.
        """              
        
        self._master_loop(do_stage=True)

        return()

    def primary_beam_correct(
        self
        ):
        """
        Apply primary beam correction to the interferometer data.
        """
 
       self._master_loop(do_pbcorr=True)

       return()

    def convolve_to_round_beam(
        self
        ):
        """
        Convolve the interferometer data to have a round beam.
        """                    

        self._master_loop(do_round=True)

        return()
        
    def stage_singledish_data(
        self
        ):
        """
        Copy the single dish data and align it to the releveant
        interferometer data set.
        """
        pass

#endregion

#region Mosaic data together.

    def convolve_for_mosaic(
        self
        ):
        """
        Convolve all interfereometer data in a mosaic to share a
        common beam.
        """
        pass

    def align_for_mosaic(
        self
        ):
        """
        Build a shared header for a mosaic and align all data to that
        new header.
        """
        pass

    def mosaic_interf_data(
        self
        ):
        """
        Mosaic the interferometric data.
        """
        pass

    def mosaic_singledish_data(
        self
        ):
        """
        Mosaic the single dish data.
        """
        pass    

#endregion

#region Feather interferometer and single dish data.

    def feather_cubes(
        self
        ):
        """
        Feather the single dish and interferometer data.
        """
        pass

#endregion

# This is all earlier...

def phangs_stage_singledish(
    gal=None, product=None, root_dir=None, 
    overwrite=False):
    """
    Copy the single dish data for further processing
    """

    if gal is None or product is None or \
            root_dir is None:
        print("Missing required input.")
        return
    
    sdk = read_singledish_key()
    if (gal in sdk.keys()) == False:
        print(gal+" not found in single dish key.")
        return
    
    this_key = sdk[gal]
    if (product in this_key.keys()) == False:
        print(product+" not found in single dish key for "+gal)
        return
    
    sdfile_in = this_key[product]
    
    sdfile_out = root_dir+'raw/'+gal+'_tp_'+product+'.image'    

    print("... importing single dish data for "+sdfile_in)

    importfits(fitsimage=sdfile_in, imagename=sdfile_out+'.temp',
               zeroblanks=True, overwrite=overwrite)

    if overwrite:
        os.system('rm -rf '+sdfile_out)
    imsubimage(imagename=sdfile_out+'.temp', outfile=sdfile_out,
               dropdeg=True)
    os.system('rm -rf '+sdfile_out+'.temp')

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# BASIC IMAGE PROCESSING STEPS
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

def phangs_primary_beam_correct(
    gal=None, array=None, product=None, root_dir=None, 
    cutoff=0.25, overwrite=False):
    """
    Construct primary-beam corrected images using PHANGS naming conventions.
    """

    input_dir = root_dir+'raw/'
    input_cube_name = input_dir+gal+'_'+array+'_'+product+'.image'
    input_pb_name = input_dir+gal+'_'+array+'_'+product+'.pb'
    output_dir = root_dir+'process/'
    output_cube_name = output_dir+gal+'_'+array+'_'+product+'_pbcorr.image'

    print("")
    print("... producing a primary beam corrected image for "+input_cube_name)
    print("")
    
    primary_beam_correct(
        infile=input_cube_name, 
        pbfile=input_pb_name,
        outfile=output_cube_name,
        cutoff=cutoff, overwrite=overwrite)

def primary_beam_correct(
    infile=None, pbfile=None, outfile=None, 
    cutoff=0.25, overwrite=False):
    """
    Construct a primary-beam corrected image.
    """

    if infile is None or pbfile is None or outfile is None:
        print("Missing required input.")
        return

    if os.path.isdir(infile) == False:
        print("Input file missing - "+infile)
        return

    if os.path.isdir(pbfile) == False:
        print("Primary beam file missing - "+pbfile)
        return

    if overwrite:
        os.system('rm -rf '+outfile)

    impbcor(imagename=infile, pbimage=pbfile, outfile=outfile, cutoff=cutoff)

def phangs_convolve_to_round_beam(
    gal=None, array=None, product=None, root_dir=None, 
    force_beam=None, overwrite=False):
    """
    Construct primary-beam corrected images using PHANGS naming
    conventions. Runs on both primary beam corrected cube and flat
    cube, forcing the same beam for both.
    """

    if gal is None or array is None or product is None or \
            root_dir is None:
        print("Missing required input.")
        return
    
    input_dir = root_dir+'process/'
    input_cube_name = input_dir+gal+'_'+array+'_'+product+'_pbcorr.image'
    output_dir = root_dir+'process/'
    output_cube_name = output_dir+gal+'_'+array+'_'+product+'_pbcorr_round.image'

    print("")
    print("... convolving to a round beam for "+input_cube_name)
    print("")

    round_beam = \
        convolve_to_round_beam(
        infile=input_cube_name,
        outfile=output_cube_name,
        force_beam=force_beam,
        overwrite=overwrite)

    print("")
    print("... found beam of "+str(round_beam)+" arcsec. Forcing flat cube to this.")
    print("")

    input_dir = root_dir+'raw/'
    input_cube_name = input_dir+gal+'_'+array+'_'+product+'.image'
    output_dir = root_dir+'process/'
    output_cube_name = output_dir+gal+'_'+array+'_'+product+'_flat_round.image'    

    convolve_to_round_beam(
        infile=input_cube_name,
        outfile=output_cube_name,
        force_beam=round_beam,
        overwrite=overwrite)
    
def convolve_to_round_beam(
    infile=None, outfile=None, force_beam=None, overwrite=False):
    """
    Convolve supplied image to have a round beam.
    """
    
    if infile is None or outfile is None:
        print("Missing required input.")
        return

    if os.path.isdir(infile) == False:
        print("Input file missing - "+infile)
        return    

    if force_beam is None:
        hdr = imhead(infile)

        if (hdr['axisunits'][0] != 'rad'):
            print("ERROR: Based on CASA experience. I expected units of radians.")
            print("I did not find this. Returning. Adjust code or investigate file "+infile)
            return
        pixel_as = abs(hdr['incr'][0]/np.pi*180.0*3600.)

        if (hdr['restoringbeam']['major']['unit'] != 'arcsec'):
            print("ERROR: Based on CASA experience. I expected units of arcseconds for the beam.")
            print("I did not find this. Returning. Adjust code or investigate file "+infile)
            return    
        bmaj = hdr['restoringbeam']['major']['value']    
        target_bmaj = np.sqrt((bmaj)**2+(2.0*pixel_as)**2)
    else:
        target_bmaj = force_beam

    imsmooth(imagename=infile,
             outfile=outfile,
             targetres=True,
             major=str(target_bmaj)+'arcsec',
             minor=str(target_bmaj)+'arcsec',
             pa='0.0deg',
             overwrite=overwrite
             )

    return target_bmaj

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# ROUTINES FOR FEATHERING THE DATA
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

def prep_for_feather(
    gal=None, array=None, product=None, root_dir=None, 
    overwrite=False):
    """
    Prepare the single dish data for feathering
    """
    
    if gal is None or array is None or product is None or \
            root_dir is None:
        print("Missing required input.")
        return    

    sdfile_in = root_dir+'raw/'+gal+'_tp_'+product+'.image'
    interf_in = root_dir+'process/'+gal+'_'+array+'_'+product+'_flat_round.image'
    pbfile_name = root_dir+'raw/'+gal+'_'+array+'_'+product+'.pb'    

    if (os.path.isdir(sdfile_in) == False):
        print("Single dish file not found: "+sdfile_in)
        return

    if (os.path.isdir(interf_in) == False):
        print("Interferometric file not found: "+interf_in)
        return

    if (os.path.isdir(pbfile_name) == False):
        print("Primary beam file not found: "+pbfile_name)
        return

    # Align the relevant TP data to the product.
    sdfile_out = root_dir+'process/'+gal+'_tp_'+product+'_align_'+array+'.image'
    imregrid(imagename=sdfile_in,
             template=interf_in,
             output=sdfile_out,
             asvelocity=True,
             axes=[-1],
             interpolation='cubic',
             overwrite=overwrite)

    # Taper the TP data by the primary beam.
    taperfile_out = root_dir+'process/'+gal+'_tp_'+product+'_taper_'+array+'.image'
    if overwrite:
        os.system('rm -rf '+taperfile_out)

    impbcor(imagename=sdfile_out, 
            pbimage=pbfile_name, 
            outfile=taperfile_out, 
            mode='multiply',
            stokes='I')

    return

def phangs_feather_data(
    gal=None, array=None, product=None, root_dir=None, 
    cutoff=-1,overwrite=False):
    """
    Feather the interferometric and total power data.
    """

    if gal is None or array is None or product is None or \
            root_dir is None:
        print("Missing required input.")
        return    

    sdfile_in = root_dir+'process/'+gal+'_tp_'+product+'_taper_'+array+'.image'
    interf_in = root_dir+'process/'+gal+'_'+array+'_'+product+'_flat_round.image'
    pbfile_name = root_dir+'raw/'+gal+'_'+array+'_'+product+'.pb' 

    if (os.path.isdir(sdfile_in) == False):
        print("Single dish file not found: "+sdfile_in)
        return
        
    if (os.path.isdir(interf_in) == False):
        print("Interferometric file not found: "+interf_in)
        return

    if (os.path.isdir(pbfile_name) == False):
        print("Primary beam file not found: "+pbfile_name)
        return

    # Feather the inteferometric and "flat" TP data.
    outfile_name = root_dir+'process/'+gal+'_'+array+'+tp_'+product+ \
        '_flat_round.image'

    if overwrite:        
        os.system('rm -rf '+outfile_name)
    os.system('rm -rf '+outfile_name+'.temp')
    feather(imagename=outfile_name+'.temp',
            highres=interf_in,
            lowres=sdfile_in,
            sdfactor=1.0,
            lowpassfiltersd=False)
    imsubimage(imagename=outfile_name+'.temp', outfile=outfile_name,
               dropdeg=True)
    os.system('rm -rf '+outfile_name+'.temp')
    infile_name = outfile_name

    # Primary beam correct the feathered data.
    outfile_name = root_dir+'process/'+gal+'_'+array+'+tp_'+product+ \
        '_pbcorr_round.image'
    
    if overwrite:        
        os.system('rm -rf '+outfile_name)

    print(infile_name)
    print(pbfile_name)
    impbcor(imagename=infile_name,
            pbimage=pbfile_name, 
            outfile=outfile_name, 
            mode='divide', cutoff=cutoff)

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# CLEANUP
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

def convert_jytok(
    infile=None, outfile=None, overwrite=False, inplace=False):
    """
    Convert a cube from Jy/beam to K.
    """

    c = 2.99792458e10
    h = 6.6260755e-27
    kb = 1.380658e-16

    if infile is None or (outfile is None and inplace==False):
        print("Missing required input.")
        return
    
    if os.path.isdir(infile) == False:
        print("Input file not found: "+infile)
        return
    
    if inplace == False:
        if overwrite:
            os.system('rm -rf '+outfile)
        
        if os.path.isdir(outfile):
            print("Output file already present: "+outfile)
            return

        os.system('cp -r '+infile+' '+outfile)
        target_file = outfile
    else:
        target_file = infile

    hdr = imhead(target_file, mode='list')
    unit = hdr['bunit']
    if unit != 'Jy/beam':
        print("Unit is not Jy/beam. Returning.")
        return

    #restfreq_hz = hdr['restfreq'][0]

    if hdr['cunit3'] != 'Hz':
        print("I expected frequency as the third axis but did not find it.")
        print("Returning.")
        return
    
    crpix3 = hdr['crpix3']
    cdelt3 = hdr['cdelt3']
    crval3 = hdr['crval3']
    naxis3 = hdr['shape'][2]
    faxis_hz = (np.arange(naxis3)+1.-crpix3)*cdelt3+crval3
    freq_hz = np.mean(faxis_hz)
    
    bmaj_unit = hdr['beammajor']['unit']
    if bmaj_unit != 'arcsec':
        print("Beam unit is not arcsec, which I expected. Returning.")
        print("Unit instead is "+bmaj_unit)
        return    
    bmaj_as = hdr['beammajor']['value']
    bmin_as = hdr['beamminor']['value']
    bmaj_sr = bmaj_as/3600.*np.pi/180.
    bmin_sr = bmin_as/3600.*np.pi/180.
    beam_in_sr = np.pi*(bmaj_sr/2.0*bmin_sr/2.0)/np.log(2)
    
    jtok = c**2 / beam_in_sr / 1e23 / (2*kb*freq_hz**2)

    myia = au.createCasaTool(iatool)
    myia.open(target_file)
    vals = myia.getchunk()
    vals *= jtok
    myia.putchunk(vals)
    myia.setbrightnessunit("K")
    myia.close()

    imhead(target_file, mode='put', hdkey='JTOK', hdvalue=jtok)

    return

def trim_cube(    
    infile=None, outfile=None, overwrite=False, inplace=False, min_pixperbeam=3):
    """
    Trim and rebin a cube to smaller size.
    """
    
    if infile is None or outfile is None:
        print("TRIM_CUBE: Missing required input.")
        return
    
    if os.path.isdir(infile) == False:
        print("TRIM_CUBE: Input file not found: "+infile)
        return

    # First, rebin if needed
    hdr = imhead(infile)
    if (hdr['axisunits'][0] != 'rad'):
        print("ERROR: Based on CASA experience. I expected units of radians.")
        print("I did not find this. Returning. Adjust code or investigate file "+infile)
        return

    pixel_as = abs(hdr['incr'][0]/np.pi*180.0*3600.)

    if (hdr['restoringbeam']['major']['unit'] != 'arcsec'):
        print("ERROR: Based on CASA experience. I expected units of arcseconds for the beam.")
        print("I did not find this. Returning. Adjust code or investigate file "+infile)
        return    
    bmaj = hdr['restoringbeam']['major']['value']    
    
    pix_per_beam = bmaj*1.0 / pixel_as*1.0
    
    if pix_per_beam > 6:
        imrebin(
            imagename=infile,
            outfile=outfile+'.temp',
            factor=[2,2,1],
            crop=True,
            dropdeg=True,
            overwrite=overwrite,
            )
    else:
        os.system('cp -r '+infile+' '+outfile+'.temp')

    # Figure out the extent of the image inside the cube
    myia = au.createCasaTool(iatool)
    myia.open(outfile+'.temp')
    mask = myia.getchunk(getmask=True)    
    myia.close()

    this_shape = mask.shape

    mask_spec_x = np.sum(np.sum(mask*1.0,axis=2),axis=1) > 0
    pad = 0
    xmin = np.max([0,np.min(np.where(mask_spec_x))-pad])
    xmax = np.min([np.max(np.where(mask_spec_x))+pad,mask.shape[0]-1])

    mask_spec_y = np.sum(np.sum(mask*1.0,axis=2),axis=0) > 0
    ymin = np.max([0,np.min(np.where(mask_spec_y))-pad])
    ymax = np.min([np.max(np.where(mask_spec_y))+pad,mask.shape[1]-1])

    mask_spec_z = np.sum(np.sum(mask*1.0,axis=0),axis=0) > 0
    zmin = np.max([0,np.min(np.where(mask_spec_z))-pad])
    zmax = np.min([np.max(np.where(mask_spec_z))+pad,mask.shape[2]-1])
    
    box_string = ''+str(xmin)+','+str(ymin)+','+str(xmax)+','+str(ymax)
    chan_string = ''+str(zmin)+'~'+str(zmax)

    print("... ... ... box selection: "+box_string)
    print("... ... ... channel selection: "+chan_string)

    if overwrite:
        os.system('rm -rf '+outfile)
        imsubimage(
        imagename=outfile+'.temp',
        outfile=outfile,
        box=box_string,
        chans=chan_string,
        )
    
    os.system('rm -rf '+outfile+'.temp')
    
def phangs_cleanup_cubes(
        gal=None, array=None, product=None, root_dir=None, 
        overwrite=False, min_pixeperbeam=3, roundbeam_tol=0.01, 
        vstring=''):
    """
    Clean up cubes.
    """

    if gal is None or array is None or product is None or \
            root_dir is None:
        print("Missing required input.")
        return

    for this_ext in ['flat', 'pbcorr']:

        root = root_dir+'process/'+gal+'_'+array+'_'+product+'_'+this_ext
        infile = root+'_round.image'
        outfile = root+'_round_k.image'
        outfile_fits = root+'_round_k.fits'
    
        if os.path.isdir(infile) == False:
            print("File does not exist: "+infile)
            print("Returning.")
            return

        # Trim the cube to a smaller size and rebin as needed

        trim_cube(infile=infile, outfile=outfile, 
                  overwrite=overwrite, inplace=False,
                  min_pixperbeam=min_pixeperbeam)

        # Convert to Kelvin

        convert_jytok(infile=outfile, inplace=True)

        # Export to FITS
    
        exportfits(imagename=outfile, fitsimage=outfile_fits,
                   velocity=True, overwrite=True, dropstokes=True, 
                   dropdeg=True, bitpix=-32)
    
        # Clean up headers

        hdu = pyfits.open(outfile_fits)

        hdr = hdu[0].header
        data = hdu[0].data

        for card in ['BLANK','DATE-OBS','OBSERVER','O_BLANK','O_BSCALE',
                     'O_BZERO','OBSRA','OBSDEC','OBSGEO-X','OBSGEO-Y','OBSGEO-Z',
                     'DISTANCE']:
            if card in hdr.keys():
                hdr.remove(card)
            
        while 'HISTORY' in hdr.keys():
            hdr.remove('HISTORY')

        hdr.add_history('This cube was produced by the PHANGS-ALMA pipeline.')
        hdr.add_history('PHANGS-ALMA Pipeline version ' + pipeVer)
        if vstring != '':
            hdr.add_history('This is part of data release '+vstring)

        hdr['OBJECT'] = dir_for_gal(gal)

        if vstring == '':
            hdr['ORIGIN'] = 'PHANGS-ALMA'
        else:
            hdr['ORIGIN'] = 'PHANGS-ALMA '+vstring

        datamax = np.nanmax(data)
        datamin = np.nanmin(data)
        hdr['DATAMAX'] = datamax
        hdr['DATAMIN'] = datamin

        # round the beam if it lies within the specified tolerance

        bmaj = hdr['BMAJ']
        bmin = hdr['BMIN']
        if bmaj != bmin:
            frac_dev = np.abs(bmaj-bmin)/bmaj
            if frac_dev <= roundbeam_tol:
                print("Rounding beam.")
                hdr['BMAJ'] = bmaj
                hdr['BMIN'] = bmaj
                hdr['BPA'] = 0.0
            else:
                print("Beam too asymmetric to round.")
                print("... fractional deviation: "+str(frac_dev))
                
        hdu.writeto(outfile_fits, clobber=True)
        
        return

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# LINEAR MOSAICKING ROUTINES
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

def phangs_common_res_for_mosaic(
    gal=None, array=None, product=None, root_dir=None, 
    overwrite=False, target_res=None):
    """
    Convolve multi-part cubes to a common res for mosaicking.
    """
    
    if gal is None or array is None or product is None or \
            root_dir is None:
        print("Missing required input.")
        return    

    # Look up parts
    this_mosaic_key = mosaic_key()
    if (gal in this_mosaic_key.keys()) == False:
        print("Galaxy "+gal+" not in mosaic key.")
        return
    parts = this_mosaic_key[gal]

    for this_ext in ['flat_round', 'pbcorr_round']:           

        infile_list = []
        outfile_list = []
        input_dir = root_dir+'process/'
        for this_part in parts:
            infile = input_dir+this_part+'_'+array+'_'+product+'_'+this_ext+'.image'
            infile_list.append(infile)
            outfile = input_dir+this_part+'_'+array+'_'+product+'_'+this_ext+'_tomerge.image'
            outfile_list.append(outfile)

        common_res_for_mosaic(
            infile_list=infile_list,
            outfile_list=outfile_list,
            overwrite=overwrite, target_res=target_res)

def common_res_for_mosaic(
    infile_list = None, outfile_list = None,
    overwrite=False, target_res=None):
    """
    Convolve multi-part cubes to a common res for mosaicking.
    """
    
    if (infile_list is None) or \
            (outfile_list is None):
        print("Missing required input.")
        return    
    
    if len(infile_list) != len(outfile_list):
        print("Mismatch in input lists.")
        return    

    for this_file in infile_list:
        if os.path.isdir(this_file) == False:
            print("File not found "+this_file)
            return
    
    # Figure out target resolution if it is not supplied

    if target_res is None:
        print("Calculating target resolution ... ")

        bmaj_list = []
        pix_list = []

        for infile in infile_list:
            print("Checking "+infile)

            hdr = imhead(infile)

            if (hdr['axisunits'][0] != 'rad'):
                print("ERROR: Based on CASA experience. I expected units of radians.")
                print("I did not find this. Returning. Adjust code or investigate file "+infile)
                return
            this_pixel = abs(hdr['incr'][0]/np.pi*180.0*3600.)

            if (hdr['restoringbeam']['major']['unit'] != 'arcsec'):
                print("ERROR: Based on CASA experience. I expected units of arcseconds for the beam.")
                print("I did not find this. Returning. Adjust code or investigate file "+infile)
                return
            this_bmaj = hdr['restoringbeam']['major']['value']

            bmaj_list.append(this_bmaj)
            pix_list.append(this_pixel)
        
        max_bmaj = np.max(bmaj_list)
        max_pix = np.max(pix_list)
        target_bmaj = np.sqrt((max_bmaj)**2+(2.0*max_pix)**2)
    else:
        target_bmaj = force_beam

    for ii in range(len(infile_list)):
        this_infile = infile_list[ii]
        this_outfile = outfile_list[ii]
        print("Convolving "+this_infile+' to '+this_outfile)
        
        imsmooth(imagename=this_infile,
             outfile=this_outfile,
             targetres=True,
             major=str(target_bmaj)+'arcsec',
             minor=str(target_bmaj)+'arcsec',
             pa='0.0deg',
             overwrite=overwrite
             )

    return target_bmaj

def build_common_header(
    infile_list = None, 
    ra_ctr = None, dec_ctr = None,
    delta_ra = None, delta_dec = None):
    """
    Build a target header to be used as a template by imregrid.
    """
    
    if infile_list is None:
        print("Missing required input.")
        return    

    # Logic to determine tuning parameters here if they aren't passed.

    if os.path.isdir(infile_list[0]) == False:
        print("File not found "+infile_list[0])
        print("Returning.")
        return None
    target_hdr = imregrid(infile_list[0], template='get')
    
    # N.B. Could put a lot of general logic here, but we are usually
    # working in a pretty specific case.

    if (target_hdr['csys']['direction0']['units'][0] != 'rad') or \
            (target_hdr['csys']['direction0']['units'][1] != 'rad'):
        print("ERROR: Based on CASA experience. I expected pixel units of radians.")
        print("I did not find this. Returning. Adjust code or investigate file "+infile_list[0])
        return

    # Put in our target values for the center after converting to radians
    ra_ctr_in_rad = ra_ctr * np.pi / 180.
    dec_ctr_in_rad = dec_ctr * np.pi / 180.

    target_hdr['csys']['direction0']['crval'][0] = ra_ctr_in_rad
    target_hdr['csys']['direction0']['crval'][1] = dec_ctr_in_rad

    # Adjust the size and central pixel
    
    ra_pix_in_as = np.abs(target_hdr['csys']['direction0']['cdelt'][0]*180./np.pi*3600.)
    dec_pix_in_as = np.abs(target_hdr['csys']['direction0']['cdelt'][1]*180./np.pi*3600.)
    ra_axis_size = np.ceil(delta_ra / ra_pix_in_as)
    new_ra_ctr_pix = ra_axis_size/2.0
    dec_axis_size = np.ceil(delta_dec / dec_pix_in_as)
    new_dec_ctr_pix = dec_axis_size/2.0
    
    target_hdr['csys']['direction0']['crpix'][0] = new_ra_ctr_pix
    target_hdr['csys']['direction0']['crpix'][1] = new_dec_ctr_pix
    
    if ra_axis_size > 1e4 or dec_axis_size > 1e4:
        print("WARNING! This is a very big image you plan to create.")
        print(ra_axis_size, " x ", dec_axis_size)
        test = raw_input("Continue? Hit [y] if so.")
        if test != 'y':
            return

    target_hdr['shap'][0] = int(ra_axis_size)
    target_hdr['shap'][1] = int(dec_axis_size)
    
    return(target_hdr)

def align_for_mosaic(
    infile_list = None, outfile_list = None,
    overwrite=False, target_hdr=None):
    """
    Align a list of files to a target coordinate system.
    """

    if infile_list is None or outfile_list is None or \
            target_hdr is None:
        print("Missing required input.")
        return    

    for ii in range(len(infile_list)):
        this_infile = infile_list[ii]
        this_outfile = outfile_list[ii]        

        if os.path.isdir(this_infile) == False:
            print("File "+this_infile+" not found. Continuing.")
            continue

        imregrid(imagename=this_infile,
                 template=target_hdr,
                 output=this_outfile,
                 asvelocity=True,
                 axes=[-1],
                 interpolation='cubic',
                 overwrite=overwrite)

    return

def phangs_align_for_mosaic(
    gal=None, array=None, product=None, root_dir=None, 
    overwrite=False, target_hdr=None):
    """
    Convolve multi-part cubes to a common res for mosaicking.
    """

    if gal is None or array is None or product is None or \
            root_dir is None:
        print("Missing required input.")
        return    
    
    # Look up parts

    this_mosaic_key = mosaic_key()
    if (gal in this_mosaic_key.keys()) == False:
        print("Galaxy "+gal+" not in mosaic key.")
        return
    parts = this_mosaic_key[gal]

    # Read the key that defines the extent and center of the mosaic
    # manually. We will use this to figure out the target header.

    multipart_key = read_multipart_key()
    if (gal in multipart_key.keys()) == False:
        print("Galaxy "+gal+" not in multipart key.")
        print("... working on a general header construction algorithm.")
        print("... for now, go enter a center and size into the multipart key:")        
        print("... multipart_fields.txt ")
        return
    this_ra_ctr = multipart_key[gal]['ra_ctr_deg']
    this_dec_ctr = multipart_key[gal]['dec_ctr_deg']
    this_delta_ra = multipart_key[gal]['delta_ra_as']
    this_delta_dec = multipart_key[gal]['delta_dec_as']

    for this_ext in ['flat_round', 'pbcorr_round']:           

        # Align data

        infile_list = []
        outfile_list = []
        input_dir = root_dir+'process/'
        output_dir = root_dir+'process/'
        for this_part in parts:
            infile = input_dir+this_part+'_'+array+'_'+product+'_'+this_ext+'_tomerge.image'
            infile_list.append(infile)
            outfile = output_dir+this_part+'_'+array+'_'+product+'_'+this_ext+'_onmergegrid.image'
            outfile_list.append(outfile)

        # Work out the target header if it does not exist yet.

        if target_hdr is None:
            target_hdr = \
                build_common_header(
                infile_list = infile_list, 
                ra_ctr = this_ra_ctr, dec_ctr = this_dec_ctr,
                delta_ra = this_delta_ra, delta_dec = this_delta_dec)
            
        align_for_mosaic(
            infile_list = infile_list, 
            outfile_list = outfile_list,
            overwrite=overwrite, target_hdr=target_hdr)

        # Align primary beam images, too, to use as weight.

        infile_list = []
        outfile_list = []
        input_dir = root_dir+'raw/'
        output_dir = root_dir+'process/'
        for this_part in parts:
            input_array = array
            if array == '7m+tp':
                input_array = '7m'
            if array == '12m+7m+tp':
                input_array = '12m+7m'
            infile = input_dir+this_part+'_'+input_array+'_'+product+'.pb'
            infile_list.append(infile)
            outfile = output_dir+this_part+'_'+array+'_'+product+'_'+this_ext+'_mergeweight.image'
            outfile_list.append(outfile)

        align_for_mosaic(
            infile_list = infile_list, 
            outfile_list = outfile_list,
            overwrite=overwrite, target_hdr=target_hdr)

def mosaic_aligned_data(
    infile_list = None, weightfile_list = None,
    outfile = None, overwrite=False):
    """
    Combine a list of aligned data with primary-beam (i.e., inverse
    noise) weights using simple linear mosaicking.
    """

    if infile_list is None or weightfile_list is None or \
            outfile is None:
        print("Missing required input.")
        return    

    sum_file = outfile+'.sum'
    weight_file = outfile+'.weight'

    if (os.path.isdir(outfile) or \
            os.path.isdir(sum_file) or \
            os.path.isdir(weight_file)) and \
            (overwrite == False):
        print("Output file present and overwrite off.")
        print("Returning.")
        return

    if overwrite:
        os.system('rm -rf '+outfile+'.temp')
        os.system('rm -rf '+outfile)
        os.system('rm -rf '+sum_file)
        os.system('rm -rf '+weight_file)
        os.system('rm -rf '+outfile+'.mask')

    imlist = infile_list[:]
    imlist.extend(weightfile_list)
    n_image = len(infile_list)
    lel_exp_sum = ''
    lel_exp_weight = ''
    first = True
    for ii in range(n_image):
        this_im = 'IM'+str(ii)
        this_wt = 'IM'+str(ii+n_image)
        this_lel_sum = '('+this_im+'*'+this_wt+'*'+this_wt+')'
        this_lel_weight = '('+this_wt+'*'+this_wt+')'
        if first:
            lel_exp_sum += this_lel_sum
            lel_exp_weight += this_lel_weight
            first=False
        else:
            lel_exp_sum += '+'+this_lel_sum
            lel_exp_weight += '+'+this_lel_weight

    immath(imagename = imlist, mode='evalexpr',
           expr=lel_exp_sum, outfile=sum_file,
           imagemd = imlist[0])
    
    myia = au.createCasaTool(iatool)
    myia.open(sum_file)
    myia.set(pixelmask=1)
    myia.close()

    immath(imagename = imlist, mode='evalexpr',
           expr=lel_exp_weight, outfile=weight_file)
    myia.open(weight_file)
    myia.set(pixelmask=1)
    myia.close()

    immath(imagename = [sum_file, weight_file], mode='evalexpr',
           expr='iif(IM1 > 0.0, IM0/IM1, 0.0)', outfile=outfile+'.temp',
           imagemd = sum_file)

    immath(imagename = weight_file, mode='evalexpr',
           expr='iif(IM0 > 0.0, 1.0, 0.0)', outfile=outfile+'.mask')

    imsubimage(imagename=outfile+'.temp', outfile=outfile,
               mask='"'+outfile+'.mask"', dropdeg=True)



    return

def phangs_mosaic_data(
    gal=None, array=None, product=None, root_dir=None, 
    overwrite=False):
    """
    Linearly mosaic multipart cubes.
    """

    if gal is None or array is None or product is None or \
            root_dir is None:
        print("Missing required input.")
        return    
    
    # Look up parts

    this_mosaic_key = mosaic_key()
    if (gal in this_mosaic_key.keys()) == False:
        print("Galaxy "+gal+" not in mosaic key.")
        return
    parts = this_mosaic_key[gal]

    for this_ext in ['flat_round', 'pbcorr_round']:           

        infile_list = []
        wtfile_list = []
        input_dir = root_dir+'process/'
        output_dir = root_dir+'process/'
        for this_part in parts:
            infile = input_dir+this_part+'_'+array+'_'+product+'_'+this_ext+'_onmergegrid.image'
            wtfile = output_dir+this_part+'_'+array+'_'+product+'_'+this_ext+'_mergeweight.image'
            if (os.path.isdir(infile) == False):
                print("Missing "+infile)
                print("Skipping.")
                continue
            if (os.path.isdir(wtfile) == False):
                print("Missing "+wtfile)
                print("Skipping.")
                continue
            infile_list.append(infile)
            wtfile_list.append(wtfile)

        if len(infile_list) == 0:
            print("No files to include in the mosaic. Returning.")
            return

        outfile = output_dir+gal+'_'+array+'_'+product+'_'+this_ext+'.image'
    
        mosaic_aligned_data(
            infile_list = infile_list, weightfile_list = wtfile_list,
            outfile = outfile, overwrite=overwrite)
