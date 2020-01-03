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

        self._feather_method = 'pbcorr'

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
        self._build_lists()

    def set_no_cont(
        self,
        no_cont = False):
        """
        Toggle the program to skip continuum products.
        """
        self._no_cont = no_cont
        self._build_lists()

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

    def set_feather_method(
        self,
        method='pbcorr'
        ):
        """
        Set the approach to feathering used in the pipeline.
        """
        valid_choices = ['pbcorr','apodize']
        if method.lower() not in valid_choices:
            if not self._quiet:
                print("Not a valid feather method: "+method)
            return(False)
        self._feather_method = method
        return(True)

#endregion

#region Behind the scenes infrastructure and book keeping.

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
        do_compress = False,
        do_convert = False,
        do_export = False,
        ):
        """
        The master loop that steps over all targets, products, and
        configurations. This is the core of the post-processing
        pipeline. It calls the various CASA routines and uses the
        keyHandler to build various file names. It's best accessed via
        the other programs.
        """              
        
        this_stub = 'POSTPROCESS MASTER LOOP: '

        if self._targets_list is None or self._interf_configs_list is None:            
            if not self._quiet:
                print(this_stub+"Need a target and interferometer configuration list.")
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

                    # Original cube and primary beam file
                    
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

                    # Original single dish file

                    orig_sd_file = self._kh.get_singledish_filename(
                        target = this_target,
                        product = this_product)

                    # Primary beam corrected file

                    pbcorr_file = self._kh.get_cube_filename(
                        target = this_target,
                        config = this_config,
                        product = this_product,
                        ext = 'pbcorr',
                        casa = True,
                        casaext = '.image')

                    # Files with round beams

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

                    # Aligned and imported single dish file

                    prepped_sd_file = self._kh.get_cube_filename(
                        target = this_target,
                        config = this_config,
                        product = this_product,
                        ext = 'singledish',
                        casa = True,
                        casaext = '.image')

                    # Compressed file with edges trimmed off

                    trimmed_file = self._kh.get_cube_filename(
                        target = this_target,
                        config = this_config,
                        product = this_product,
                        ext = 'trimmed',
                        casa = True,
                        casaext = '.image')

                    pbcorr_trimmed_file = self._kh.get_cube_filename(
                        target = this_target,
                        config = this_config,
                        product = this_product,
                        ext = 'pbcorr_trimmed',
                        casa = True,
                        casaext = '.image')

                    trimmed_pb_file = self._kh.get_cube_filename(
                        target = this_target,
                        config = this_config,
                        product = this_product,
                        ext = 'trimmed',
                        casa = True,
                        casaext = '.pb')

                    # Files converted to Kelvin and FITS counterparts

                    trimmed_k_file = self._kh.get_cube_filename(
                        target = this_target,
                        config = this_config,
                        product = this_product,
                        ext = 'trimmed_k',
                        casa = True,
                        casaext = '.image')

                    trimmed_k_fits = self._kh.get_cube_filename(
                        target = this_target,
                        config = this_config,
                        product = this_product,
                        ext = 'trimmed_k',
                        casa = False)

                    pbcorr_trimmed_k_file = self._kh.get_cube_filename(
                        target = this_target,
                        config = this_config,
                        product = this_product,
                        ext = 'pbcorr_trimmed_k',
                        casa = True,
                        casaext = '.image')

                    pbcorr_trimmed_k_fits = self._kh.get_cube_filename(
                        target = this_target,
                        config = this_config,
                        product = this_product,
                        ext = 'pbcorr_trimmed_k',
                        casa = False)

                    trimmed_pb_fits = self._kh.get_cube_filename(
                        target = this_target,
                        config = this_config,
                        product = this_product,
                        ext = 'trimmed_pb',
                        casa = False)

                    # Put together some flags indicating whether these
                    # data have single dish, are a mosaic, are JUST a
                    # mosaic (i.e., have no MS files / original
                    # images), etc.

                    has_imaging = os.path.isdir(imaging_dir + orig_file)
                    has_singledish = (orig_sd_file is not None)
                    is_mosaic = self._kh.is_target_linear_mosaic(this_target)
                    if is_mosaic:                        
                        mosaic_parts = self._kh.get_parts_for_linear_mosaic(this_target)
                    else:
                        mosaic_parts = None

                    # Copy the data from the original location to the
                    # postprocessing directories.

                    if do_stage and config_type == 'interf' and \
                            has_imaging:

                        for fname in [orig_file, pb_file]:
                            
                            infile = imaging_dir + fname
                            outfile = postprocess_dir + fname

                            if not self._quiet:
                                print(this_stub+" ")
                                print(this_stub+"Staging data: ")
                                print(this_stub+"... using ccr.copy_dropdeg")
                                print(this_stub+"... from "+infile)
                                print(this_stub+"... to "+outfile)
                                print(this_stub+" ")
                            
                            if not self._dry_run:
                                ccr.copy_dropdeg(
                                    infile=infile,
                                    outfile=outfile,
                                    overwrite=True,
                                    quiet=self._quiet)

                    # Apply the primary beam correction to the data.

                    if do_pbcorr and config_type == 'interf' and \
                            has_imaging:

                        infile = postprocess_dir + orig_file
                        outfile = postprocess_dir + pbcorr_file
                        pbfile = postprocess_dir + pb_file

                        if not self._quiet:
                            print(this_stub+" ")
                            print(this_stub+"Primary beam correcting:")
                            print(this_stub+"... using ccr.primary_beam_correct")
                            print(this_stub+"... from "+infile)
                            print(this_stub+"... to "+outfile)
                            print(this_stub+"... using "+pbfile)
                            print(this_stub+" ")

                        if not self._dry_run:
                            ccr.primary_beam_correct(
                                infile=infile,
                                outfile=outfile,
                                pbfile=pbfile,
                                overwrite=True,
                                quiet=self._quiet)

                    # Convolve the data to have a round beam.

                    if do_round and config_type == 'interf' and \
                            has_imaging:
                        
                        infile = postprocess_dir + pbcorr_file
                        outfile = postprocess_dir + pbcorr_round_file

                        if not self._quiet:
                            print(this_stub+" ")
                            print(this_stub+"Convolving to have a round beam:")
                            print(this_stub+"... using ccr.convolve_to_round_beam")
                            print(this_stub+"... from "+infile)
                            print(this_stub+"... to "+outfile)
                            print(this_stub+" ")
                        
                        if not self._dry_run:
                            ccr.convolve_to_round_beam(
                                infile=infile,
                                outfile=outfile,
                                overwrite=True,
                                quiet=self._quiet)

                    # Stage the singledish data for feathering

                    if do_singledish and config_type == 'interf' and \
                            has_singledish and has_imaging:

                        template = postprocess_dir + pbcorr_round_file
                        infile = orig_sd_file
                        outfile = postprocess_dir + prepped_sd_file

                        if not self._quiet:
                            print(this_stub+" ")
                            print(this_stub+"Prepping single dish target for feathering:")
                            print(this_stub+"... using cfr.prep_sd_for_feather")
                            print(this_stub+"... original file "+infile)
                            print(this_stub+"... prepped file "+outfile)
                            print(this_stub+"... template "+template)
                            print(this_stub+" ")
                        
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

                    if do_feather and config_type == 'interf' and \
                            has_singledish and has_imaging:

                        interf_file = postprocess_dir + pbcorr_round_file
                        sd_file = postprocess_dir + prepped_sd_file

                        corresponding_feather_config = self._kh.get_feather_config_for_interf_config(
                            interf_config=this_config
                            )
                            
                        feather_file = postprocess_dir + self._kh.get_cube_filename(                            
                            target = this_target,
                            config = corresponding_feather_config,
                            product = this_product,
                            ext = 'pbcorr_round',
                            casa = True,
                            casaext = '.image'
                            )

                        if not self._quiet:
                            print(this_stub+" ")
                            print(this_stub+"Feathering interferometric and single dish data:")
                            print(this_stub+"... using cfr.feather_two_cubes")
                            print(this_stub+"... interferometric data "+interf_file)
                            print(this_stub+"... single dish data "+sd_file)
                            print(this_stub+"... output "+feather_file)
                            print(this_stub+"... feather method: "+self._feather_method)

                        # Feather has a couple of algorithmic choices
                        # associated with it. Run the method that the
                        # user has selected.

                        if self._feather_method == 'apodize':
                                
                            if not self._quiet:
                                print(this_stub+"... apodizing using file "+apod_file)
                                print(this_stub+" ")

                            if not self._dry_run:
                                cfr.feather_two_cubes(
                                    interf_file=interf_file,
                                    sd_file=sd_file,
                                    out_file=feather_file,
                                    apodize=True,
                                    apod_file=apod_file,
                                    apod_cutoff=0.0,
                                    blank=True,
                                    overwrite=True,
                                    quiet=self._quiet)

                        if self._feather_method == 'pbcorr':

                            if not self._quiet:
                                print(this_stub+" ")
                                
                            if not self._dry_run:
                                cfr.feather_two_cubes(
                                    interf_file=interf_file,
                                    sd_file=sd_file,
                                    out_file=feather_file,
                                    apodize=False,
                                    apod_file=None,
                                    apod_cutoff=-1.0,
                                    blank=True,
                                    overwrite=True,
                                    quiet=self._quiet)
                                 
                    # Compress, reducing cube volume.

                    if do_compress:

                        infile = postprocess_dir + pbcorr_round_file
                        outfile = postprocess_dir + pbcorr_trimmed_file

                        if not self._quiet:
                            print(this_stub+" ")
                            print(this_stub+"Reducing cube volume:")
                            print(this_stub+"... using ccr.trim_cube")
                            print(this_stub+"... original file "+infile)
                            print(this_stub+"... output file "+outfile)
                            print(this_stub+" ")

                        if not self._dry_run:
                            ccr.trim_cube(
                                infile=infile,
                                outfile=outfile,
                                overwrite=True,
                                inplace=False,
                                min_pixperbeam=3,
                                quiet=self._quiet)

                        infile_pb = postprocess_dir + pb_file
                        outfile_pb = postprocess_dir + trimmed_pb_file
                        template = postprocess_dir + pbcorr_trimmed_file

                        if not self._quiet:
                            print(this_stub+" ")
                            print(this_stub+"Aligning primary beam image to new astrometry:")
                            print(this_stub+"... using ccr.align_to_target")
                            print(this_stub+"... original file "+infile_pb)
                            print(this_stub+"... output file "+outfile_pb)
                            print(this_stub+"... template "+template)
                            print(this_stub+" ")

                        if not self._dry_run:
                            ccr.align_to_target(
                                infile=infile_pb,
                                outfile=outfile_pb,
                                template=template,
                                interpolation='cubic',
                                overwrite=True,
                                quiet=self._quiet
                                )

                    # Change units from Jy/beam to Kelvin.

                    if do_convert:

                        infile = postprocess_dir + pbcorr_trimmed_file
                        outfile = postprocess_dir + pbcorr_trimmed_k_file
                        
                        if not self._quiet:
                            print(this_stub+" ")
                            print(this_stub+"Converting cube units:")
                            print(this_stub+"... using ccr.convert_jytok")
                            print(this_stub+"... original file "+infile)
                            print(this_stub+"... output file "+outfile)
                            print(this_stub+" ")

                        if not self._dry_run:
                            ccr.convert_jytok(
                                infile=infile,
                                outfile=outfile,
                                overwrite=True,
                                inplace=False,
                                quiet=self._quiet)

                    # Export to FITS and clean up output

                    if do_export:

                        infile = postprocess_dir + pbcorr_trimmed_k_file
                        outfile = postprocess_dir + pbcorr_trimmed_k_fits

                        pb_infile = postprocess_dir + trimmed_pb_file
                        pb_outfile = postprocess_dir + trimmed_pb_fits

                        if not self._quiet:
                            print(this_stub+" ")
                            print(this_stub+"Export to FITS and clean up header: ")
                            print(this_stub+"... using ccr.export_and_cleanup")
                            print(this_stub+"... input cube "+infile)
                            print(this_stub+"... output cube "+outfile)
                            print(this_stub+"... input primary beam "+pb_infile)
                            print(this_stub+"... output primary beam "+pb_outfile)
                            print(this_stub+" ")

                        if not self._dry_run:
                            ccr.export_and_cleanup(
                                infile=infile,
                                outfile=outfile,
                                overwrite=True,    
                                remove_cards=[],
                                add_cards=[],
                                add_history=[],
                                zap_history=True,
                                round_beam=True,
                                roundbeam_tol=0.01,
                                quiet=self._quiet)

                            ccr.export_and_cleanup(
                                infile=pb_infile,
                                outfile=pb_outfile,
                                overwrite=True,    
                                remove_cards=[],
                                add_cards=[],
                                add_history=[],
                                zap_history=True,
                                round_beam=False,
                                roundbeam_tol=0.01,
                                quiet=self._quiet)
                            
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

        self._master_loop(do_singledish=True)

        return()

    def feather(
        self
        ):
        """
        Feather the singledish and interferometer data.
        """

        self._master_loop(do_feather=True)

        return()

    def compress(
        self
        ):
        """
        Reduce cube volume.
        """

        self._master_loop(do_compress=True)

        return()

    def convert(
        self
        ):
        """
        Convert cubes from Jy/beam to K.
        """

        self._master_loop(do_convert=True)

        return()

    def export(
        self
        ):
        """
        Export to FITS and clean up headers.
        """

        self._master_loop(do_export=True)

        return()

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
