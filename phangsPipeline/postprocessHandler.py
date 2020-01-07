"""
Parts of the PHANGS pipeline to handle post-processing cubes.
"""

import os
import glob
import casaCubeRoutines as ccr
import casaMosaicRoutines as cmr
import casaFeatherRoutines as cfr

import logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

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
        dochecks = True
        ):

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
            logger.error("Not a valid feather method: "+method)
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
            logger.error("Cannot build lists without a keyHandler.")
            return(None)

        self._targets_list = self._kh.get_targets(            
            only = self._targets_only,
            skip = self._targets_skip,
            first = self._targets_first,
            last = self._targets_last,
            )

        self._mosaics_list = self._kh.get_linmos_targets(            
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

#region Master loop and master file name routines

    def _master_loop(
        self,
        do_stage = False,
        do_pbcorr = False,
        do_round = False,
        do_sd = False,
        do_conv_for_mosaic = False,
        do_align_for_mosaic = False,
        do_linmos = False,
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
            logger.error("Need a target and interferometer configuration list.")
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

                    orig_sd_file = self._kh.get_sd_filename(
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
                    has_sd = (orig_sd_file is not None)
                    is_mosaic = self._kh.is_target_linmos(this_target)
                    if is_mosaic:                        
                        mosaic_parts = self._kh.get_parts_for_linmos(this_target)
                    else:
                        mosaic_parts = None

                    # Skip out of the loop if we are in a feather
                    # configuration and working with a target that
                    # does not have single dish data and is not a
                    # mosaic.

                    if (config_type == 'feather') and \
                            (not is_mosaic) and (not has_sd):
                        continue

                    # Copy the data from the original location to the
                    # postprocessing directories.

                    if do_stage and config_type == 'interf' and \
                            has_imaging:

                        for fname in [orig_file, pb_file]:
                            
                            infile = imaging_dir + fname
                            outfile = postprocess_dir + fname

                            logger.info("Staging "+fname+" using ccr.copy_dropdeg")

                            if not self._dry_run:
                                ccr.copy_dropdeg(
                                    infile=infile,
                                    outfile=outfile,
                                    overwrite=True)

                    # Apply the primary beam correction to the data.

                    if do_pbcorr and config_type == 'interf' and \
                            has_imaging:

                        infile = postprocess_dir + orig_file
                        outfile = postprocess_dir + pbcorr_file
                        pbfile = postprocess_dir + pb_file

                        logger.info("Correcting to "+pbcorr_file+" using ccr.primary_beam_correct")
                        logger.debug("Correcting from "+orig_file)
                        logger.debug("Correcting using "+pb_file)

                        if not self._dry_run:
                            ccr.primary_beam_correct(
                                infile=infile,
                                outfile=outfile,
                                pbfile=pbfile,
                                overwrite=True)

                    # Convolve the data to have a round beam.

                    if do_round and config_type == 'interf' and \
                            has_imaging:
                        
                        infile = postprocess_dir + pbcorr_file
                        outfile = postprocess_dir + pbcorr_round_file

                        logger.info("Convolving to "+pbcorr_round_file+" using ccr.convolve_to_round_beam")
                        logger.debug("Convolving from "+pbcorr_file)
                        
                        if not self._dry_run:
                            ccr.convolve_to_round_beam(
                                infile=infile,
                                outfile=outfile,
                                overwrite=True)

                    # Stage the singledish data for feathering

                    if do_sd and config_type == 'interf' and \
                            has_sd and has_imaging:

                        template = postprocess_dir + pbcorr_round_file
                        infile = orig_sd_file
                        outfile = postprocess_dir + prepped_sd_file
                        
                        logger.info("Prepping "+prepped_sd_file+" for using cfr.prep_sd_for_feather.")
                        logger.debug("Original file "+orig_sd_file)
                        logger.debug("Using interferometric template "+pbcorr_round_file)
                        
                        if not self._dry_run:
                            cfr.prep_sd_for_feather(
                                sdfile_in=infile,
                                sdfile_out=outfile,
                                interf_file=template,
                                do_import=True,
                                do_dropdeg=True,
                                do_align=True,
                                do_checkunits=True,                                
                                overwrite=True)

                    # Prepare data for linear mosaicking

                    if do_conv_for_mosaic and config_type == 'interf' and \
                            is_mosaic:

                        # Build the list of input files

                        infile_list = []
                        outfile_list = []

                        # TBD - build list of files

                        logger.info("Convolving "+this_target+" for linear mosaicking using cmr.common_res_for_mosaic.")
                        logger.debug("Convolving original files "+str(infile_list))
                        logger.debug("Convolving to convolved output "+str(outfile_list))

                        # Allow overrides for the pixel padding (the
                        # number of pixels added to the greatest
                        # common beam for calculating the target
                        # resolution) and the target resolution.

                        pixel_padding = 2.0
                        target_res = None
                        
                        # TBD - check override dict.

                        if not self._dry_run:
                            cmr.common_res_for_mosaic(
                                infile_list = infile_list,
                                outfile_list = outfile_list,
                                doconvolve = True,
                                target_res = target_res,
                                pixel_padding = pixel_padding,
                                overwrite=True,
                                )

                    # Generate a header and align data to this for the
                    # linear mosaic

                    if do_align_for_mosaic and config_type == 'interf' and \
                            is_mosaic:

                        # Build the list of input files

                        infile_list = []
                        outfile_list = []

                        logger.info("Aligning "+this_target+" for linear mosaicking using cmr._for_mosaic.")
                        logger.debug("Using original files "+str(infile_list))
                        logger.debug("Ending with aligned files "+str(outfile_list))

                        # TBD

                    # Execute linear mosaicking

                    if do_linmos and config_type == 'interf' and \
                            is_mosaic:

                        pass
                            
                    # Feather the single dish and interferometer data

                    if do_feather and config_type == 'interf' and \
                            has_sd and has_imaging:

                        interf_file = postprocess_dir + pbcorr_round_file
                        sd_file = postprocess_dir + prepped_sd_file

                        corresponding_feather_config = self._kh.get_feather_config_for_interf_config(
                            interf_config=this_config
                            )
                            
                        feather_file = self._kh.get_cube_filename(                            
                            target = this_target,
                            config = corresponding_feather_config,
                            product = this_product,
                            ext = 'pbcorr_round',
                            casa = True,
                            casaext = '.image'
                            )

                        outfile = postprocess_dir + feather_file

                        logger.info("Feathering "+feather_file+" using cfr.feather_two_cubes.")
                        logger.debug("Feathering interferometric data "+pbcorr_round_file)
                        logger.debug("Feathering single dish data "+prepped_sd_file)
                        logger.debug("Feathering method: "+self._feather_method)

                        # Feather has a couple of algorithmic choices
                        # associated with it. Run the method that the
                        # user has selected.

                        if self._feather_method == 'apodize':
                                
                            logger.debug("Apodizing using file "+apod_file)

                            if not self._dry_run:
                                cfr.feather_two_cubes(
                                    interf_file=interf_file,
                                    sd_file=sd_file,
                                    out_file=outfile,
                                    do_blank=True,
                                    do_apodize=True,
                                    apod_file=apod_file,
                                    apod_cutoff=0.0,
                                    overwrite=True)

                        if self._feather_method == 'pbcorr':
                                
                            if not self._dry_run:
                                cfr.feather_two_cubes(
                                    interf_file=interf_file,
                                    sd_file=sd_file,
                                    out_file=outfile,
                                    do_blank=True,
                                    do_apodize=False,
                                    apod_file=None,
                                    apod_cutoff=-1.0,
                                    overwrite=True)

                    # Compress, reducing cube volume.

                    if do_compress:

                        infile = postprocess_dir + pbcorr_round_file
                        outfile = postprocess_dir + pbcorr_trimmed_file

                        logger.info("Producing "+pbcorr_trimmed_file+" using ccr.trim_cube.")
                        logger.debug("Trimming from original file "+pbcorr_round_file)

                        if not self._dry_run:
                            ccr.trim_cube(
                                infile=infile,
                                outfile=outfile,
                                overwrite=True,
                                inplace=False,
                                min_pixperbeam=3)

                        infile_pb = postprocess_dir + pb_file
                        outfile_pb = postprocess_dir + trimmed_pb_file
                        template = postprocess_dir + pbcorr_trimmed_file

                        logger.info("Aligning primary beam image to new astrometry using ccr.align_to_target.")
                        logger.debug("Aligning original file "+pb_file)
                        logger.debug("Aligning to produce output file "+trimmed_pb_file)
                        logger.debug("Aligning to template "+pbcorr_trimmed_file)

                        if not self._dry_run:
                            ccr.align_to_target(
                                infile=infile_pb,
                                outfile=outfile_pb,
                                template=template,
                                interpolation='cubic',
                                overwrite=True,
                                )

                    # Change units from Jy/beam to Kelvin.

                    if do_convert:

                        infile = postprocess_dir + pbcorr_trimmed_file
                        outfile = postprocess_dir + pbcorr_trimmed_k_file
                        
                        logger.info("Creating "+pbcorr_trimmed_k_file+" using ccr.convert_jytok")
                        logger.debug("Converting from original file "+pbcorr_trimmed_file)

                        if not self._dry_run:
                            ccr.convert_jytok(
                                infile=infile,
                                outfile=outfile,
                                overwrite=True,
                                inplace=False,
                                )

                    # Export to FITS and clean up output

                    if do_export:

                        infile = postprocess_dir + pbcorr_trimmed_k_file
                        outfile = postprocess_dir + pbcorr_trimmed_k_fits

                        pb_infile = postprocess_dir + trimmed_pb_file
                        pb_outfile = postprocess_dir + trimmed_pb_fits

                        logger.info("Export to "+pbcorr_trimmed_k_fits+" using ccr.export_and_cleanup")
                        logger.debug("Writing from input cube "+pbcorr_trimmed_k_file)
                        logger.debug("Writing from primary beam "+trimmed_pb_file)
                        logger.debug("Writing output primary beam "+trimmed_pb_fits)

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
                                )

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
                                )
                            
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
        
    def stage_sd_data(
        self
        ):
        """
        Copy the single dish data and align it to the releveant
        interferometer data set.
        """

        self._master_loop(do_sd=True)

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
