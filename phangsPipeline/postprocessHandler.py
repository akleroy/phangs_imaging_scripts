"""
The PHANGS pipeline to handle post-processing of cubes. Works through
a single big class (the PostProcessHandler) that needs to be attached
to a keyHandler. Then it calls the standalone routines.
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

    def _fname_dict(
        self,
        target=None,
        config=None,
        product=None,
        ):
        """
        Make the file name dictionary for all postprocess files given
        some target, config, product configuration.
        """

        fname_dict = {}

        # Original cube and primary beam file
                    
        tag = 'orig'
        orig_file = self._kh.get_cube_filename(
            target = target,
            config = config,
            product = product,
            ext = None,
            casa = True,
            casaext = '.image')
        fname_dict[tag] = orig_file
        
        tag = 'pb'
        pb_file = self._kh.get_cube_filename(
            target = target,
            config = config,
            product = product,
            ext = None,
            casa = True,
            casaext = '.pb')
        fname_dict[tag] = pb_file

        # Original single dish file (note that this comes with a
        # directory)

        tag = 'orig_sd'
        orig_sd_file = self._kh.get_sd_filename(
            target = target,
            product = product)
        fname_dict[tag] = orig_sd_file

        # Primary beam corrected file

        tag = 'pbcorr'
        pbcorr_file = self._kh.get_cube_filename(
            target = target,
            config = config,
            product = product,
            ext = 'pbcorr',
            casa = True,
            casaext = '.image')
        fname_dict[tag] = pbcorr_file

        # Files with round beams

        tag = 'round'
        round_file = self._kh.get_cube_filename(
            target = target,
            config = config,
            product = product,
            ext = 'round',
            casa = True,
            casaext = '.image')
        fname_dict[tag] = round_file

        tag = 'pbcorr_round'
        pbcorr_round_file = self._kh.get_cube_filename(
            target = target,
            config = config,
            product = product,
            ext = 'pbcorr_round',
            casa = True,
            casaext = '.image')
        fname_dict[tag] = pbcorr_round_file

        # Weight file for use in linear mosaicking

        tag = 'weight'
        weight_file = self._kh.get_cube_filename(
            target = target,
            config = config,
            product = product,
            ext = 'weight',
            casa = True,
            casaext = '.image')
        fname_dict[tag] = weight_file

        tag = 'weight_aligned'
        weight_file = self._kh.get_cube_filename(
            target = target,
            config = config,
            product = product,
            ext = 'weight_aligned',
            casa = True,
            casaext = '.image')
        fname_dict[tag] = weight_file

        # Common resolution parts for mosaic

        tag = 'linmos_commonres'
        commonres_file = self._kh.get_cube_filename(
            target = target,
            config = config,
            product = product,
            ext = 'linmos_commonres',
            casa = True,
            casaext = '.image')
        fname_dict[tag] = commonres_file

        # Aligned parts for mosaic

        tag = 'linmos_aligned'
        aligned_file = self._kh.get_cube_filename(
            target = target,
            config = config,
            product = product,
            ext = 'linmos_aligned',
            casa = True,
            casaext = '.image')
        fname_dict[tag] = aligned_file

        # Imported single dish file aligned to the interfometer data

        tag = 'prepped_sd'
        prepped_sd_file = self._kh.get_cube_filename(
            target = target,
            config = config,
            product = product,
            ext = 'singledish',
            casa = True,
            casaext = '.image')
        fname_dict[tag] = prepped_sd_file

        # Singledish weight for use in linear mosaicking

        tag = 'sd_weight'
        sd_weight_file = self._kh.get_cube_filename(
            target = target,
            config = config,
            product = product,
            ext = 'singledish_weight',
            casa = True,
            casaext = '.image')
        fname_dict[tag] = sd_weight_file 

        # Singledish data aliged to a common grid for mosaicking

        tag = 'sd_aligned'
        sd_align_file = self._kh.get_cube_filename(
            target = target,
            config = config,
            product = product,
            ext = 'singledish_aligned',
            casa = True,
            casaext = '.image')
        fname_dict[tag] = sd_align_file 

        # Singledish weight for use in linear mosaicking now on a
        # common astrometric grid

        tag = 'sd_weight_aligned'
        sd_weight_aligned_file = self._kh.get_cube_filename(
            target = target,
            config = config,
            product = product,
            ext = 'singledish_weight_aligned',
            casa = True,
            casaext = '.image')
        fname_dict[tag] = sd_weight_aligned_file 

        # Compressed files with edges trimmed off and smallest
        # reasonable pixel size.

        tag = 'trimmed'
        trimmed_file = self._kh.get_cube_filename(
            target = target,
            config = config,
            product = product,
            ext = 'trimmed',
            casa = True,
            casaext = '.image')
        fname_dict[tag] = trimmed_file
        
        tag = 'pbcorr_trimmed'
        pbcorr_trimmed_file = self._kh.get_cube_filename(
            target = target,
            config = config,
            product = product,
            ext = 'pbcorr_trimmed',
            casa = True,
            casaext = '.image')
        fname_dict[tag] = pbcorr_trimmed_file
        
        tag = 'trimmed_pb'
        trimmed_pb_file = self._kh.get_cube_filename(
            target = target,
            config = config,
            product = product,
            ext = 'trimmed',
            casa = True,
            casaext = '.pb')
        fname_dict[tag] = trimmed_pb_file

        # Files converted to Kelvin, including FITS output files

        tag = 'trimmed_k'
        trimmed_k_file = self._kh.get_cube_filename(
            target = target,
            config = config,
            product = product,
            ext = 'trimmed_k',
            casa = True,
            casaext = '.image')
        fname_dict[tag] = trimmed_k_file

        tag = 'trimmed_k_fits'
        trimmed_k_fits = self._kh.get_cube_filename(
            target = target,
            config = config,
            product = product,
            ext = 'trimmed_k',
            casa = False)
        fname_dict[tag] = trimmed_k_fits
        
        tag = 'pbcorr_trimmed_k'
        pbcorr_trimmed_k_file = self._kh.get_cube_filename(
            target = target,
            config = config,
            product = product,
            ext = 'pbcorr_trimmed_k',
            casa = True,
            casaext = '.image')
        fname_dict[tag] = pbcorr_trimmed_k_file

        tag = 'pbcorr_trimmed_k_fits'
        pbcorr_trimmed_k_fits = self._kh.get_cube_filename(
            target = target,
            config = config,
            product = product,
            ext = 'pbcorr_trimmed_k',
            casa = False)
        fname_dict[tag] = pbcorr_trimmed_k_fits

        tag = 'trimmed_pb_fits'
        trimmed_pb_fits = self._kh.get_cube_filename(
            target = target,
            config = config,
            product = product,
            ext = 'trimmed_pb',
            casa = False)
        fname_dict[tag] = trimmed_pb_fits

        # Return
        
        return(fname_dict)

    def _master_loop(
        self,
        do_stage = False,
        do_pbcorr = False,
        do_round = False,
        do_sd = False,
        do_weight = False,
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

                    fname_dict = self._fname_dict(
                        target=this_target,
                        product=this_product,
                        config=this_config)

                    # Put together some flags indicating whether these
                    # data have single dish, are a mosaic, are JUST a
                    # mosaic (i.e., have no MS files / original
                    # images), etc.

                    has_imaging = os.path.isdir(imaging_dir + fname_dict['orig'])
                    has_sd = (fname_dict['orig_sd'] is not None)
                    is_mosaic = self._kh.is_target_linmos(this_target)
                    if is_mosaic:
                        mosaic_parts = self._kh.get_parts_for_linmos(this_target)
                    else:
                        mosaic_parts = None
                    is_part_of_mosaic = self._kh.is_target_in_mosaic(this_target)
                    
                    # If we are looking at an interferometric
                    # configuration, get the corresponding feather
                    # configuration and also build a file name
                    # dictionary for that.

                    if (config_type == 'interf'):
                        corresponding_feather_config = self._kh.get_feather_config_for_interf_config(
                            interf_config=this_config
                            )
                        feather_fname_dict = self._fname_dict(
                            target=this_target,
                            product=this_product,
                            config=corresponding_feather_config)
                        
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

                        for this_fname in [fname_dict['orig'], fname_dict['pb']]:

                            indir = imaging_dir
                            outdir = postprocess_dir
                            
                            infile = this_fname
                            outfile = this_fname

                            logger.info("Staging "+outfile+" using ccr.copy_dropdeg")

                            if not self._dry_run:
                                ccr.copy_dropdeg(
                                    infile=indir+infile,
                                    outfile=outdir+outfile,
                                    overwrite=True)

                    # Apply the primary beam correction to the data.

                    if do_pbcorr and config_type == 'interf' and \
                            has_imaging:

                        indir = postprocess_dir
                        outdir = postprocess_dir

                        infile = fname_dict['orig']
                        outfile = fname_dict['pbcorr']
                        pbfile = fname_dict['pb']

                        logger.info("Correcting to "+outfile+" using ccr.primary_beam_correct")
                        logger.debug("Correcting from "+infile)
                        logger.debug("Correcting using "+pbfile)

                        if not self._dry_run:
                            ccr.primary_beam_correct(
                                infile=indir+infile,
                                outfile=outdir+outfile,
                                pbfile=indir+pbfile,
                                overwrite=True)

                    # Convolve the data to have a round beam.

                    if do_round and config_type == 'interf' and \
                            has_imaging:

                        indir = postprocess_dir
                        outdir = postprocess_dir
                        
                        infile = fname_dict['pbcorr']
                        outfile = fname_dict['pbcorr_round']

                        logger.info("Convolving to "+outfile+" using ccr.convolve_to_round_beam")
                        logger.debug("Convolving from "+infile)
                        
                        if not self._dry_run:
                            ccr.convolve_to_round_beam(
                                infile=indir+infile,
                                outfile=outdir+outfile,
                                overwrite=True)

                    # Stage the singledish data for feathering

                    if do_sd and config_type == 'interf' and \
                            has_sd and has_imaging:

                        indir = ''
                        outdir = postprocess_dir
                        tempdir = postprocess_dir

                        template = fname_dict['pbcorr_round']
                        infile = fname_dict['orig_sd']
                        outfile = fname_dict['prepped_sd']
                        
                        logger.info("Prepping "+outfile+" for using cfr.prep_sd_for_feather.")
                        logger.debug("Original file "+infile)
                        logger.debug("Using interferometric template "+template)
                        
                        if not self._dry_run:
                            cfr.prep_sd_for_feather(
                                sdfile_in=indir+infile,
                                sdfile_out=outdir+outfile,
                                interf_file=tempdir+template,
                                do_import=True,
                                do_dropdeg=True,
                                do_align=True,
                                do_checkunits=True,                                
                                overwrite=True)

                    # Create a weight file for targets that are part
                    # of a linear mosaic

                    if do_weight and config_type == 'interf' and \
                            is_part_of_mosaic and has_imaging:
                        
                        indir = postprocess_dir
                        outdir = postprocess_dir
                        
                        image_file = fname_dict['pbcorr_round']
                        infile = fname_dict['pb']
                        outfile = fname_dict['weight']

                        logger.info("Making weight file "+outfile+" for using cmr.generate_weight_file.")
                        logger.debug("Based off of primary beam file "+infile)
                        logger.debug("Measuring noise from file "+image_file)
                        
                        if not self._dry_run:
                            cmr.generate_weight_file(
                                image_file = indir+image_file,
                                input_file = indir+infile,
                                input_type = 'pb',
                                outfile = indir + outfile,
                                scale_by_noise = True,
                                overwrite=True)

                    # Make a weight file for single dish targets that
                    # are part of a linear mosaic

                    if do_weight and config_type == 'interf' and \
                            is_part_of_mosaic and has_sd and has_imaging:
                        
                        indir = postprocess_dir
                        outdir = postprocess_dir
                        
                        image_file = fname_dict['prepped_sd']
                        outfile = fname_dict['sd_weight']

                        logger.info("Making weight file "+outfile+" for using cmr.generate_weight_file.")                        
                        logger.debug("Measuring noise from file "+image_file)
                        
                        if not self._dry_run:
                            cmr.generate_weight_file(
                                image_file = indir+image_file,
                                input_value = 1.0,
                                input_type = 'weight',
                                outfile = indir + outfile,
                                scale_by_noise = True,
                                overwrite=True)

                    # Convolve data that are part of a linear mosaic
                    # to a common resolution

                    if do_conv_for_mosaic and config_type == 'interf' and \
                            is_mosaic:

                        indir = postprocess_dir
                        outdir = postprocess_dir

                        # Build the list of input files

                        infile_list = []
                        outfile_list = []

                        for this_part in mosaic_parts:
                            
                            this_part_dict = self._fname_dict(
                                target=this_part,
                                config=this_config,
                                product=this_product,
                                )

                            infile_list.append(indir+this_part_dict['pbcorr_round'])
                            outfile_list.append(outdir+this_part_dict['linmos_commonres'])

                        logger.info("Convolving "+this_target+" for using cmr.common_res_for_mosaic.")
                        logger.debug("Convolving original files "+str(infile_list))
                        logger.debug("Convolving to convolved output "+str(outfile_list))

                        # Allow overrides for the pixel padding (the
                        # number of pixels added to the greatest
                        # common beam for calculating the target
                        # resolution) and the target resolution.

                        pixel_padding = 2.0
                        target_res = None
                        
                        # TBD - check override dict for target
                        # resolution and (I guess?) pixel padding.

                        if not self._dry_run:
                            cmr.common_res_for_mosaic(
                                infile_list = infile_list,
                                outfile_list = outfile_list,
                                do_convolve = True,
                                target_res = target_res,
                                pixel_padding = pixel_padding,
                                overwrite=True,
                                )

                    # Generate a header and align both data and
                    # weights to this new for use in linear
                    # mosaicking. This step is for the interferometric
                    # data.

                    if do_align_for_mosaic and config_type == 'interf' and \
                            is_mosaic:

                        indir = postprocess_dir
                        outdir = postprocess_dir

                        # Build the list of input files

                        infile_list = []
                        outfile_list = []

                        # Get the input and output files for
                        # individual parts. Also include the weights
                        # here.

                        for this_part in mosaic_parts:
                            
                            this_part_dict = self._fname_dict(
                                target=this_part,
                                config=this_config,
                                product=this_product,
                                )

                            infile_list.append(indir+this_part_dict['linmos_commonres'])
                            outfile_list.append(outdir+this_part_dict['linmos_aligned'])

                            infile_list.append(indir+this_part_dict['weight'])
                            outfile_list.append(outdir+this_part_dict['weight_aligned'])

                        logger.info("Aligning "+this_target+" using cmr.common_grid_for_mosaic.")
                        logger.debug("Convolving original files "+str(infile_list))
                        logger.debug("Convolving to convolved output "+str(outfile_list))

                        # TBD implement overrides

                        ra_ctr = None 
                        dec_ctr = None
                        delta_ra = None 
                        delta_dec = None
                        
                        if not self._dry_run:
                            cmr.common_grid_for_mosaic(
                                infile_list = infile_list,
                                outfile_list = outfile_list,
                                ra_ctr = ra_ctr, 
                                dec_ctr = dec_ctr,
                                delta_ra = delta_ra, 
                                delta_dec = delta_dec,
                                allow_big_image = False,
                                too_big_pix=1e4,   
                                asvelocity=True,
                                interpolation='cubic',
                                axes=[-1],
                                overwrite=True,
                                )

                    # Generate a header and align both data and
                    # weights to this new for use in linear
                    # mosaicking. This step is for the single dish
                    # data.

                    if do_align_for_mosaic and config_type == 'interf' and \
                            is_mosaic and has_singleidsh:

                        indir = postprocess_dir
                        outdir = postprocess_dir

                        # Build the list of input files

                        infile_list = []
                        outfile_list = []

                        # Get the input and output files for
                        # individual parts. Also include the weights
                        # here.

                        for this_part in mosaic_parts:
                            
                            this_part_dict = self._fname_dict(
                                target=this_part,
                                config=this_config,
                                product=this_product,
                                )

                            infile_list.append(indir+this_part_dict['prepped_sd'])
                            outfile_list.append(outdir+this_part_dict['sd_aligned'])

                            infile_list.append(indir+this_part_dict['sd_weight'])
                            outfile_list.append(outdir+this_part_dict['sd_weight_aligned'])

                        logger.info("Aligning "+this_target+" using cmr.common_grid_for_mosaic.")
                        logger.debug("Convolving original files "+str(infile_list))
                        logger.debug("Convolving to convolved output "+str(outfile_list))

                        # TBD implement overrides

                        ra_ctr = None 
                        dec_ctr = None
                        delta_ra = None 
                        delta_dec = None
                        
                        if not self._dry_run:
                            cmr.common_grid_for_mosaic(
                                infile_list = infile_list,
                                outfile_list = outfile_list,
                                ra_ctr = ra_ctr, 
                                dec_ctr = dec_ctr,
                                delta_ra = delta_ra, 
                                delta_dec = delta_dec,
                                allow_big_image = False,
                                too_big_pix=1e4,   
                                asvelocity=True,
                                interpolation='cubic',
                                axes=[-1],
                                overwrite=True,
                                )

                    # Execute linear mosaicking for the interferometer data

                    if do_linmos and config_type == 'interf' and \
                            is_mosaic:

                        indir = postprocess_dir
                        outdir = postprocess_dir

                        infile_list = []
                        weightfile_list = []

                        # TBD

                        logger.info("Creating "+outfile+" using cmr.mosaic_aligned_data.")
                        logger.debug("Mosaicking original files "+str(infile_list))
                        logger.debug("Weighting by "+str(weightfile_list))

                        if not self._dry_run:
                            cmr.mosaic_aligned_data(
                                infile_list = infile_list,
                                weightfile_list = weightfile_list,
                                outfile = outfile,
                                overwrite=True)

                    # Execute linear mosaicking for the single dish data

                    if do_linmos and config_type == 'interf' and \
                            is_mosaic:

                        indir = postprocess_dir
                        outdir = postprocess_dir

                        infile_list = []
                        weightfile_list = []

                        # TBD

                        logger.info("Creating "+outfile+" using cmr.mosaic_aligned_data.")
                        logger.debug("Mosaicking original files "+str(infile_list))
                        logger.debug("Weighting by "+str(weightfile_list))

                        if not self._dry_run:
                            cmr.mosaic_aligned_data(
                                infile_list = infile_list,
                                weightfile_list = weightfile_list,
                                outfile = outfile,
                                overwrite=True)
                            
                    # TBD execute linear mosaicking for the feathered
                    # data (need to get the name scheme right here)

                    # Feather the single dish and interferometer data

                    if do_feather and config_type == 'interf' and \
                            has_sd and has_imaging:

                        indir = postprocess_dir

                        interf_file = fname_dict['pbcorr_round']
                        sd_file = fname_dict['prepped_sd']

                        outdir = postprocess_dir
                            
                        outfile = self._kh.get_cube_filename(                            
                            target = this_target,
                            config = corresponding_feather_config,
                            product = this_product,
                            ext = 'pbcorr_round',
                            casa = True,
                            casaext = '.image'
                            )

                        logger.info("Feathering "+outfile+" using cfr.feather_two_cubes.")
                        logger.debug("Feathering interferometric data "+interf_file)
                        logger.debug("Feathering single dish data "+sd_file)
                        logger.debug("Feathering method: "+self._feather_method)

                        # Feather has a couple of algorithmic choices
                        # associated with it. Run the method that the
                        # user has selected.

                        if self._feather_method == 'apodize':
                                
                            logger.debug("Apodizing using file "+apod_file)

                            if not self._dry_run:
                                cfr.feather_two_cubes(
                                    interf_file=indir+interf_file,
                                    sd_file=indir+sd_file,
                                    out_file=outdir+outfile,
                                    do_blank=True,
                                    do_apodize=True,
                                    apod_file=apod_file,
                                    apod_cutoff=0.0,
                                    overwrite=True)

                        if self._feather_method == 'pbcorr':
                                
                            if not self._dry_run:
                                cfr.feather_two_cubes(
                                    interf_file=indir+interf_file,
                                    sd_file=indir+sd_file,
                                    out_file=outdir+outfile,
                                    do_blank=True,
                                    do_apodize=False,
                                    apod_file=None,
                                    apod_cutoff=-1.0,
                                    overwrite=True)

                    # Compress, reducing cube volume.

                    if do_compress:

                        indir = postprocess_dir
                        outdir = postprocess_dir

                        infile = fname_dict['pbcorr_round']
                        outfile = fname_dict['pbcorr_trimmed']

                        logger.info("Producing "+outfile+" using ccr.trim_cube.")
                        logger.debug("Trimming from original file "+infile)

                        if not self._dry_run:
                            ccr.trim_cube(
                                infile=indir+infile,
                                outfile=outdir+outfile,
                                overwrite=True,
                                inplace=False,
                                min_pixperbeam=3)

                        infile_pb = fname_dict['pb']
                        outfile_pb = fname_dict['trimmed_pb']
                        template = fname_dict['pbcorr_trimmed']

                        logger.info("Aligning primary beam image to new astrometry using ccr.align_to_target.")
                        logger.debug("Aligning original file "+infile_pb)
                        logger.debug("Aligning to produce output file "+outfile_pb)
                        logger.debug("Aligning to template "+template)

                        if not self._dry_run:
                            ccr.align_to_target(
                                infile=indir+infile_pb,
                                outfile=outdir+outfile_pb,
                                template=indir+template,
                                interpolation='cubic',
                                overwrite=True,
                                )

                    # Change units from Jy/beam to Kelvin.

                    if do_convert:

                        indir = postprocess_dir
                        outdir = postprocess_dir

                        infile = fname_dict['pbcorr_trimmed']
                        outfile = fname_dict['pbcorr_trimmed_k']
                        
                        logger.info("Creating "+outfile+" using ccr.convert_jytok")
                        logger.debug("Converting from original file "+infile)

                        if not self._dry_run:
                            ccr.convert_jytok(
                                infile=indir+infile,
                                outfile=outdir+outfile,
                                overwrite=True,
                                inplace=False,
                                )

                    # Export to FITS and clean up output

                    if do_export:

                        indir = postprocess_dir
                        outdir = postprocess_dir

                        infile = fname_dict['pbcorr_trimmed_k']
                        outfile = fname_dict['pbcorr_trimmed_k_fits']

                        infile_pb = fname_dict['trimmed_pb']
                        outfile_pb = fname_dict['trimmed_pb_fits']

                        logger.info("Export to "+outfile+" using ccr.export_and_cleanup")
                        logger.debug("Writing from input cube "+infile)
                        logger.debug("Writing from primary beam "+infile_pb)
                        logger.debug("Writing output primary beam "+outfile_pb)

                        if not self._dry_run:
                            ccr.export_and_cleanup(
                                infile=indir+infile,
                                outfile=outdir+outfile,
                                overwrite=True,    
                                remove_cards=[],
                                add_cards=[],
                                add_history=[],
                                zap_history=True,
                                round_beam=True,
                                roundbeam_tol=0.01,
                                )

                            ccr.export_and_cleanup(
                                infile=indir+infile_pb,
                                outfile=outdir+outfile_pb,
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

    def make_weight(
        self
        ):
        """
        Make weights for the interferometric and single dish images
        that are part of a mosaic.
        """

        self._master_loop(do_weight=True)

        return()

    def conv_for_mosaic(
        self
        ):
        """
        Convolve the parts of mosaic to a common resolution.
        """

        self._master_loop(do_conv_for_mosaic=True)

        return()

    def align_for_mosaic(
        self
        ):
        """
        Convolve the parts of mosaic to a common resolution.
        """

        self._master_loop(do_align_for_mosaic=True)

        return()

    def linmos(
        self
        ):
        """
        Convolve the parts of mosaic to a common resolution.
        """

        self._master_loop(do_linmos=True)

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
