"""
PostProcessHandler

The PHANGS pipeline to handle post-processing of cubes. Works through
a single big class (the PostProcessHandler). This needs to be attached
to a keyHandler to access the target, product, and configuration
keys.

There should not be any direct calls to CASA in here. Eventually, this
should be able to run without CASA enabled (though it won't be able to
call any of the CASA-specific routines). Right now, just avoid direct
calls to CASA from this class.
"""

import os, sys, re, shutil
import glob
import numpy as np
import handlerTemplate

import logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

casa_enabled = (sys.argv[0].endswith('start_casa.py'))
if casa_enabled:
    import casaCubeRoutines as ccr
    import casaMosaicRoutines as cmr
    import casaFeatherRoutines as cfr

class PostProcessHandler(handlerTemplate.HandlerTemplate):
    """
    Class to handle post-processing of ALMA data. Post-processing here
    begins with the results of imaging and proceeds through reduced,
    science-ready data cubes.
    """

    def __init__(
        self,
        key_handler = None,
        dry_run = False,
        ):
        # Can't use super and keep python2/3 agnostic
        handlerTemplate.HandlerTemplate.__init__(self,key_handler = key_handler, dry_run = dry_run)

#region File name routines

    ###########################################
    # Defined file names for various products #
    ###########################################

    def _fname_dict(
        self,
        target=None,
        config=None,
        product=None,
        extra_ext='',
        ):
        """
        Make the file name dictionary for all postprocess files. This
        will give the a big dictionary of names where one can look up
        each type of file (e.g., primary beam corrected, single dish
        aligned, etc.) given some target, config, and product. This
        routine has a lot of hard-coded knowledge about our
        postprocessing conventions.
        """

        # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
        # Error checking
        # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

        if target is None:
            logger.error("Need a target.")
            return()

        if product is None:
            logger.error("Need a product.")
            return()

        if config is None:
            logger.error("Need a config.")
            return()

        # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
        # Initialize
        # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

        fname_dict = {}

        # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
        # Original files
        # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

        # Original cube
                    
        tag = 'orig'
        orig_file = self._kh.get_cube_filename(
            target = target, config = config, product = product,
            ext = None,
            casa = True,
            casaext = '.image')
        fname_dict[tag] = orig_file
        
        # Original primary beam file

        tag = 'pb'
        pb_file = self._kh.get_cube_filename(
            target = target, config = config, product = product,
            ext = None,
            casa = True,
            casaext = '.pb')
        fname_dict[tag] = pb_file

        # Original single dish file (note that this comes with a
        # directory)

        has_sd = self._kh.has_singledish(target=target, product=product)
        tag = 'orig_sd'
        if has_sd:
            orig_sd_file = self._kh.get_sd_filename(
                target = target, product = product)            
            fname_dict[tag] = orig_sd_file
        else:
            fname_dict[tag] = None

        # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
        # Processed files (apply the extra_ext tag here)
        # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

        # Primary beam corrected file

        tag = 'pbcorr'
        pbcorr_file = self._kh.get_cube_filename(
            target = target, config = config, product = product,
            ext = 'pbcorr'+extra_ext,
            casa = True,
            casaext = '.image')
        fname_dict[tag] = pbcorr_file

        # Files with round beams

        tag = 'round'
        round_file = self._kh.get_cube_filename(
            target = target, config = config, product = product,
            ext = 'round'+extra_ext,
            casa = True,
            casaext = '.image')
        fname_dict[tag] = round_file

        tag = 'pbcorr_round'
        pbcorr_round_file = self._kh.get_cube_filename(
            target = target, config = config, product = product,
            ext = 'pbcorr_round'+extra_ext,
            casa = True,
            casaext = '.image')
        fname_dict[tag] = pbcorr_round_file

        # Weight file for use in linear mosaicking

        tag = 'weight'
        weight_file = self._kh.get_cube_filename(
            target = target, config = config, product = product,
            ext = 'weight'+extra_ext,
            casa = True,
            casaext = '.image')
        fname_dict[tag] = weight_file

        tag = 'weight_aligned'
        weight_file = self._kh.get_cube_filename(
            target = target, config = config, product = product,
            ext = 'weight_aligned'+extra_ext,
            casa = True,
            casaext = '.image')
        fname_dict[tag] = weight_file

        # Common resolution parts for mosaic

        tag = 'linmos_commonres'
        commonres_file = self._kh.get_cube_filename(
            target = target, config = config, product = product,
            ext = 'linmos_commonres'+extra_ext,
            casa = True,
            casaext = '.image')
        fname_dict[tag] = commonres_file

        # Aligned parts for mosaic

        tag = 'linmos_aligned'
        aligned_file = self._kh.get_cube_filename(
            target = target, config = config, product = product,
            ext = 'linmos_aligned'+extra_ext,
            casa = True,
            casaext = '.image')
        fname_dict[tag] = aligned_file

        # Imported single dish file aligned to the interfometer data

        tag = 'prepped_sd'
        prepped_sd_file = self._kh.get_cube_filename(
            target = target, config = config, product = product,
            ext = 'singledish'+extra_ext,
            casa = True,
            casaext = '.image')
        fname_dict[tag] = prepped_sd_file

        # Singledish weight for use in linear mosaicking

        tag = 'sd_weight'
        sd_weight_file = self._kh.get_cube_filename(
            target = target, config = config, product = product,
            ext = 'singledish_weight'+extra_ext,
            casa = True,
            casaext = '.image')
        fname_dict[tag] = sd_weight_file 

        # Singledish data aliged to a common grid for mosaicking

        tag = 'sd_aligned'
        sd_align_file = self._kh.get_cube_filename(
            target = target, config = config, product = product,
            ext = 'singledish_aligned'+extra_ext,
            casa = True,
            casaext = '.image')
        fname_dict[tag] = sd_align_file 

        # Singledish weight for use in linear mosaicking now on a
        # common astrometric grid

        tag = 'sd_weight_aligned'
        sd_weight_aligned_file = self._kh.get_cube_filename(
            target = target, config = config, product = product,
            ext = 'singledish_weight_aligned'+extra_ext,
            casa = True,
            casaext = '.image')
        fname_dict[tag] = sd_weight_aligned_file 

        # Compressed files with edges trimmed off and smallest
        # reasonable pixel size.

        tag = 'trimmed'
        trimmed_file = self._kh.get_cube_filename(
            target = target, config = config, product = product,
            ext = 'trimmed'+extra_ext,
            casa = True,
            casaext = '.image')
        fname_dict[tag] = trimmed_file
        
        tag = 'pbcorr_trimmed'
        pbcorr_trimmed_file = self._kh.get_cube_filename(
            target = target, config = config, product = product,
            ext = 'pbcorr_trimmed'+extra_ext,
            casa = True,
            casaext = '.image')
        fname_dict[tag] = pbcorr_trimmed_file
        
        tag = 'trimmed_pb'
        trimmed_pb_file = self._kh.get_cube_filename(
            target = target, config = config, product = product,
            ext = 'trimmed'+extra_ext,
            casa = True,
            casaext = '.pb')
        fname_dict[tag] = trimmed_pb_file

        # Files converted to Kelvin, including FITS output files

        tag = 'trimmed_k'
        trimmed_k_file = self._kh.get_cube_filename(
            target = target, config = config, product = product,
            ext = 'trimmed_k'+extra_ext,
            casa = True,
            casaext = '.image')
        fname_dict[tag] = trimmed_k_file

        tag = 'trimmed_k_fits'
        trimmed_k_fits = self._kh.get_cube_filename(
            target = target, config = config, product = product,
            ext = 'trimmed_k'+extra_ext,
            casa = False)
        fname_dict[tag] = trimmed_k_fits
        
        tag = 'pbcorr_trimmed_k'
        pbcorr_trimmed_k_file = self._kh.get_cube_filename(
            target = target, config = config, product = product,
            ext = 'pbcorr_trimmed_k'+extra_ext,
            casa = True,
            casaext = '.image')
        fname_dict[tag] = pbcorr_trimmed_k_file

        tag = 'pbcorr_trimmed_k_fits'
        pbcorr_trimmed_k_fits = self._kh.get_cube_filename(
            target = target, config = config, product = product,
            ext = 'pbcorr_trimmed_k'+extra_ext,
            casa = False)
        fname_dict[tag] = pbcorr_trimmed_k_fits

        tag = 'trimmed_pb_fits'
        trimmed_pb_fits = self._kh.get_cube_filename(
            target = target, config = config, product = product,
            ext = 'trimmed_pb'+extra_ext,
            casa = False)
        fname_dict[tag] = trimmed_pb_fits

        # Return
        
        return(fname_dict)

#endregion

#region "Tasks" : Individual postprocessing steps

    def task_stage_interf_data(
        self,
        target = None,
        product = None,
        config = None,
        extra_ext_in = '',
        extra_ext_out = '',
        check_files = True,
        ):
        """
        For one target, product, config combination copy the
        interferometric cube and primary beam file to the working
        postprocessing directory.
        """

        # Generate file names

        indir = self._kh.get_imaging_dir_for_target(target)
        outdir = self._kh.get_postprocess_dir_for_target(target)
        fname_dict_in = self._fname_dict(
            target=target, config=config, product=product, extra_ext=extra_ext_in)
        fname_dict_out = self._fname_dict(
            target=target, config=config, product=product, extra_ext=extra_ext_out)
                
        # Copy the primary beam and the interferometric imaging
        
        for this_tag in ['orig', 'pb']:
            
            infile = fname_dict_in[this_tag]
            outfile = fname_dict_out[this_tag]
        
            # Check input file existence
            if check_files:
                if not (os.path.isdir(indir+infile)):
                    logger.warning("Missing "+indir+infile)
                    continue
    
            logger.info("")
            logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%&%")
            logger.info("Staging data for:")
            logger.info(str(target)+" , "+str(product)+" , "+str(config))
            logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%&%")
            logger.info("")
            logger.info("Using ccr.copy_dropdeg.")
            logger.info("Staging "+outfile)
            
            if (not self._dry_run) and casa_enabled:
                ccr.copy_dropdeg(
                    infile=indir+infile,
                    outfile=outdir+outfile,
                    overwrite=True)

        return()

    def task_pbcorr(
        self,
        target = None,
        product = None,
        config = None,
        in_tag = 'orig',
        out_tag = 'pbcorr',
        extra_ext_in = '',
        extra_ext_out = '',
        check_files = True,
        ):
        """
        For one target, product, config combination primary beam
        correct the interferometer data.
        """

        # Generate file names

        indir = self._kh.get_postprocess_dir_for_target(target)
        outdir = self._kh.get_postprocess_dir_for_target(target)
        fname_dict_in = self._fname_dict(
            target=target, config=config, product=product, extra_ext=extra_ext_in)
        fname_dict_out = self._fname_dict(
            target=target, config=config, product=product, extra_ext=extra_ext_out)

        infile = fname_dict_in[in_tag]
        outfile = fname_dict_out[out_tag]
        pbfile = fname_dict_in['pb']

        # Check input file existence
         
        if check_files:
            if not (os.path.isdir(indir+infile)):
                logger.warning("Missing "+indir+infile)
                return()
            if not (os.path.isdir(indir+pbfile)):
                logger.warning("Missing "+indir+pbfile)
                return()

        # Apply the primary beam correction to the data.
        
        logger.info("")
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%&%")
        logger.info("Primary beam correction for:")
        logger.info(str(target)+" , "+str(product)+" , "+str(config))
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%&%")
        logger.info("")
        
        logger.info("Using ccr.primary_beam_correct")
        logger.info("Correcting to "+outfile)
        logger.info("Correcting from "+infile)
        logger.info("Correcting using "+pbfile)
        
        if (not self._dry_run) and casa_enabled:
            ccr.primary_beam_correct(
                infile=indir+infile,
                outfile=outdir+outfile,
                pbfile=indir+pbfile,
                overwrite=True)

        return()

    def task_round_beam(
        self,
        target = None,
        product = None,
        config = None,
        in_tag = 'pbcorr',
        out_tag = 'pbcorr_round',
        extra_ext_in = '',
        extra_ext_out = '',
        force_beam_as = None,
        check_files = True,
        ):
        """
        For one target, product, config combination, convolve the cube
        to have a round beam. Note that via the force_beam_as keyword
        this task can also be used to convolve data to a fixed (round)
        angular resolution.
        """

        # Generate file names

        indir = self._kh.get_postprocess_dir_for_target(target)
        outdir = self._kh.get_postprocess_dir_for_target(target)
        fname_dict_in = self._fname_dict(
            target=target, config=config, product=product, extra_ext=extra_ext_in)
        fname_dict_out = self._fname_dict(
            target=target, config=config, product=product, extra_ext=extra_ext_out)
        
        infile = fname_dict_in[in_tag]
        outfile = fname_dict_out[out_tag]

        # Check input file existence        

        if check_files:
            if not (os.path.isdir(indir+infile)):
                logger.warning("Missing "+infile)
                return()

        # Convolve the data to have a round beam.
        
        logger.info("")
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%")
        logger.info("Convolving to a round beam for:")
        logger.info(str(target)+" , "+str(product)+" , "+str(config))
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%")
        logger.info("")
        
        logger.info("Using ccr.convolve_to_round_beam")
        logger.info("Convolving from "+infile)
        logger.info("Convolving to "+outfile)
        if force_beam_as is not None:
            logger.info("Forcing beam to "+str(force_beam_as))
        
        if (not self._dry_run) and casa_enabled:
            ccr.convolve_to_round_beam(
                infile=indir+infile,
                outfile=outdir+outfile,
                force_beam=force_beam_as,
                overwrite=True)

        return()

    def task_stage_singledish(
        self,
        target = None,
        product = None,
        config = None,
        template_tag = 'pbcorr_round',
        out_tag = 'prepped_sd',
        extra_ext_in = '',
        extra_ext_out = '',
        check_files = True,
        ):
        """
        For one target, product, config combination, copy the single
        dish data and align it to the interferometric grid.
        """

        # Generate file names

        indir = ''
        outdir = self._kh.get_postprocess_dir_for_target(target)
        tempdir = self._kh.get_postprocess_dir_for_target(target)
        fname_dict_in = self._fname_dict(
            target=target, config=config, product=product, extra_ext=extra_ext_in)
        fname_dict_out = self._fname_dict(
            target=target, config=config, product=product, extra_ext=extra_ext_out)

        template = fname_dict_in[template_tag]
        infile = fname_dict_in['orig_sd']
        outfile = fname_dict_out[out_tag]

        # Check input file existence        
        
        if check_files:
            if (not (os.path.isdir(indir+infile))) and \
                    (not (os.path.isfile(indir+infile))):
                logger.warning("Missing "+infile)
                return()
            if not (os.path.isdir(tempdir+template)):
                logger.warning("Missing "+tempdir+template)
                return()

        # Stage the singledish data for feathering

        logger.info("")
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%")
        logger.info("Preparing single dish data for:")
        logger.info(str(target)+" , "+str(product)+" , "+str(config))
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%")
        logger.info("")
        
        logger.info("Using cfr.prep_sd_for_feather.")
        logger.info("Prepping "+outfile)
        logger.info("Original file "+infile)
        logger.info("Using interferometric template "+template)
        
        if (not self._dry_run) and casa_enabled:
            cfr.prep_sd_for_feather(
                sdfile_in=indir+infile,
                sdfile_out=outdir+outfile,
                interf_file=tempdir+template,
                do_import=True,
                do_dropdeg=True,
                do_align=True,
                do_checkunits=True,                                
                overwrite=True)

        return()

    def task_make_interf_weight(
        self,
        target = None,
        product = None,
        config = None,
        image_tag = 'pbcorr_round',
        in_tag = 'pb',
        input_type = 'pb',
        scale_by_noise = True,
        out_tag = 'weight',
        extra_ext_in = '',
        extra_ext_out = '',
        check_files = True,
        ):
        """
        For one target, product, config combination, make a 'weight'
        image for use in linearly mosaicking the cube with other,
        overlapping cubes. This task targets interferometric dish
        data.
        """

        # Generate file names

        indir = self._kh.get_postprocess_dir_for_target(target)
        outdir = self._kh.get_postprocess_dir_for_target(target)
        fname_dict_in = self._fname_dict(
            target=target, config=config, product=product, extra_ext=extra_ext_in)
        fname_dict_out = self._fname_dict(
            target=target, config=config, product=product, extra_ext=extra_ext_out)
        
        image_file = fname_dict_in[image_tag]
        infile = fname_dict_in[in_tag]
        outfile = fname_dict_out[out_tag]        

        # Check input file existence        
        
        if check_files:
            if not (os.path.isdir(indir+infile)):
                logger.warning("Missing "+infile)
                return()
            if not (os.path.isdir(indir+image_file)):
                logger.warning("Missing "+image_file)
                return()

        # Create a weight image for use linear mosaicking targets that
        # are part of a linear mosaic

        logger.info("")
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&")
        logger.info("Making weight file for:")
        logger.info(str(target)+" , "+str(product)+" , "+str(config))
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&")
        logger.info("")
        
        logger.info("Using cmr.generate_weight_file.")
        logger.info("Making weight file "+outfile)
        logger.info("Based off of primary beam file "+infile)
        logger.info("Measuring noise from file "+image_file)
                        
        if (not self._dry_run) and casa_enabled:
            cmr.generate_weight_file(
                image_file = indir+image_file,
                input_file = indir+infile,
                input_type = input_type,
                outfile = indir + outfile,
                scale_by_noise = scale_by_noise,
                overwrite=True)

        return()

    def task_make_singledish_weight(
        self,
        target = None,
        product = None,
        config = None,
        image_tag = 'prepped_sd',
        out_tag = 'sd_weight',
        extra_ext_in = '',
        extra_ext_out = '',
        check_files = True,
        ):
        """
        For one target, product, config combination, make a 'weight'
        image for use in linearly mosaicking the cube with other,
        overlapping cubes. This task targets single dish data.
        """

        # Generate file names

        indir = self._kh.get_postprocess_dir_for_target(target)
        outdir = self._kh.get_postprocess_dir_for_target(target)
        fname_dict_in = self._fname_dict(
            target=target, config=config, product=product, extra_ext=extra_ext_in)
        fname_dict_out = self._fname_dict(
            target=target, config=config, product=product, extra_ext=extra_ext_out)
                        
        image_file = fname_dict_in[image_tag]
        outfile = fname_dict_out[out_tag]

        # Check input file existence        
    
        if check_files:
            if not (os.path.isdir(indir+image_file)):
                logger.warning("Missing "+image_file)
                return()

        # Make a weight file for single dish targets that
        # are part of a linear mosaic

        logger.info("")
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&")
        logger.info("Making single dish weight file for:")
        logger.info(str(target)+" , "+str(product)+" , "+str(config))
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&")
        logger.info("")
        
        logger.info("Using cmr.generate_weight_file.")
        logger.info("Making weight file "+outfile)
        logger.info("Measuring noise from file "+image_file)
            
        if (not self._dry_run) and casa_enabled:
            cmr.generate_weight_file(
                image_file = indir+image_file,
                input_value = 1.0,
                input_type = 'weight',
                outfile = indir + outfile,
                scale_by_noise = True,
                overwrite=True)
                
        return()

    def task_feather(
        self,
        target = None,
        product = None,
        config = None,
        interf_tag = 'pbcorr_round',
        sd_tag = 'prepped_sd',
        out_tag = 'pbcorr_round',
        extra_ext_in = '',
        extra_ext_out = '',
        apodize = False,
        apod_ext = 'pb',
        copy_weights = True,
        check_files = True,       
        ):
        """
        For one target, product, config combination, feather together
        a single dish and interferometric data set. Note that
        apodization is exposed as an option. Also note that the
        configuration of the input and output will differ (an
        interferometric configuration comes in, a feather
        configuration comes out). Optionally, propagate the weights
        from the interferometric side to become the weights for the
        new feathered data.
        """

        # Generate file names

        indir = self._kh.get_postprocess_dir_for_target(target)
        outdir = self._kh.get_postprocess_dir_for_target(target)
        fname_dict_in = self._fname_dict(
            target=target, config=config, product=product, extra_ext=extra_ext_in)

        # Note that feather changes the config

        feather_config = self._kh.get_feather_config_for_interf_config(
            interf_config=config)

        fname_dict_out = self._fname_dict(
            target=target, config=feather_config, product=product, 
            extra_ext=extra_ext_out)
                        
        interf_file = fname_dict_in[interf_tag]
        sd_file = fname_dict_in[sd_tag]
        outfile = fname_dict_out[out_tag]

        # Error checking

        # Check input file existence        
    
        # Feather the single dish and interferometer data
                
        logger.info("")
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%")
        logger.info("Feathering interferometer and single dish data for:")
        logger.info(str(target)+" , "+str(product)+" , "+str(config))
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%")
        logger.info("")
        
        logger.info("Using cfr.feather_two_cubes.")
        logger.info("Feathering "+outfile)
        logger.info("Feathering interferometric data "+interf_file)
        logger.info("Feathering single dish data "+sd_file)

        # Feather has a couple of algorithmic choices
        # associated with it. Run the method that the
        # user has selected.
        
        if apodize:

            apod_file = fname_dict_in[apod_ext]

            logger.info("Apodizing using file "+apod_file)

            if (not self._dry_run) and casa_enabled:
                cfr.feather_two_cubes(
                    interf_file=indir+interf_file,
                    sd_file=indir+sd_file,
                    out_file=outdir+outfile,
                    do_blank=True,
                    do_apodize=True,
                    apod_file=indir+apod_file,
                    apod_cutoff=0.0,
                    overwrite=True)
                
        else:
            
            if (not self._dry_run) and casa_enabled:                                
                cfr.feather_two_cubes(
                    interf_file=indir+interf_file,
                    sd_file=indir+sd_file,
                    out_file=outdir+outfile,
                    do_blank=True,
                    do_apodize=False,
                    apod_file=None,
                    apod_cutoff=-1.0,
                    overwrite=True)

        if copy_weights:
                
            interf_weight_exists = False
            interf_weight_file = fname_dict_in['weight']
            if os.path.isdir(indir+interf_weight_file):
                interf_weight_exists = True
            else:
                logger.info("Interferometric weight file not found "+interf_weight_file)

            if interf_weight_exists:
                logger.info("")
                logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%")
                logger.info("Copying weights for:")
                logger.info(str(target)+" , "+str(product)+" , "+str(config))
                logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%")
                logger.info("")
                
                out_weight_file=fname_dict_out['weight']

                logger.info("Copying from "+interf_weight_file)
                logger.info("Copying to "+out_weight_file)
                if (not self._dry_run) and casa_enabled:
                    ccr.copy_dropdeg(infile=indir+interf_weight_file, 
                                     outfile=outdir+out_weight_file, 
                                     overwrite=True)
        return()

    def task_compress(
        self,
        target = None,
        product = None,
        config = None,
        in_tag = 'pbcorr_round',
        out_tag = 'pbcorr_trimmed',
        do_pb_too = True,
        in_pb_tag = 'pb',
        out_pb_tag = 'pb_trimmed',
        extra_ext_in = '',
        extra_ext_out = '',
        check_files = True,
        ):
        """
        For one target, product, config combination, compress the cube
        to the smallest reasonable volume. Also align the primary beam
        file out onto this grid.
        """

        # Generate file names

        indir = self._kh.get_postprocess_dir_for_target(target)
        outdir = self._kh.get_postprocess_dir_for_target(target)
        fname_dict_in = self._fname_dict(
            target=target, config=config, product=product, extra_ext=extra_ext_in)
        fname_dict_out = self._fname_dict(
            target=target, config=config, product=product, extra_ext=extra_ext_out)

        infile = fname_dict_in['pbcorr_round']
        outfile = fname_dict_out['pbcorr_trimmed']

        infile_pb = fname_dict_in['pb']
        outfile_pb = fname_dict_out['trimmed_pb']

        # Check input file existence        

        if check_files:
            if not (os.path.isdir(indir+infile)):
                logger.warning("Missing "+infile)
                return()

        # Compress, reducing cube volume.

        logger.info("")
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%&%")
        logger.info("Trimming cube for:")
        logger.info(str(target)+" , "+str(product)+" , "+str(config))
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%&%")
        logger.info("")

        logger.info("Producing "+outfile+" using ccr.trim_cube.")
        logger.info("Trimming from original file "+infile)
        
        if (not self._dry_run) and casa_enabled:
            ccr.trim_cube(
                infile=indir+infile,
                outfile=outdir+outfile,
                overwrite=True,
                inplace=False,
                min_pixperbeam=3)

        if do_pb_too is False:
            return()

        if check_files:
            if not (os.path.isdir(indir+infile_pb)):
                logger.warning("Missing "+infile_pb)
                return()
            
        template = fname_dict_out['pbcorr_trimmed']

        if check_files:
            if not (os.path.isdir(outdir+template)):
                logger.warning("Missing "+template)
                return()

        logger.info("Aligning primary beam image to new astrometry")
        logger.info("Using ccr.align_to_target.")
        logger.info("Aligning original file "+infile_pb)
        logger.info("Aligning to produce output file "+outfile_pb)
        logger.info("Aligning to template "+template)

        if (not self._dry_run) and casa_enabled:
            ccr.align_to_target(
                infile=indir+infile_pb,
                outfile=outdir+outfile_pb,
                template=outdir+template,
                interpolation='cubic',
                overwrite=True,
                )

        return()

    def task_convert_units(
        self,
        target = None,
        product = None,
        config = None,
        in_tag = 'pbcorr_trimmed',
        out_tag = 'pbcorr_trimmed_k',
        extra_ext_in = '',
        extra_ext_out = '',
        check_files = True,
        ):
        """
        For one target, config, product combination convert the units
        from Jy/beam to Kelvin.
        """

        # Generate file names

        indir = self._kh.get_postprocess_dir_for_target(target)
        outdir = self._kh.get_postprocess_dir_for_target(target)
        fname_dict_in = self._fname_dict(
            target=target, config=config, product=product, extra_ext=extra_ext_in)
        fname_dict_out = self._fname_dict(
            target=target, config=config, product=product, extra_ext=extra_ext_out)
            
        infile = fname_dict_in[in_tag]
        outfile = fname_dict_out[out_tag]

        # Check input file existence        

        if check_files:
            if not (os.path.isdir(indir+infile)):
                logger.warning("Missing "+infile)
                return()

        # Change units from Jy/beam to Kelvin.
                        
        logger.info("")
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%&%")
        logger.info("Converting units for:")
        logger.info(str(target)+" , "+str(product)+" , "+str(config))
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%&%")
        logger.info("")
        
        logger.info("Using ccr.convert_jytok")
        logger.info("Creating "+outfile)
        logger.info("Converting from original file "+infile)
        
        if (not self._dry_run) and casa_enabled:
            ccr.convert_jytok(
                infile=indir+infile,
                outfile=outdir+outfile,
                overwrite=True,
                inplace=False,
                )

        return()

    def task_export_to_fits(
        self,
        target = None,
        product = None,
        config = None,
        in_tag = 'pbcorr_trimmed_k',
        out_tag = 'pbcorr_trimmed_k_fits',
        do_pb_too = True,
        in_pb_tag = 'trimmed_pb',
        out_pb_tag = 'trimmed_pb_fits',
        extra_ext_in = '',
        extra_ext_out = '',
        check_files = True,
        ):
        """
        For one target, config, product combination export to
        FITS. Optionally also export the primary beam files.
        """

        # Generate file names

        indir = self._kh.get_postprocess_dir_for_target(target)
        outdir = self._kh.get_postprocess_dir_for_target(target)
        fname_dict_in = self._fname_dict(
            target=target, config=config, product=product, extra_ext=extra_ext_in)
        fname_dict_out = self._fname_dict(
            target=target, config=config, product=product, extra_ext=extra_ext_out)
        
        infile = fname_dict_in[in_tag]
        outfile = fname_dict_out[out_tag]
        
        # Check input file existence        

        if check_files:
            if not (os.path.isdir(indir+infile)):
                logger.warning("Missing "+infile)
                return()

        # Export to FITS and clean up output
        
        logger.info("")
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%")
        logger.info("Exporting data to FITS and cleaning up cubes for:")
        logger.info(str(target)+" , "+str(product)+" , "+str(config))
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%")
        logger.info("")
        
        logger.info("Using ccr.export_and_cleanup.")
        logger.info("Export to "+outfile)
        logger.info("Writing from input cube "+infile)

        if (not self._dry_run) and casa_enabled:
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

        if do_pb_too is False:
            return()

        # Check input file existence        

        infile_pb = fname_dict_in[in_pb_tag]
        outfile_pb = fname_dict_out[out_pb_tag]

        if check_files:
            if not (os.path.isdir(indir+infile_pb)):
                logger.warning("Missing "+infile_pb)
                return()

        logger.info("Writing from primary beam "+infile_pb)
        logger.info("Writing output primary beam "+outfile_pb)
        
        if (not self._dry_run) and casa_enabled:
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

    def task_convolve_parts_for_mosaic(
        self,
        target = None,
        product = None,
        config = None,
        in_tag = 'pbcorr_round',
        out_tag = 'linmos_commonres',
        extra_ext_in = '',
        extra_ext_out = '',
        check_files = True,
        ):
        """
        For one target, config, product combination that is a linear
        mosaic, convolve all of the parts of the mosaic to share a
        common angular resolution, appropriate for gridding together
        into a single image.
        """

        # Generate file names

        indir = self._kh.get_postprocess_dir_for_target(target)
        outdir = self._kh.get_postprocess_dir_for_target(target)

        mosaic_parts = self._kh.get_parts_for_linmos(target)

        infile_list = []
        outfile_list = []
        
        for this_part in mosaic_parts:
                            
            this_part_dict_in = self._fname_dict(
                target=this_part, config=config, product=product,
                extra_ext=extra_ext_in,
                )

            this_part_dict_out = self._fname_dict(
                target=this_part, config=config, product=product,
                extra_ext=extra_ext_out,
                )

            infile_list.append(indir+this_part_dict_in[in_tag])
            outfile_list.append(outdir+this_part_dict_out[out_tag])
            
        logger.info("")
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%")
        logger.info("Convolving for mosaic for:")
        logger.info(str(target)+" , "+str(product)+" , "+str(config))
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%")
        logger.info("")

        logger.info("Using cmr.common_res_for_mosaic.")
        logger.info("Convolving "+target)
        logger.info("Convolving original files "+str(infile_list))
        logger.info("Convolving to convolved output "+str(outfile_list))
        
        # Allow overrides for the pixel padding (the
        # number of pixels added to the greatest
        # common beam for calculating the target
        # resolution) and the target resolution.
            
        pixel_padding = 2.0
        target_res = None
                        
        # TBD - check override dict for target
        # resolution and (maybe?) pixel padding.

        if (not self._dry_run) and casa_enabled:
            cmr.common_res_for_mosaic(
                infile_list = infile_list,
                outfile_list = outfile_list,
                do_convolve = True,
                target_res = target_res,
                pixel_padding = pixel_padding,
                overwrite=True,
                )

        return()


    def task_align_for_mosaic(
        self,
        target = None,
        product = None,
        config = None,
        in_tags = ['linmos_commonres', 'weight', 'prepped_sd', 'sd_weight'],
        out_tags = ['linmos_aligned', 'weight_aligned', 'sd_aligned', 'sd_weight_aligned'],
        extra_ext_in = '',
        extra_ext_out = '',
        check_files = True,
        ):
        """
        For one target, config, product combination that is a linear
        mosaic, align all parts of the mosaic to a common astrometric
        grid for combination into a single image.
        """
        
        # Map the input and output tags to one another in a dictionary

        if (type(in_tags) != type([])) or type(out_tags) != type([]):
            logger.error("Input and output tag lists must be lists.")
            return(None)
            
        if len(in_tags) != len(out_tags):
            logger.error("Mismatch in input and output list tag list.")
            return(None)

        out_tag_dict = {}
        for ii in range(len(in_tags)):
            out_tag_dict[in_tags[ii]] = out_tags[ii]

        # Generate file names

        indir = self._kh.get_postprocess_dir_for_target(target)
        outdir = self._kh.get_postprocess_dir_for_target(target)

        mosaic_parts = self._kh.get_parts_for_linmos(target)

        infile_list = []
        outfile_list = []
        
        for this_part in mosaic_parts:
            
            this_part_dict_in = self._fname_dict(
                target=this_part, config=config, product=product,
                extra_ext=extra_ext_in,
                )
            
            this_part_dict_out = self._fname_dict(
                target=this_part, config=config, product=product,
                extra_ext=extra_ext_out,
                )
            
            for this_tag_in in in_tags:
                
                this_tag_out = out_tag_dict[this_tag_in]
                infile_list.append(indir+this_part_dict_in[this_tag_in])
                outfile_list.append(outdir+this_part_dict_out[this_tag_out])

        logger.info("")
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%")
        logger.info("Aligning for mosaic for:")
        logger.info(str(target)+" , "+str(product)+" , "+str(config))
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%")
        logger.info("")
            
        logger.info("Using cmr.common_grid_for_mosaic.")
        logger.info("Aligning "+target)
        logger.info("Convolving original files "+str(infile_list))
        logger.info("Convolving to convolved output "+str(outfile_list))
        
        # TBD implement overrides
            
        ra_ctr = None 
        dec_ctr = None
        delta_ra = None 
        delta_dec = None
        
        if (not self._dry_run) and casa_enabled:
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

        return()

    def task_linear_mosaic(
        self,
        target = None,
        product = None,
        config = None,
        image_tag = 'linmos_aligned', # 'sd_aligned'
        weight_tag = 'weight_aligned', # 'sd_weight_aligned'
        out_tag = 'pbcorr_round', # 'prepped_sd'
        extra_ext_in = '',
        extra_ext_out = '',
        check_files = True,
        ):
        """
        For one target, config, product combination that is a linear
        mosaic and has already been aligned and convolved, execute the
        linear mosaic. Needs to be run separately for single dish and
        interferometer data.
        """

        # Set input and output directories and define output file

        indir = self._kh.get_postprocess_dir_for_target(target)
        outdir = self._kh.get_postprocess_dir_for_target(target)

        fname_dict_out = self._fname_dict(
            target=target, config=config, product=product, 
            extra_ext=extra_ext_out)

        outfile = fname_dict_out[out_tag]        

        mosaic_parts = self._kh.get_parts_for_linmos(target)

        infile_list = []
        weightfile_list = []

        # Get the input and weight files for  individual parts.
        
        for this_part in mosaic_parts:

            this_part_dict_in = self._fname_dict(
                target=this_part, config=config, product=product,
                extra_ext=extra_ext_in,
                )
            
            infile_list.append(indir+this_part_dict_in[image_tag])
            weightfile_list.append(indir+this_part_dict_in[weight_tag])

        logger.info("")
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%")
        logger.info("Executing linear mosaic for:")
        logger.info(str(target)+" , "+str(product)+" , "+str(config))
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%")
        logger.info("")
            
        logger.info("Using cmr.mosaic_aligned_data.")
        logger.info("Creating "+outfile)
        logger.info("Mosaicking original files "+str(infile_list))
        logger.info("Weighting by "+str(weightfile_list))

        if not self._dry_run:
            cmr.mosaic_aligned_data(
                infile_list = infile_list,
                weightfile_list = weightfile_list,
                outfile = outdir+outfile,
                overwrite=True)
                
        return()

#endregion

#region Recipes execute a set of linked tasks for one data set.
    
    def recipe_prep_one_target(
        self,
        target = None,
        product = None,
        config = None,
        check_files = True,
        ):
        """
        Recipe that takes data from imaging through all steps that
        come before feathering and/or mosaicking. This means copying
        the data, primary beam correction, convolution to a round
        beam, importing and aligning the single dish data, and making
        weight files for targets that are part of linear mosaicks.

        The recipe assumes a lot of file name conventions so just
        needs the target/product/config defined.
        """

        # Work out file names and note whether the target is part of a
        # mosaic, has single dish data, etc.

        fname_dict = self._fname_dict(
            target=target, product=product, config=config)

        imaging_dir = self._kh.get_imaging_dir_for_target(target)
        has_imaging = os.path.isdir(imaging_dir + fname_dict['orig'])
        has_singledish = self._kh.has_singledish(target=target, product=product)        
        is_part_of_mosaic = self._kh.is_target_in_mosaic(target)

        if not has_imaging:
            logger.warning("No imaging for "+fname_dict['orig']+". Returning.")
            return()

        # Call tasks

        self.task_stage_interf_data(
            target=target, config=config, product=product,
            check_files=check_files
            )

        self.task_pbcorr(
            target=target, config=config, product=product,
            check_files=check_files
            )

        self.task_round_beam(
            target=target, config=config, product=product,
            check_files=check_files
            )

        if has_singledish:
            self.task_stage_singledish(
                target=target, config=config, product=product,
                check_files=check_files
                )

        if is_part_of_mosaic:
            self.task_make_interf_weight(
                target=target, config=config, product=product,
                check_files=check_files, scale_by_noise=True,
                )

        if is_part_of_mosaic and has_singledish:
            self.task_make_singledish_weight(
                target=target, config=config, product=product,
                check_files=check_files,
                )

        return()

    def recipe_mosaic_one_target(
        self,
        target = None,
        product = None,
        config = None,
        check_files = True,
        extra_ext_in = '',
        extra_ext_out = '',
        ):
        """
        Linearly mosaic a single target, performing the convolution,
        alignment, and mosaicking steps.
        """

        # Check that the target is a mosaic

        is_mosaic = self._kh.is_target_linmos(target)

        if not is_mosaic:
            logger.warning("Not a mosaic, returning.")
            return()

        mosaic_parts = self._kh.get_parts_for_linmos(target)
            
        # Check if the individual parts have single dish data. If they
        # do, flip the single dish flag to true.

        parts_have_singledish = False

        for this_part in mosaic_parts:
            
            this_part_has_sd = self._kh.has_singledish(target=this_part, product=product) 

            if this_part_has_sd:
                parts_have_singledish = True
            
        # Check if this is a feather configuration. If so, then flip
        # the single dish flag to false. This overrides the presence
        # of data - we don't treat the singledish for feathered data.

        if config in self.get_feather_configs():

            parts_have_singledish = False
    
        self.task_convolve_parts_for_mosaic(
            target = target,
            product = product,
            config = config,
            in_tag = 'pbcorr_round',
            out_tag = 'linmos_commonres',
            extra_ext_in = extra_ext_in,
            extra_ext_out = extra_ext_in,
            check_files = check_files,
            )

        in_tag_list = ['linmos_commonres', 'weight']
        out_tag_list = ['linmos_aligned', 'weight_aligned']

        if parts_have_singledish:
            in_tag_list.append('prepped_sd')
            in_tag_list.append('sd_weight')
            out_tag_list.append('sd_aligned')
            out_tag_list.append('sd_weight_aligned')

        self.task_align_for_mosaic(
            target = target,
            product = product,
            config = config,
            in_tags = in_tag_list,
            out_tags = out_tag_list,
            extra_ext_in = extra_ext_in,
            extra_ext_out = extra_ext_in,
            check_files = check_files,
            )
            
        self.task_linear_mosaic(
            target = target,
            product = product,
            config = config,
            image_tag = 'linmos_aligned',
            weight_tag = 'weight_aligned',
            out_tag = 'pbcorr_round',
            extra_ext_in = extra_ext_in,
            extra_ext_out = extra_ext_out,
            check_files = check_files,
            )

        if parts_have_singledish:
            self.task_linear_mosaic(
                target = target,
                product = product,
                config = config,
                image_tag = 'sd_aligned',
                weight_tag = 'sd_weight_aligned',
                out_tag = 'prepped_sd',
                extra_ext_in = extra_ext_in,
                extra_ext_out = extra_ext_out,
                check_files = check_files,
                )

        return()
 
    def recipe_cleanup_one_target(
        self,
        target = None,
        product = None,
        config = None,
        check_files = True,
        ext_ext = '',
        ):
        """
        Recipe that cleans up the output for one target, converting to
        Kelvin, compressing and trimming the cube and then exporting
        as a FITS file.
        """

        self.task_compress(
            target=target, config=config, product=product,
            check_files=check_files, do_pb_too=True,
            extra_ext_in=ext_ext, extra_ext_out=ext_ext,
            )

        self.task_convert_units(
            target=target, config=config, product=product,
            check_files=check_files,
            extra_ext_in=ext_ext, extra_ext_out=ext_ext,
            )

        self.task_export_to_fits(
            target=target, config=config, product=product,
            check_files=check_files, do_pb_too=True,
            extra_ext_in=ext_ext, extra_ext_out=ext_ext,
            )

        return()

    def recipe_convolve_to_scale(
        self,
        target = None,
        product = None,
        config = None,
        ext_ext = '',
        export_to_fits = True,
        check_files = True,
        ):
        """
        Convolve a target, product, config combination to a succession
        of angulars scale using the task that convolves to a round
        beam.
        """        

        res_list = self._kh.get_res_for_config(config)
        if res_list is None:
            logger.error("No target resolutions found for config "+config)
            return()

        for this_res in res_list:
            
            res_tag = self._kh.get_tag_for_res(this_res)
            
            # Check if the requested beam is smaller than the current one

            self.task_round_beam(
                target=target, config=config, product=product,
                in_tag = 'pbcorr_trimmed_k', out_tag = 'pbcorr_trimmed_k',
                extra_ext_in=ext_ext, extra_ext_out=ext_ext+'_res'+res_tag,
                force_beam_as=this_res,
                check_files=check_files
                )

            if export_to_fits:
                self.task_export_to_fits(
                    target=target, config=config, product=product,
                    check_files=check_files, do_pb_too=True,
                    extra_ext_in=ext_ext+'_res'+res_tag, extra_ext_out=ext_ext+'_res'+res_tag,
                    )

#endregion

#region Loops

    def loop_postprocess(
        self,
        do_prep=False,
        do_feather=False,
        do_mosaic=False,
        do_cleanup=False,
        do_convolve=False,
        feather_apod=False,
        feather_noapod=False,
        feather_before_mosaic=False,
        feather_after_mosaic=False,
        make_directories=True,
        ):
        """
        Loops over the full set of targets, products, and
        configurations to run the postprocessing. Toggle the parts of
        the loop using the do_XXX booleans. Other choices affect the
        algorithms used.
        """

        if len(self.get_targets()) == 0:            
            logger.error("Need a target list.")
            return(None)
 
        if len(self.get_all_products()) == 0:            
            logger.error("Need a products list.")
            return(None)

        if make_directories:
            self._kh.make_missing_directories(postprocess=True)

        # Prepare the interferometer data that has imaging. Includes
        # staging the single dish data, making weights, etc.
        
        if do_prep:

            for this_target, this_product, this_config in \
                    self.looper(do_targets=True,do_products=True,do_configs=True):
                           
                fname_dict = self._fname_dict(
                    target=this_target, product=this_product, config=this_config)
                        
                imaging_dir = self._kh.get_imaging_dir_for_target(this_target)
                has_imaging = os.path.isdir(imaging_dir + fname_dict['orig'])
                has_singledish = self._kh.has_singledish(target=this_target, product=this_product)        
                is_part_of_mosaic = self._kh.is_target_in_mosaic(this_target)

                if not has_imaging:
                    logger.debug("Skipping "+this_target+" because it lacks imaging.")
                    logger.debug(imaging_dir+fname_dict['orig'])
                    continue

                self.recipe_prep_one_target(
                    target = this_target, product = this_product, config = this_config,
                    check_files = True)

        # Feather the interferometer configuration data that has
        # imaging. We'll return to feather mosaicked intereferometer
        # and single dish data in the next steps.
                        
        if do_feather:

            for this_target, this_product, this_config in \
                    self.looper(do_targets=True,do_products=True,do_configs=True,just_interf=True):

                fname_dict = self._fname_dict(
                    target=this_target, product=this_product, config=this_config)
                        
                imaging_dir = self._kh.get_imaging_dir_for_target(this_target)
                has_imaging = os.path.isdir(imaging_dir + fname_dict['orig'])
                has_singledish = self._kh.has_singledish(target=this_target, product=this_product)        
                
                is_part_of_mosaic = self._kh.is_target_in_mosaic(this_target)
                if is_part_of_mosaic and not feather_before_mosaic:
                    logger.debug("Skipping "+this_target+" because feather_before_mosaic is False.")
                    continue
                            
                if not has_imaging:
                    logger.debug("Skipping "+this_target+" because it lacks imaging.")
                    logger.debug(imaging_dir+fname_dict['orig'])
                    continue
                            
                if not has_singledish:
                    logger.debug("Skipping "+this_target+" because it lacks single dish.")
                    continue

                if feather_apod:                            
                    self.task_feather(
                        target = this_target, product = this_product, config = this_config,
                        apodize=True, apod_ext='pb',extra_ext_out='_apod',check_files=True,
                        copy_weights=True,
                        )

                if feather_noapod:
                    self.task_feather(
                        target = this_target, product = this_product, config = this_config,
                        apodize=False, extra_ext_out='',check_files=True, 
                        copy_weights=True,
                        )
                    

        # Mosaic the intereferometer, single dish, and feathered data.

        if do_mosaic:
            
            # Loop over interferometer configurations
            
            for this_target, this_product, this_config in \
                    self.looper(do_targets=True,do_products=True,do_configs=True,just_interf=True):

                is_mosaic = self._kh.is_target_linmos(this_target)
                if not is_mosaic:
                    continue

                # Mosaic the interferometer data and the
                # single dish data (need to verify if parts
                # have single dish, enforce the same
                # astrometric grid).
                
                self.recipe_mosaic_one_target(
                    target = this_target, product = this_product, config = this_config,
                    check_files = True,
                    extra_ext_in = '',
                    extra_ext_out = '',
                    )

            # Loop over feather configurations
            
            for this_target, this_product, this_config in \
                    self.looper(do_targets=True,do_products=True,do_configs=True,just_feather=True):
                    
                # Mosaic the previously feathered data.

                if feather_apod:
                    self.recipe_mosaic_one_target(
                        target = this_target, product = this_product, config = this_config,
                        check_files = True,
                        extra_ext_in = '_apod',
                        extra_ext_out = '_prefeather_apod',
                        )

                if feather_noapod:
                    self.recipe_mosaic_one_target(
                        target = this_target, product = this_product, config = this_config,
                        check_files = True,
                        extra_ext_in = '',
                        extra_ext_out = '_prefeather',
                        )                            
                        
        # This round of feathering targets only mosaicked data. All
        # other data have been feathered above already.

        if do_feather:

            for this_target, this_product, this_config in \
                    self.looper(do_targets=True,do_products=True,do_configs=True,just_interf=True):
            
                is_mosaic = self._kh.is_target_linmos(this_target)
                if not is_mosaic:
                    continue

                if feather_apod:                  
                    self.task_feather(
                        target = this_target, product = this_product, config = this_config,
                        apodize=True, apod_ext='pb',extra_ext_out='_apod',check_files=True,
                        )
                            
                if feather_noapod:
                    self.task_feather(
                        target = this_target, product = this_product, config = this_config,
                        apodize=False, extra_ext_out='',check_files=True,
                        )
                        
        if do_cleanup:

            for this_target, this_product, this_config in \
                    self.looper(do_targets=True,do_products=True,do_configs=True):
                                        
                self.recipe_cleanup_one_target(
                    target = this_target, product = this_product, config = this_config,
                    check_files = True)
                
                        # TBD - different feather and mosaic ordering outputs? Maybe?                        

        if do_convolve:

            for this_target, this_product, this_config in \
                    self.looper(do_targets=True,do_products=True,do_configs=True):
                        
                self.recipe_convolve_to_scale(
                    target = this_target, product = this_product, config = this_config,
                    check_files = True, export_to_fits = True)

                
#endregion
