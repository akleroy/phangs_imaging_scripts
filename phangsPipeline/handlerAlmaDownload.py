"""AlmaDownloaderHandler

This module will query the ALMA archive, download and calibrate data as appropriate.
Since it uses astroquery, should be used outside monolithic CASA versions.

Example:
    $ ipython
    from phangsPipeline import handlerKeys as kh
    from phangsPipeline import handlerAlmaDownload as adl
    this_kh = kh.KeyHandler(master_key = 'phangsalma_keys/master_key.txt')
    this_adl = adl.DerivedHandler(key_handler=this_kh)
    this_adl.loop_alma_download()

Note that as this module will rewrite keys, you should re-initialise the key handler
after running through this module.

"""
import copy
import fnmatch
import glob
import logging
import multiprocessing as mp
import os
import re
import tarfile
from functools import partial

# Tags to identify antenna arrays
ANTENNA_ARRAY_SETUP = {'12m': ['DV', 'DA'],
                       '7m': ['CM'],
                       'tp': ['PM']}

ALLOWED_MS_GROUPINGS = ['mosaic', 'join', 'separate']


def query_target(target, max_query_failures=10):
    """Light wrapper around astroquery object query to allow for HTTP errors"""
    query_success = False
    query_failures = 0
    while not query_success and query_failures <= max_query_failures:
        try:
            observations = Alma.query_object(target)
            return observations
        except requests.exceptions.HTTPError:
            query_failures += 1

    # If we can't reach the server, fail out
    if not query_success:
        logger.warning('Maximum ALMA archive queries reached. Server is probably down.')
        raise Exception('Maximum ALMA archive queries reached. Server is probably down.')


def astroquery_download(row, cache_location=None):
    """Download ALMA archive files"""
    myAlma = Alma()

    if cache_location is not None:
        myAlma.cache_location = cache_location

    myAlma.download_files([row['access_url']], cache=True)
    return True


def get_casa_version_from_qa_report(qa_report):
    """Pulls CASA pipeline version from HTML QA report.

    Taken from Daizhong's A3COSMOS stuff, thanks Daizhong!

    """

    casa_pipeline_version = None

    with open(qa_report, 'r') as html:
        html_content = html.read()

        soup = BeautifulSoup(html_content, 'html.parser')

        for soup_th in soup.findAll('th', text='CASA Version'):
            for soup_td in soup_th.find_next_siblings():
                if len(soup_td.text) > 0:
                    if soup_td.text.find('.') > 0:
                        if soup_td.text.find(' ') > 0:
                            casa_pipeline_version = soup_td.text.split()[0]
                        else:
                            casa_pipeline_version = soup_td.text
                if casa_pipeline_version is not None:
                    break

            if casa_pipeline_version is not None:
                break

        # Try again to find span that includes 'CASA version:'
        if casa_pipeline_version is None:
            for soup_span in soup.findAll('span'):
                if soup_span.text.startswith('CASA version:'):
                    casa_pipeline_version = re.sub(r'CASA version: ([0-9.]+).*', r'\1', soup_span.text)
                    if casa_pipeline_version is not None:
                        break

        # Try again to find span that include 'CASA version' and next is 'Report date'
        if casa_pipeline_version is None:
            for soup_span in soup.findAll('span'):
                if soup_span.text == 'CASA version':
                    soup_tr = soup_span.find_parent('tr')
                    soup_tr_td = soup_tr.find_all('td')
                    soup_tr_td_text = [t.text.strip() for t in soup_tr_td]
                    soup_tr_td_index = soup_tr_td_text.index('CASA version')
                    soup_tr_next = soup_tr.find_next_sibling('tr')
                    soup_tr_next_td = soup_tr_next.find_all('td')
                    soup_tr_next_td_text = [t.text.strip() for t in soup_tr_next_td]
                    found_text = soup_tr_next_td_text[soup_tr_td_index]
                    if re.match(r'^[0-9.-]+$', found_text):
                        casa_pipeline_version = found_text
                    if casa_pipeline_version is not None:
                        break

    return casa_pipeline_version


def get_casa_version_from_weblog(weblog):
    """Pulls CASA pipeline version from weblog

    Again, taken from Daizhong's A3COSMOS stuff, thanks Daizhong!

    """

    casa_pipeline_version = None

    weblog_obj = tarfile.open(weblog)

    for weblog_item in weblog_obj.getmembers():
        if os.path.basename(weblog_item.name) == 'index.html':
            weblog_index_html = weblog_obj.extractfile(weblog_item)
            weblog_index_html_content = weblog_index_html.read()

            soup = BeautifulSoup(weblog_index_html_content, 'html.parser')

            for soup_th in soup.findAll('th', text='CASA Version'):
                for soup_td in soup_th.find_next_siblings():
                    if len(soup_td.text) > 0:
                        if soup_td.text.find('.') > 0:
                            if soup_td.text.find(' ') > 0:
                                casa_pipeline_version = soup_td.text.split()[0]
                            else:
                                casa_pipeline_version = soup_td.text
                    if casa_pipeline_version is not None:
                        break

                if casa_pipeline_version is not None:
                    break

            if casa_pipeline_version is not None:
                break

    return casa_pipeline_version


def get_casa_version_from_calibration_script(script):
    """Pulls CASA pipeline version from calibration script

    Again, taken from Daizhong's A3COSMOS stuff, thanks Daizhong!

    """

    casa_pipeline_version = None

    with open(script, 'r') as fp:

        for line in fp:
            if line.find('PLEASE USE THE SAME VERSION OF CASA') >= 0:
                if re.match(r'^.*PLEASE USE THE SAME VERSION OF CASA.*:[ \t]*([0-9.]+)\b.*$', line):
                    casa_pipeline_version = re.sub(r'^.*PLEASE USE THE SAME VERSION OF CASA.*:[ \t]*([0-9.]+)\b.*$',
                                                   r'\1', line)
                    # Trim any whitespace
                    casa_pipeline_version = casa_pipeline_version.strip()
                    break

    return casa_pipeline_version


logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

# Check casa environment by importing CASA-only packages
from .casa_check import is_casa_installed

casa_enabled = is_casa_installed()

if casa_enabled:
    logger.debug('casa_enabled = True')
else:
    logger.debug('casa_enabled = False')

# import phangs pipeline stuff
from . import utilsLines
from . import handlerTemplate
from . import handlerKeys

try:
    import astropy
    import astroquery
    from bs4 import BeautifulSoup
    import requests

    has_imports = True
except ImportError:
    logger.debug("Some required packages not installed.")
    has_imports = False

if has_imports:
    import astropy.units as u
    import numpy as np
    from astropy import table
    from astropy.coordinates import SkyCoord
    from astroquery.alma import Alma
    from astroquery.alma.utils import parse_frequency_support

    class AlmaDownloadHandler(handlerTemplate.HandlerTemplate):
        """
        Class to automate downloading and calibrating ALMA data.

        N.B. Since clean masks are generated externally, this code will not modify that key file. It will, however,
        modify ms_key, dir_key, linmos_key, singledish_key (TODO), and target_key.

        This script may also produce weirdness if you have a multiple files for each key type.

        TODO:
            Allow for TP
            Allow for proprietary data
            Add another key file for granular control over what we want to download. Something like
                {target}    {project_code:XXX, spectral_resolution:YYY} or so
            If target is not resolved then query will fail. Fall back to RA/Dec coordinates?
            Ongoing checks to catch different ways to pick up CASA pipeline versions.

            Just some thoughts...

            for the TP, maybe we do put that into a line directory just to make sorting out
            singledish_key easier later

        """

        ############
        # __init__ #
        ############

        def __init__(
                self,
                key_handler=None,
                master_key=None,
                restore_previous_target_keys=True,
                dry_run=False,
        ):
            # inherit template class
            handlerTemplate.HandlerTemplate.__init__(self,
                                                     key_handler=key_handler,
                                                     dry_run=dry_run)

            # If we've run this before, restore the original key file and reinitialise
            key_root = self._kh._key_dir
            target_key_file = self._kh._target_keys[0]
            target_key_path = os.path.join(key_root, target_key_file)
            pre_dl_target_key_path = target_key_path + '_pre_AlmaDownload'

            if os.path.exists(pre_dl_target_key_path) and restore_previous_target_keys:
                logger.info('This is a rerun, restoring old target keys')

                os.system('rm -rf %s' % target_key_path)
                os.system('mv %s %s' % (pre_dl_target_key_path, target_key_path))

                key_handler = handlerKeys.KeyHandler(master_key=master_key)
                handlerTemplate.HandlerTemplate.__init__(self,
                                                         key_handler=key_handler,
                                                         dry_run=dry_run)

        def loop_alma_download(
                self,
                do_all=False,
                make_directories=True,
                do_download=False,
                do_calibrate=False,
                do_build_file_keys=False,
                do_tp=False,
                split_ms='mosaic',
                overwrite_download=False,
                overwrite_calibrate=False,
                overwrite_build_file_keys=True,
                overwrite_all=False,
        ):
            """Download and calibrate ALMA data.

            Major steps:

            (1) Query archive for target/line/antenna config and download
            (2) Figure out which version of CASA to run the scriptForPI in, and run that
            (3) Build the MS file key from the downloaded, calibrated files

            N.B. After running this you will need to reinitialise the key handler, since this makes edits to various
            key files

            Args:
                do_tp (bool, optional): If True, will also download TP data and sort out for the TP pipeline.
                    TODO: CURRENTLY NOT IMPLEMENTED
                split_ms (str, optional): Can be one of 'mosaic', 'join', 'separate'. Determines how the MS file key is
                    set up for later staging and imaging. If mosaic, will attempt to be smart and split observations by
                    project ID/science goal. If 'join', all observations will be joined together into a potentially
                    huge map. If 'separate', then each calibrated MS will be assigned to a separate target. Defaults to
                    'mosaic', which mimics the way the PHANGS-ALMA targets are set up.

            """

            if do_all:
                make_directories = True
                do_download = True
                do_calibrate = True
                do_build_file_keys = True
            if overwrite_all:
                overwrite_download = True
                overwrite_calibrate = True
                overwrite_build_file_keys = True

            # Error checking

            if len(self.get_targets()) == 0:
                logger.error("Need a target list.")
                return None

            if len(self.get_all_products()) == 0:
                logger.error("Need a products list.")
                return None

            # If requested, make the directories

            if make_directories:
                self._kh.make_missing_directories(ms_root=True)

            if do_download:

                logger.info("")
                logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&")
                logger.info("Beginning download of ALMA data")
                logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&")
                logger.info("")

                for this_target, this_product, this_config in \
                        self.looper(do_targets=True, do_products=True, do_configs=True, just_interf=True):
                    uids = self.task_query(target=this_target,
                                           product=this_product,
                                           config=this_config)
                    self.task_download(target=this_target,
                                       product=this_product,
                                       config=this_config,
                                       uids=uids,
                                       overwrite=overwrite_download)

                if do_tp:

                    for this_target, this_product in \
                            self.looper(do_targets=True, do_products=True, do_configs=False):
                        uids = self.task_query(target=this_target,
                                               product=this_product,
                                               config='tp')
                        self.task_download(target=this_target,
                                           product=this_product,
                                           config='tp',
                                           uids=uids,
                                           overwrite=overwrite_download)

            if do_calibrate:

                logger.info("")
                logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&")
                logger.info("Beginning calibration of ALMA data")
                logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&")
                logger.info("")

                for this_target, this_config in \
                        self.looper(do_targets=True, do_products=False, do_configs=True, just_interf=True):
                    self.task_run_scriptforpi(target=this_target,
                                              config=this_config,
                                              overwrite=overwrite_calibrate)

            if do_build_file_keys:
                logger.info("")
                logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&")
                logger.info("Building file keys")
                logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&")
                logger.info("")

                self.task_build_file_keys(split_ms=split_ms,
                                          overwrite=overwrite_build_file_keys)

        def task_query(self,
                       target=None,
                       product=None,
                       config=None,
                       max_query_failures=10,
                       ):
            """Query ALMA archive.

            Uses astroquery to search the archive for data which matches (a) the target, (b) the product, (c) the
            antenna config (12m or 7m), and then velocity resolution.

            Args:
                max_query_failures (int, optional): Number of times to try querying the database before failing out.
                    Defaults to 10
            Returns:
                list of UIDs to download

            """

            logger.info("")
            logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&")
            logger.info("Querying ALMA database for:")
            logger.info('{0}, {1}, {2}'.format(target, product, config))
            logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&")
            logger.info("")

            if target is None:
                logger.warning('Require a target')
                return None
            if product is None:
                logger.warning('Require a product')
                return None
            if config is None or config not in ANTENNA_ARRAY_SETUP.keys():
                logger.warning('Require a valid config (%s)' % list(ANTENNA_ARRAY_SETUP.keys()))
                return None

            # Pull out line, channel width, and central frequency for this combination of target, product, config
            line = self._kh.get_line_tag_for_line_product(product)
            channel_kms = self._kh._config_dict['line_product'][product]['channel_kms']
            vsys, vwidth = self._kh.get_system_velocity_and_velocity_width_for_target(target, check_parent=False)
            line_low_ghz, line_high_ghz = utilsLines.get_ghz_range_for_line(line=line,
                                                                            vsys_kms=vsys,
                                                                            vwidth_kms=vwidth,
                                                                            )
            line_ghz = ((line_high_ghz + line_low_ghz) / 2) * u.GHz

            # Perform a query.
            observations = query_target(target, max_query_failures=max_query_failures)

            parsed_obs = table.Table()
            for observation in observations:

                # Do checks on velocity resolution
                velocity_res = observation['velocity_resolution'] / 1000
                if velocity_res > channel_kms:
                    continue

                # Checks on arrays we want
                antenna_arrays = observation['antenna_arrays']

                array_wanted = False

                array_setup_tags = ANTENNA_ARRAY_SETUP[config]
                for array_setup_tag in array_setup_tags:
                    if array_setup_tag in antenna_arrays:
                        array_wanted = True

                        # Sometimes it seems like TP can use other antenna (maybe this is only early science?) so double
                        # check here
                        if config != 'tp':
                            for tp_array_tag in ANTENNA_ARRAY_SETUP['tp']:
                                if tp_array_tag in antenna_arrays:
                                    array_wanted = False

                if not array_wanted:
                    continue

                # Check the line we want is in the frequency range

                freqs = parse_frequency_support(observation['frequency_support'])

                freq_wanted = False

                for freq in freqs:
                    if freq[0] <= line_ghz <= freq[1]:
                        freq_wanted = True

                if not freq_wanted:
                    continue

                # TODO: Include custom switches

                # I think that's everything we want, append that observation row to a new table
                parsed_obs = table.vstack([parsed_obs, observation])

            uids = np.unique(parsed_obs['member_ous_uid'])
            logger.info("")
            logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&")
            logger.info('Found %d suitable UIDs:' % len(uids))
            for uid in uids:
                logger.info('-> %s' % uid)
            logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&")
            logger.info("")

            return uids

        def task_download(self,
                          target=None,
                          config=None,
                          product=None,
                          uids=None,
                          n_simultaneous=5,
                          overwrite=False,
                          ):
            """Downloads queried UIDs.

            Hooks the queried UIDs to download raw/ancillary/readme data. Will go into a folder structure as
                ms_root/target/config/(product for TP), to make navigation a little easier

            Args:
                n_simultaneous (int, optional): Number of download processes to spawn, to speed up download time

            """

            logger.info("")
            logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&")
            logger.info("Downloading data for:")
            logger.info('{0}, {1}, {2}'.format(target, product, config))
            logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&")
            logger.info("")

            if target is None:
                logger.warning('Require a target')
                return None
            if config is None:
                logger.warning('Require a config')
                return None
            if product is None and config == 'tp':
                logger.warning('Require a product for TP')
                return None
            if uids is None:
                logger.warning('Require UIDs to download')
                return None

            # Setup download location, this will be ms_root/{target}/{config}
            ms_root = self._kh._ms_roots[0]

            if config == 'tp':
                dl_dir = os.path.join(ms_root, target, config, product)
            else:
                dl_dir = os.path.join(ms_root, target, config)

            # Delete download folder if we're overwriting
            if overwrite:
                os.system('rm -rf %s' % dl_dir)

            if not os.path.exists(dl_dir):
                os.makedirs(dl_dir)

            logger.info("")
            logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&")
            logger.info('Downloading ALMA data for %d UIDs:' % len(uids))
            for uid in uids:
                logger.info('-> %s' % uid)
            logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&")
            logger.info("")

            link_list = Alma.get_data_info(uids, expand_tarfiles=False)

            # Trim out any pointless files
            rows_to_remove = []
            for i, link in enumerate(link_list):
                if 'asdm.sdm' not in link['access_url'] and \
                        'README' not in link['access_url'] and \
                        'auxiliary' not in link['access_url']:
                    rows_to_remove.append(i)
            link_list.remove_rows(rows_to_remove)

            # Download files. Allow multiple downloads via pool
            with mp.Pool(n_simultaneous) as p:
                map_result = p.map_async(partial(astroquery_download, cache_location=dl_dir), link_list)
                map_result.wait()

            # Extract tar files
            original_dir = os.getcwd()
            os.chdir(dl_dir)

            tar_files = sorted(glob.glob('*.tar'))

            logger.info("")
            logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&")
            logger.info('Untarring %d files:' % len(tar_files))
            for tar_file in tar_files:
                logger.info('-> %s' % tar_file)
            logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&")
            logger.info("")

            os.system('cat *.tar | tar --skip-old-files -xf - -i')
            os.chdir(original_dir)

        def task_run_scriptforpi(self,
                                 target=None,
                                 config=None,
                                 suppress_casa_output=True,
                                 overwrite=False):
            """Runs scriptForPI on downloaded datasets.

            Figures out the pipeline version of CASA to run (currently, just from the QA report), and then will execute
            that script to produce a calibrated measurement set

            Args:
                suppress_casa_output (bool, optional): If True, will suppress most CASA output while running
                    scriptForPI. Mostly useful for debugging, and will still put the whole log into the /script
                    directory. Defaults to True.

            """

            logger.info("")
            logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&")
            logger.info("Calibrating data for:")
            logger.info('{0}, {1}'.format(target, config))
            logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&")
            logger.info("")

            if target is None:
                logger.warning('Require a target')
                return None
            if config is None or config not in ANTENNA_ARRAY_SETUP.keys():
                logger.warning('Require a valid config (%s)' % list(ANTENNA_ARRAY_SETUP.keys()))
                return None

            ms_root = self._kh._ms_roots[0]
            dl_dir = os.path.join(ms_root, target, config)

            original_dir = os.getcwd()
            os.chdir(dl_dir)

            # First, search for QA reports
            casa_version_files = {}
            for root, dirnames, filenames in os.walk(os.getcwd()):
                for filename in fnmatch.filter(filenames, '*.qa2_report.html'):
                    casa_version_files[os.path.join(root, filename)] = 'qa_report'

            # Search for weblogs
            for root, dirnames, filenames in os.walk(os.getcwd()):
                for filename in fnmatch.filter(filenames, '*.weblog.tgz'):
                    casa_version_files[os.path.join(root, filename)] = 'weblog'

            # Search for calibration scripts
            calibration_dirs = []
            for root, dirnames, filenames in os.walk(os.getcwd()):
                for filename in fnmatch.filter(filenames, '*.scriptForCalibration.py'):
                    # There may be multiple calibration scripts so make sure we only use one per directory
                    if root not in calibration_dirs:
                        casa_version_files[os.path.join(root, filename)] = 'calibration_script'
                        calibration_dirs.append(root)

            for casa_version_file in casa_version_files.keys():

                casa_version_file_type = casa_version_files[casa_version_file]

                member_dir = os.path.sep + os.path.join(*casa_version_file.split(os.path.sep)[:-2])
                calibrated_dir = os.path.join(member_dir, 'calibrated')
                script_dir = os.path.join(member_dir, 'script')

                if os.path.exists(calibrated_dir) and not overwrite:
                    logger.info("")
                    logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&")
                    logger.info('Data already calibrated, will not rerun for %s' % member_dir.replace(ms_root, ''))
                    logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&")
                    logger.info("")
                    continue

                # Clear out the calibrated directory
                os.system('rm -rf %s' % calibrated_dir)

                if casa_version_file_type == 'qa_report':
                    casa_pipeline_version = get_casa_version_from_qa_report(casa_version_file)
                elif casa_version_file_type == 'weblog':
                    casa_pipeline_version = get_casa_version_from_weblog(casa_version_file)
                elif casa_version_file_type == 'calibration_script':
                    casa_pipeline_version = get_casa_version_from_calibration_script(casa_version_file)
                else:
                    logger.warning('Version file type %s not understood' % casa_version_file_type)
                    raise Exception('Version file type %s not understood' % casa_version_file_type)

                if casa_pipeline_version is None:
                    logger.warning('Could not find a CASA version to run scriptForPI')
                    raise Exception('Could not find a CASA version to run scriptForPI')

                # Get the CASA path
                casa_path = self._kh.get_path_for_casaversion(casa_pipeline_version)

                if casa_path is None:
                    logger.warning('No CASA path defined for %s' % casa_pipeline_version)
                    raise Exception('No CASA path defined for %s' % casa_pipeline_version)

                # Make sure we're properly pointing at CASA
                if not casa_path.endswith('casa'):
                    casa_path = os.path.join(casa_path, 'casa')

                os.chdir(script_dir)
                if os.path.exists('scriptForPI.py'):
                    script_names = ['scriptForPI.py']
                else:
                    script_names = glob.glob('*scriptForPI.py')

                logger.info("")
                logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&")
                logger.info('Calibrating %s' % member_dir.replace(ms_root, ''))
                logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&")
                logger.info("")

                for script_name in script_names:
                    logger.info("")
                    logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&")
                    logger.info('Running %s in CASA %s' % (script_name, casa_pipeline_version))
                    logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&")
                    logger.info("")
                    cmd = '%s --pipeline --nologger -c %s' % (casa_path, script_name)
                    if suppress_casa_output:
                        cmd += ' >/dev/null 2>&1'
                    os.system(cmd)

                os.chdir(dl_dir)

            os.chdir(dl_dir)

            # Check everything now has a calibrated dir, which indicates the script has run
            calibrated_dirs = []
            for path, dirs, files in os.walk(os.getcwd()):
                if 'calibration' in dirs:
                    calibrated_dirs.append(os.path.join(path, 'calibrated'))

            missing_dirs = False
            for calibrated_dir in calibrated_dirs:
                if not os.path.exists(calibrated_dir):
                    # TODO: DEBUG AS THIS COMES UP
                    logger.warning('Unexpected missing calibrated directory! %s' % calibrated_dir)
                    missing_dirs = True

            if missing_dirs:
                raise Exception('Missing some calibrated directories!')

            # Move back to the original directory
            os.chdir(original_dir)

        def task_build_file_keys(self,
                                 split_ms='mosaic',
                                 overwrite=False):
            """Builds MS file key from calibrated measurement sets.

            Recursively search through calibrated measurement sets and build up into a key file to be read in for later
            imaging.

            """

            if split_ms not in ALLOWED_MS_GROUPINGS:
                logger.warning('split_ms=%s not allowed. Should be one of %s' % (split_ms, ALLOWED_MS_GROUPINGS))
                return None

            key_root = self._kh._key_dir
            ms_root = self._kh._ms_roots[0]

            # Pull out file names
            dir_key_file = self._kh._dir_keys[0]
            linmos_key_file = self._kh._linmos_keys[0]
            ms_key_file = self._kh._ms_keys[0]
            target_key_file = self._kh._target_keys[0]

            dir_key_file_name = os.path.join(key_root, dir_key_file)
            linmos_key_file_name = os.path.join(key_root, linmos_key_file)
            ms_file_name = os.path.join(key_root, ms_key_file)
            target_file_name = os.path.join(key_root, target_key_file)

            file_names_stacked = [dir_key_file_name,
                                  linmos_key_file_name,
                                  ms_file_name]

            if any([os.path.exists(file_name) for file_name in file_names_stacked]) and not overwrite:
                logger.info('Files already exists and overwrite is False. Will not overwrite')
                return None

            os.system('rm -rf %s' % ms_file_name)

            ms_file = open(ms_file_name, 'w+')

            mosaic_info = {}

            targets = self._kh.get_targets()

            for target in targets:

                observations = query_target(target)

                dl_dir = os.path.join(ms_root, target)

                original_dir = os.getcwd()
                os.chdir(dl_dir)

                target_ms_dict = {}
                # Find target MSs
                all_configs = []
                all_project_ids = []
                all_science_goals = []
                for root, dirnames, filenames in os.walk(os.getcwd()):
                    for dirname in fnmatch.filter(dirnames, '*.ms.split.cal'):
                        # TODO: This will miss not-split MSs, so double check
                        full_dir = os.path.join(root, dirname).split(ms_root)[1]
                        full_dir_split = full_dir.split(os.path.sep)

                        # Pull out config, project ID, science goal, member UID
                        config = full_dir_split[1]
                        project_id = full_dir_split[2]
                        science_goal = full_dir_split[3]

                        # Use member uid to get at coordinates of the observation

                        member_uid = full_dir_split[5]
                        member_uid = member_uid.replace('___', '://').replace('_', '/').replace('member.', '')
                        observations_ms = observations[member_uid == observations['member_ous_uid']][0]
                        ra, dec = observations_ms['s_ra'], observations_ms['s_dec']

                        coord = SkyCoord(ra=ra * u.deg, dec=dec * u.deg)
                        ra_str, dec_str = coord.to_string('hmsdms').split()

                        if config not in all_configs:
                            all_configs.append(config)
                        if project_id not in all_project_ids:
                            all_project_ids.append(project_id)
                        if science_goal not in all_science_goals:
                            all_science_goals.append(science_goal)

                        target_ms_dict[full_dir] = [config, project_id, science_goal, ra_str, dec_str]

                # Start writing these things out
                mosaic_no = 1
                if split_ms == 'mosaic':

                    # For mosaics, group up by project ID/science goal. This captures mosaics at different times and
                    # multiple mosaics within a single project

                    using_mosaic = False

                    if len(all_project_ids) == 1 and len(all_science_goals) == 1:
                        target_key = copy.deepcopy(target)
                    else:
                        target_key = '%s_%d' % (target, mosaic_no)
                        using_mosaic = True

                    for project_id in all_project_ids:
                        for science_goal in all_science_goals:
                            observation_number = {config: 1 for config in all_configs}
                            match_found = False
                            match_key = None
                            for key in target_ms_dict.keys():
                                if target_ms_dict[key][1] == project_id and \
                                        target_ms_dict[key][2] == science_goal:
                                    match_found = True
                                    match_key = copy.deepcopy(key)
                                    config = target_ms_dict[key][0]
                                    ms_file.write('%s\t%s\tall\t%s\t%d\t%s\n' %
                                                  (target_key, project_id, config, observation_number[config], key))
                                    observation_number[config] += 1
                            if match_found:

                                ra, dec = target_ms_dict[match_key][3], target_ms_dict[match_key][4]

                                if using_mosaic:
                                    mosaic_info[target_key] = {}
                                    mosaic_info[target_key]['original_target'] = target
                                    mosaic_info[target_key]['ra'] = ra
                                    mosaic_info[target_key]['dec'] = dec

                                mosaic_no += 1
                                target_key = '%s_%d' % (target, mosaic_no)
                                ms_file.write('\n')

                elif split_ms == 'join':

                    # For join, collapse all observations for a single target down so imaging will clean everything at
                    # once

                    target_key = copy.deepcopy(target)
                    observation_number = {config: 1 for config in all_configs}
                    for key in target_ms_dict.keys():
                        config = target_ms_dict[key][0]
                        project_id = target_ms_dict[key][1]
                        ms_file.write('%s\t%s\tall\t%s\t%d\t%s\n' %
                                      (target_key, project_id, config, observation_number[config], key))
                        observation_number[config] += 1

                    ms_file.write('\n')

                elif split_ms == 'separate':

                    # For separate, every observation will be imaged separately

                    target_key = '%s_%d' % (target, mosaic_no)
                    observation_number = 1
                    for key in target_ms_dict.keys():
                        config = target_ms_dict[key][0]
                        project_id = target_ms_dict[key][1]

                        ms_file.write('%s\t%s\tall\t%s\t%d\t%s\n' %
                                      (target_key, project_id, config, observation_number, key))
                        ms_file.write('\n')

                        ra, dec = target_ms_dict[key][3], target_ms_dict[key][4]

                        mosaic_info[target_key] = {}
                        mosaic_info[target_key]['original_target'] = target
                        mosaic_info[target_key]['ra'] = ra
                        mosaic_info[target_key]['dec'] = dec

                        mosaic_no += 1
                        target_key = '%s_%d' % (target, mosaic_no)

                os.chdir(original_dir)

            ms_file.close()

            # Write out dir key

            os.system('rm -rf %s' % dir_key_file_name)
            dir_file = open(dir_key_file_name, 'w+')

            for key in mosaic_info.keys():

                target_name = mosaic_info[key]['original_target']
                dir_file.write('%s\t%s\n' % (key, target_name))

            dir_file.write('\n')
            dir_file.close()

            # Write out linmos keys

            os.system('rm -rf %s' % linmos_key_file_name)
            linmos_file = open(linmos_key_file_name, 'w+')

            for key in mosaic_info.keys():
                target_name = mosaic_info[key]['original_target']
                linmos_file.write('%s\t%s\n' % (target_name, key))

            linmos_file.write('\n')
            linmos_file.close()

            # Write out target definitions. Start by renaming the original, so we can pull it back in if this gets rerun

            os.system('rm -rf %s_pre_AlmaDownload' % target_file_name)
            os.system('mv %s %s_pre_AlmaDownload' % (target_file_name, target_file_name))

            target_file = open(target_file_name, 'w+')

            for target in targets:

                ra, dec = self._kh.get_phasecenter_for_target(target=target)

                # Pull out velocity and width
                vel, vel_width = self._kh.get_system_velocity_and_velocity_width_for_target(target=target)

                target_file.write('%s\t%s\t%s\t%s\t%s\n' % (target, ra, dec, vel, vel_width))

                for target_mosaic in mosaic_info.keys():
                    if mosaic_info[target_mosaic]['original_target'] == target:

                        ra = mosaic_info[target_mosaic]['ra']
                        dec = mosaic_info[target_mosaic]['dec']

                        target_file.write('%s\t%s\t%s\t%s\t%s\n' % (target_mosaic, ra, dec, vel, vel_width))

            target_file.close()
