"""AlmaDownloaderHandler

This module will query the ALMA archive, download and calibrate data as appropriate.
Since it uses astroquery, should be used outside of monolithic CASA versions.

Example:
    $ ipython
    from phangsPipeline import handlerKeys as kh
    from phangsPipeline import handlerAlmaDownload as adl
    this_kh = kh.KeyHandler(master_key = 'phangsalma_keys/master_key.txt')
    this_dh = dh.DerivedHandler(key_handler = this_kh)
    this_dh.set_targets(only = ['ngc0628', 'ngc2997', 'ngc4321'])
    this_dh.set_interf_configs(only = ['7m'])
    this_dh.set_line_products(only = ['co21'])
    this_dh.loop_alma_download()

"""

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
from . import utilsResolutions
from . import utilsFilenames
from . import utilsLines
from . import handlerTemplate

try:
    import astropy
    import astroquery

    has_astropy = True
except ImportError:
    logger.debug("astropy not installed.")
    has_astropy = False

if has_astropy:
    import astropy.units as u
    from astropy import table
    from astroquery.alma import Alma
    from astroquery.alma.utils import parse_frequency_support


    class DerivedHandler(handlerTemplate.HandlerTemplate):
        """
        TODO
        """

        ############
        # __init__ #
        ############

        def __init__(
                self,
                key_handler=None,
                dry_run=False,
        ):
            # inherit template class
            handlerTemplate.HandlerTemplate.__init__(self,
                                                     key_handler=key_handler,
                                                     dry_run=dry_run)

        def loop_alma_download(
                self,
                do_all=False,
                make_directories=True,
                do_download=False,
                do_calibrate=False,
                overwrite=False,
        ):

            if do_all:
                make_directories = True
                do_download = True
                do_calibrate = True

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


# ANTENNA_ARRAY_SETUP = {'12m': ['DV', 'DA'],
#                        '7m': ['CM'],
#                        'tp': ['PM']}
#
# observations = Alma.query_object('NGC5068')
#
# min_vel_res = 2.5
# target_freq = 230 * u.GHz
# arrays = ['12m']
#
# parsed_obs = table.Table()
#
# for observation in observations:
#
#     # Do checks on velocity resolution
#     velocity_res = observation['velocity_resolution'] / 1000
#     if velocity_res > min_vel_res:
#         continue
#
#     # Checks on arrays we want
#     antenna_arrays = observation['antenna_arrays']
#
#     array_wanted = False
#
#     for array in arrays:
#         array_setup_tags = ANTENNA_ARRAY_SETUP[array]
#         for array_setup_tag in array_setup_tags:
#             if array_setup_tag in antenna_arrays:
#                 array_wanted = True
#
#     if not array_wanted:
#         continue
#
#     # Check the line we want is in the frequency range
#
#     freqs = parse_frequency_support(observation['frequency_support'])
#
#     freq_wanted = False
#
#     for freq in freqs:
#         if freq[0] <= target_freq <= freq[1]:
#             freq_wanted = True
#
#     if not freq_wanted:
#         continue
#
#     # I think that's everything we want, append that observation row to a new table
#     parsed_obs = table.vstack([parsed_obs, observation])
#
# uids = np.unique(parsed_obs['member_ous_uid'])
#
# link_list = Alma.get_data_info(uids, expand_tarfiles=True)
#
# # myAlma = Alma()
# # myAlma.cache_location = '/big/external/drive/'
# # myAlma.download_files(link_list, cache=True)
# # CASA -c scriptforpi.py or something