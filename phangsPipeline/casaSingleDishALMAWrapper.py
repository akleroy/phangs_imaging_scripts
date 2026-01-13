

import os
import glob
import tarfile
import shutil

import logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

from . import casaLegacySingleDishRoutines as csdr

# path constants
path_script = '../script/'       # Path to the script folder.
path_raw    = '../raw/'          # Path to the raw folder.
path_dataproduct = '../data/'    # Path to data products.


def extractJyperK():
    logger.info("Extracting Jy per K conversion factor")

    file_script = 'jyperk.csv'

    if not os.path.isfile(file_script):
        filetgz = glob.glob("*auxproducts.tgz")
        if len(filetgz) > 0:
            tar = tarfile.open(filetgz[0])
            tar.extractall()
            tar.close()

        # Check if it's still not there or search for jyperk_query.csv
        if not os.path.exists(file_script):
            alt_file_script = 'jyperk_query.csv'
            if not os.path.exists(alt_file_script):
                raise ValueError("Unable to find the jyperk.csv or jyperk_query.csv files.")
            else:
                file_script = alt_file_script

    # Copying the jyperk.csv file to the working directory
    shutil.copyfile(file_script, 'jyperk.csv')


def runALMAPipeline(path_galaxy,
                    source='all',
                    baseline_fit_func='poly',
                    baseline_fit_order=1,
                    baseline_linewindowmode='replace',
                    baseline_linewindow=None,
                    product_dict=None):
    '''
    Wraps the ALMA SD pipeline recipes into a single function. This it to pass
    custom parameters to the pipeline, primarily the baseline correction.

    Parameters
    ----------
    path_galaxy : str
        Path to the galaxy folder.
    baseline_fit_func : str
        Baseline fit function.
    baseline_fit_order : int
        Baseline fit order.
    baseline_linewindowmode : str
        Baseline line window mode.
    baseline_linewindow : dict
        Dictionary of baseline line windows. e.g. {"23": ['230.62GHz', '230.78GHz']}

    '''


    # Change to working directory:
    logger.info("> Changing directory to "+os.path.join(path_galaxy,'raw'))
    os.chdir(os.path.join(path_galaxy,'calibration'))   # Working on the calibration folder of the current galaxy


    # Defining Execution Blocks (EBS) names
    EBsnames = [f for f in os.listdir(".") if f.endswith('.asdm.sdm')]

    # Retrieve the k2jy scaling from the auxiliary calibration products.
    extractJyperK()

    # Drop the '.asdm.sdm' extension. This matches with the naming in the JyperK file
    # and the format in the provided pipeline script.
    newEBnames = []
    for EBname in EBsnames:
        newEBname = EBname.split('.asdm.sdm')[0]

        os.rename(EBname, newEBname)

        newEBnames.append(newEBname)

    EBsnames = newEBnames

    # Create baseline dict with freq ranges to mask:
    if baseline_linewindow is None:
        import astropy.units as u

        # Construct the line window via spw:low~high strings based on the
        # product_dict.
        baseline_linewindow = []
        for this_product in product_dict:

            freq_rest = product_dict[this_product]['freq_rest_MHz']
            vsys = product_dict[this_product]['vsys']
            vel_line_mask = product_dict[this_product]['vel_line_mask']

            # Convert velocity range to frequency range
            freq_line_mask = (vel_line_mask * u.km / u.s).to(u.Hz, u.doppler_optical(freq_rest * u.MHz)).value

            # Giving list of [low, high] frequency ranges. The pipeline will apply the ranges to all valid SPWs.
            baseline_linewindow.append(list(freq_line_mask))

    # Run pipeline recipe.
    context = h_init()
    try:
        hsd_importdata(vis=EBsnames)

        hsd_flagdata(pipelinemode="automatic")
        h_tsyscal(pipelinemode="automatic")
        hsd_tsysflag(pipelinemode="automatic")
        hsd_skycal(pipelinemode="automatic")
        hsd_k2jycal(dbservice=False)
        hsd_applycal(pipelinemode="automatic")
        hsd_atmcor(pipelinemode="automatic")

        # NOTE: single baseline fit only. The pipeline repeats after 1 flagging step.
        # May need to revisit if flagging is missed.
        hsd_baseline(pipelinemode="automatic",
                     fit_func=baseline_fit_func,
                     fit_order=baseline_fit_order,
                     linewindowmode=baseline_linewindowmode,
                     linewindow=baseline_linewindow
                     )
        hsd_blflag(pipelinemode="automatic")


    finally:
        h_save()


    # TODO: Add imaging with requested spw, channel and range.
    # TODO: if possible, we should just loop through all required products in 1 go.
    # Per product for SD imaging is cheap and fast.

    # TODO: Export the final imaging products
