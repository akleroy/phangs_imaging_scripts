
Scripts in this folder produce shuffled cubes and flat maps of the PHANGS-ALMA CO cubes,
using an external velocity field and a velocity-integration mask that combines a fixed
velocity window with a CO-based line emission mask.

### WORKFLOW FOR MOST USERS

If you want to produce shuffled cubes and flat maps of PHANSG-ALMA CO line cubes, you either produce a velocity field using `prepare_velocity_field.py` or input your on velocity field and directly run `shuffled_cubes_pipeline.py`.

0. If no velocity field exists beforehand: Execute `prepare_velocity_field.py` from the terminal with python3 to produce a velocity field by providing several additional velocity fields, e.g. CO, Halpha, HI moment-1 maps, or modelled velocity fields. The script can combine various velocity maps into a master velocity map.
1. Execute `shuffled_cubes_pipeline.py` from the terminal with python3 to produce shuffled cubes and flat maps.

### PIPELINE PRODUCTS

The pipeline produces shuffled cubes, noise cubes, a set of flat maps with corresponding uncertainties and velocity-integration masks:

Products:

* cubes:
    1. shuffled cube
    2. shuffled noise cube
    3. shuffled masks (two versions: narrow_strict, wide_broad)
    4. re-shuffled masks (same as shuffled masks but shuffled back to the original velocity field)
* maps: 
    1. narrow (+-50 km/s) fixed window flat maps
    2. wide (+-100 km/s) fixed window flat maps
    3. narrow (+-20-200 km/s; adapted to each galaxy) fixed window combined with a strict CO line emission mask
    4. narrow (+-20-200 km/s; adapted to each galaxy) fixed window combined with a broad CO line emission mask
    5. wide (+-100 km/s) fixed window combined with a strict CO line emission mask
    6. wide (+-100 km/s) fixed window combined with a broad CO line emission mask

### PYTHON PACKAGES

* Python 3.9 or later
    * [numpy](https://numpy.org)
    * [astropy](https://www.astropy.org)
    * [matplotlib](https://matplotlib.org)
    * [pandas](https://pandas.pydata.org)
    * [reproject](https://reproject.readthedocs.io)

* Run with: 
    * python 3.9
    * numpy 1.26.4
    * astropy 5.2.2
    * matplotlib 3.5.3 
    * pandas 2.2.3
    * reproject 0.13.0

### SCRIPTS

- `prepare_velocity_field.py` is used to create a master velocity field from various tracers (e.g. CO, Halpha, HI, model) that are combined in a hirarchhical order
- `shuffled_cubes_pipeline.py` produces the shuffled cubes and flat maps
- `ancillary_functions.py` contains ancillary functions used in the main pipeline script
