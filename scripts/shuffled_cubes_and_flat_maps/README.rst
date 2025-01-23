
Scripts in this folder produce shuffled cubes and flat maps of the PHANGS-ALMA CO cubes,
using an external velocity field and a velocity-integration mask that combines a fixed
velocity window with a CO-based line emission mask.

- `prepare_velocity_field.py` is used to create a master velocity field from various tracers (e.g. CO, Halpha, HI, model) that are combined in a hirarchhical order
- `shuffled_cubes_pipeline.py` produces the shuffled cubes and flat maps
- `ancillary_functions.py` contains ancillary functions used in the main pipeline script