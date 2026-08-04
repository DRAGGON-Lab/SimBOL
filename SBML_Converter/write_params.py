"""
generate_params.py

Builds the params.json consumed by sbml_to_cellmodeller.py.

Meant to be run cell-by-cell in a notebook: edit the values in the
PARAMS dict below to match your circuit, run the cell, and a params.json
file is written next to this script (path is set in OUTPUT_PATH).

Every field has a default already applied inside sbml_to_cellmodeller.py
if you omit it — so you only need to fill in what matters for your model.
Delete keys you don't need; you don't have to keep every section.
"""

import json

# WHERE TO WRITE THE FILE

OUTPUT_PATH = "params.json"


# PARAMETERS
# Edit these to match your SBML model and desired simulation setup.

PARAMS = {

    # Overall simulation settings
    "simulation": {
        "max_cells":    5000,     # pre-allocated array size; raise if the colony will grow larger
        "jitter_z":     False,    # True = let cells jiggle out of the 2D plane a bit
        "gamma":        100.0,    # physical stiffness/friction constant for CLBacterium
        "pickle_steps": 10,       # save a snapshot every N simulation steps
        "random_seed":  42,       # set to null / omit for a different run each time
    },

    # One entry per cell type / strain in the simulation.
    # "initial_concentrations" keys are SBML species IDs (not the safe_name
    # Python variable names) — use the *_ids_ from your SBML file, e.g. "GFP".
    "cell_types": [
        {
            "display_name":         "WT",
            "color":                [0.1, 0.9, 0.1],   # RGB, 0-1
            "division_length":      3.5,
            "growth_rate":          1.0,
            "division_noise":       0.005,
            "initial_pos":          [0.0, 0.0, 0.0],
            "initial_dir":          [1.0, 0.0, 0.0],
            "initial_concentrations": {
                # "GFP": 0.0,
            },
        },
    ],

    # Intercellular signalling (GridDiffusion + CLCrankNicIntegrator).
    # Only species living in a compartment classified as "diffusible" ever
    # reach this — see sbml_mapping below to control that classification.
    "signaling": {
        "enabled":             True,
        "grid_len":             60,        # grid cells along x and y
        "grid_z_cells":         3,         # grid cells along z (small = thin/2D colony)
        "grid_size":            4.0,       # microns per grid cell (must be isotropic)
        # "grid_origin":       [x, y, z],  # defaults to centering the grid on (0,0,0)
        "boundary_condition":   "reflect", # one of: reflect, constant, wrap, mirror, nearest
        "signals": [
            # Per-signal overrides. "name" must match the SBML species id.
            # Anything omitted here falls back to a diffusion_rate of 0.1
            # and a membrane_exchange_rate equal to the diffusion_rate.
            # {
            #     "name":                   "AHL",
            #     "diffusion_rate":         0.1,
            #     "membrane_exchange_rate": 0.1,
            #     "initial_concentration":  0.0,
            # },
        ],
    },

    # Overrides for how SBML compartments/species get classified.
    "sbml_mapping": {
        # "diffusible_compartments": ["extracellular"],
        # "local_compartments":      ["cytoplasm"],
        # "ignore_species_ids":      ["H2O", "ATP"],
    },

    # Fixed numeric values for boundaryCondition/constant SBML species
    # (keyed by SBML species id). Omit to use each species' SBML
    # initialConcentration/initialAmount as the fixed value.
    "species_overrides": {
        # "Inducer": 1.0,
    },

    # Optional physical walls for the simulation (e.g. a microfluidic trap).
    "walls": [
        # {"point": [0, -20, 0], "normal": [0, 1, 0], "coeff": 1.0},
    ],
}


# WRITE THE FILE

def write_params(params=PARAMS, path=OUTPUT_PATH):
    with open(path, "w") as f:
        json.dump(params, f, indent=2)
    print(f"Wrote {path}")
    return path


if __name__ == "__main__":
    write_params()
