"""
params.py

Builds the `parameters` dict consumed by `generate_script()` (and `main()`)
in cellmodeller_converter.py, and provides fluorescent-protein colour hints
used by the UI form in ui_params.py.
"""

from cellmodeller_converter import (
    parse_json,
    detect_signaling_topology,
)


# Fluorescent-protein keyword / BioBrick to colour-name lookup.
#
# This is *not* used by generate_script() itself (it has its own RGB lookup
# in _FP_KEYWORDS / _FP_BIOBRICK for per-cell shading by expression level).
# It exists purely so the UI can suggest a sensible default strain colour
# when a reporter protein is detected in the circuit.

PROTEIN_REPORTER_COLOR_NAMES = {
    # Green
    "BBa_E0040": "green", "BBa_K2148009": "green", "BBa_K2560042": "green",
    "BBa_K2009820": "green", "BBa_K4946001": "green", "BBa_K4159005": "green",
    "GFP": "green", "gfp": "green", "EGFP": "green", "eGFP": "green",
    "eGFPuv": "green", "mNeonGreen": "green", "mEmerald": "green",
    "mGFP": "green", "esmGFP": "green", "mGFP2": "green", "avGFP": "green",

    # Red
    "BBa_E1010": "red", "BBa_K1323009": "red", "BBa_K1688019": "red",
    "BBa_K3128008": "red", "BBa_K1399000": "red", "BBa_K1399001": "red",
    "BBa_K1399002": "red", "BBa_K3841014": "red",
    "RFP": "red", "rfp": "red", "DsRED": "red", "dsRed": "red",
    "mCherry": "red", "mStrawberry": "red", "mCherry2": "red",
    "mRFP1": "red", "mRFP2": "red", "mNeptune": "red", "mRuby": "red",
    "mKate": "red", "mNeptune2": "red", "mStrawberry2": "red",
    "mStrawberry3": "red", "mKate2": "red", "mRuby2": "red",
    "mOrange": "orange",

    # Yellow
    "BBa_E0030": "yellow", "BBa_K592101": "yellow", "BBa_K2656020": "yellow",
    "BBa_K3427000": "yellow", "BBa_K2656021": "yellow", "BBa_K1323010": "yellow",
    "BBa_K165005": "yellow",
    "YFP": "yellow", "yfp": "yellow", "EYFP": "yellow", "eYFP": "yellow",
    "mCitrine": "yellow", "mVenus": "yellow", "mYFP": "yellow",
    "mVenus2": "yellow", "YFP_LVA": "yellow", "YFP_LOV": "yellow",

    # Cyan
    "BBa_E0020": "cyan", "CFP": "cyan", "cfp": "cyan", "ECFP": "cyan",
    "eCFP": "cyan", "mTurquoise": "cyan", "mTFP1": "cyan", "mTurq2": "cyan",
    "mTFP2": "cyan", "mCFP": "cyan", "mTurquoise3": "cyan", "mTFP3": "cyan",
}

# Colour name to RGB triple (0-1 floats), used when seeding a cell_type's
# default "color" field from a detected reporter.

COLOR_NAME_TO_RGB = {
    "green":  [0.0, 1.0, 0.2],
    "red":    [1.0, 0.1, 0.1],
    "yellow": [0.9, 0.9, 0.0],
    "cyan":   [0.1, 0.8, 1.0],
    "orange": [1.0, 0.5, 0.0],
}


def detect_reporter_color(protein_display_id):
    """
    Best-effort colour-name guess for a reporter protein display_id.
    Tries an exact match first (BioBrick ID or common name), then a
    case-insensitive substring match (so e.g. "GFP_LVA" still matches
    "gfp"). Returns None if nothing matches.
    """
    if protein_display_id in PROTEIN_REPORTER_COLOR_NAMES:
        return PROTEIN_REPORTER_COLOR_NAMES[protein_display_id]
    lowered = protein_display_id.lower()
    for keyword, color in PROTEIN_REPORTER_COLOR_NAMES.items():
        if keyword.lower() in lowered:
            return color
    return None


# Defaults, matching the fallback values used internally

DEFAULT_SIMULATION = {
    "max_cells":    10000,
    "jitter_z":     False,
    "gamma":        100.0,
    "pickle_steps": 10,
    "random_seed":  None,
}

DEFAULT_CELL_TYPE = {
    "display_name":    "Strain 0",
    "color":            [1.0, 0.3, 0.3],
    "division_length":  3.5,
    "growth_rate":      1.0,
    "division_noise":   0.005,
    "initial_pos":      [0.0, 0.0, 0.0],
    "initial_dir":      [1.0, 0.0, 0.0],
    "initial_concentrations": {},
}

DEFAULT_KINETICS = {
    "production_rate":        2.0,
    "max_production_rate":    2.0,
    "degradation_rate":       0.1,
    "hill_coefficient":       4.0,
    "activation_threshold":   2.0,
    "repression_threshold":   2.0,
    "signal_production_rate": 0.1,
}

DEFAULT_SIGNAL_RATES = {
    "diffusion_rate":   0.1,
    "degradation_rate": 0.05,
}

DEFAULT_SIGNALING = {
    "enabled":            True,
    "grid_len":           100,
    "grid_size":          4.0,
    "boundary_condition": "fixed",
}


def detect_signals(sbol_data, ignore_component_ids=None):
    """
    Parse `sbol_data` and return the same signalling topology that
    `generate_script()` will auto-detect: which chemicals diffuse, which
    proteins produce each one, and which promoters they activate. Use this
    to know up front which signals need diffusion/degradation rate
    controls in a UI, or simply to sanity-check a circuit.
    """
    proteins, modules, interactions, component_map, ed_chemicals = parse_json(
        sbol_data, ignore_ids=ignore_component_ids
    )
    diffusible_signals, signal_producers, signal_activated_promoters = \
        detect_signaling_topology(proteins, interactions, component_map, ed_chemicals)

    return {
        "proteins":                   proteins,
        "modules":                    modules,
        "interactions":               interactions,
        "component_map":              component_map,
        "ed_chemicals":               ed_chemicals,
        "diffusible_signals":         diffusible_signals,
        "signal_producers":           signal_producers,
        "signal_activated_promoters": signal_activated_promoters,
    }


def _deep_merge_parameters(parameters, overrides):
    """Merge a user-supplied `overrides` dict into `parameters`, in place."""
    for key in ("simulation", "kinetics", "sbol_mapping"):
        if key in overrides:
            parameters[key].update(overrides[key])

    if "cell_types" in overrides:
        parameters["cell_types"] = overrides["cell_types"]

    if "signaling" in overrides:
        sig_over = overrides["signaling"]
        for key in ("enabled", "grid_len", "grid_size", "boundary_condition"):
            if key in sig_over:
                parameters["signaling"][key] = sig_over[key]
        if "signals" in sig_over:
            by_name = {s["name"]: s for s in parameters["signaling"]["signals"]}
            for s_over in sig_over["signals"]:
                existing = by_name.get(s_over["name"])
                if existing is not None:
                    existing.update(s_over)
                else:
                    parameters["signaling"]["signals"].append(dict(s_over))


def prepare_parameters_and_data(sbol_data, overrides=None):
    """
    Build the complete `parameters` dict consumed by
    `generate_script(proteins, modules, interactions, component_map,
    parameters, ed_chemicals)` / `main()` in cellmodeller_converter.py,
    starting from sensible defaults, auto-filling one signalling entry per
    diffusible signal detected in the SBOL circuit, and applying any
    `overrides` collected from a UI (see ui_params.display_form) or passed
    in directly.

    Args:
      sbol_data (dict): SBOL JSON circuit description (keys "ED",
                        "interactions", "hierarchy", "components").
      overrides (dict):  Partial parameters dict to merge on top of the
                        defaults / auto-detected signals. Recognised
                        top-level keys: "simulation", "cell_types",
                        "signaling" (with a "signals" list of
                        {"name", "diffusion_rate", "degradation_rate"}),
                        "kinetics", "sbol_mapping".

    Returns:
      tuple:
        - parameters (dict): ready to pass straight into
          `generate_script(...)`, or to `json.dump()` for the CLI.
        - topology (dict): the auto-detected signalling topology (see
          `detect_signals`) — handy for printing a summary or building a UI.
    """
    overrides = overrides or {}
    ignore_ids = overrides.get("sbol_mapping", {}).get("ignore_component_ids", [])
    topology = detect_signals(sbol_data, ignore_component_ids=ignore_ids)

    parameters = {
        "simulation":   dict(DEFAULT_SIMULATION),
        "cell_types":   [dict(DEFAULT_CELL_TYPE)],
        "signaling":    {**DEFAULT_SIGNALING, "signals": []},
        "kinetics":     dict(DEFAULT_KINETICS),
        "sbol_mapping": {"ignore_component_ids": list(ignore_ids)},
    }

    # One signal entry per auto-detected diffusible signal, so
    # generate_script() always has explicit rates for it (rather than
    # silently falling back to its own internal defaults).
    for sid in sorted(topology["diffusible_signals"]):
        parameters["signaling"]["signals"].append({
            "name": sid,
            **DEFAULT_SIGNAL_RATES,
        })

    # Only default signalling "on" if a diffusible signal actually exists
    # avoids implying a circuit has quorum sensing etc. when it doesn't.
    parameters["signaling"]["enabled"] = bool(topology["diffusible_signals"])

    # If a reporter protein is detected, suggest its colour for strain 0.
    for protein in topology["proteins"]:
        color_name = detect_reporter_color(protein["display_id"])
        if color_name and color_name in COLOR_NAME_TO_RGB:
            parameters["cell_types"][0]["color"] = COLOR_NAME_TO_RGB[color_name]
            parameters["cell_types"][0]["display_name"] = (
                f"Strain 0 ({protein['display_id']})"
            )
            break

    _deep_merge_parameters(parameters, overrides)

    return parameters, topology
