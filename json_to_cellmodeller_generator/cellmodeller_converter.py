"""
SBOL JSON to CellModeller Script Converter

Topology (who produces what, what activates what) is auto-detected from SBOL.
Numbers (rates, grid dimensions) come from the parameters dict.

CLI usage:
    python json_to_cellmodeller.py --sbol circuit.json --params params.json --output sim.py
"""

import json
import argparse
from datetime import datetime


# INTERACTION TYPE SETS 

PRODUCTION_TYPES       = {"Genetic Production", "Production"}
INHIBITION_TYPES       = {"Inhibition", "Repression", "Genetic Inhibition"}
STIMULATION_TYPES      = {"Stimulation", "Activation", "Genetic Activation"}
DEGRADATION_TYPES      = {"Degradation"}
BIOCHEM_REACTION_TYPES = {"Biochemical Reaction", "Non-Covalent Binding", "Control"}


# FLUORESCENT PROTEIN COLOUR LOOKUP
# Keyword (substring of display_id, lowercased) → (R, G, B)
_FP_KEYWORDS = {
    "gfp":        (0.0, 1.0, 0.2),
    "egfp":       (0.0, 1.0, 0.2),
    "mgfp":       (0.0, 1.0, 0.2),
    "mneongreen": (0.0, 1.0, 0.2),
    "rfp":        (1.0, 0.1, 0.1),
    "dsred":      (1.0, 0.1, 0.1),
    "mrfp":       (1.0, 0.1, 0.1),
    "mcherry":    (0.9, 0.1, 0.3),
    "mstrawberry":(0.9, 0.2, 0.1),
    "mkate":      (0.9, 0.1, 0.2),
    "mruby":      (0.85, 0.1, 0.1),
    "mneptune":   (0.7, 0.0, 0.3),
    "cfp":        (0.1, 0.8, 1.0),
    "ecfp":       (0.1, 0.8, 1.0),
    "mturquoise": (0.1, 0.9, 0.8),
    "mtfp":       (0.1, 0.9, 0.8),
    "yfp":        (0.9, 0.9, 0.0),
    "eyfp":       (0.9, 0.9, 0.0),
    "mvenus":     (0.9, 0.9, 0.0),
    "mcitrine":   (0.9, 0.9, 0.0),
    "morange":    (1.0, 0.5, 0.0),
}

# Common BioBrick registry IDs → (R, G, B)
_FP_BIOBRICK = {
    "BBa_E0040": (0.0, 1.0, 0.2),  # GFP
    "BBa_K2148009": (0.0, 1.0, 0.2),
    "BBa_K2560042": (0.0, 1.0, 0.2),
    "BBa_K4946001": (0.0, 1.0, 0.2),
    "BBa_E1010": (1.0, 0.1, 0.1),  # RFP
    "BBa_K1323009": (1.0, 0.1, 0.1),
    "BBa_K3128008": (1.0, 0.1, 0.1),
    "BBa_E0020": (0.1, 0.8, 1.0),  # CFP
    "BBa_E0030": (0.9, 0.9, 0.0),  # YFP
    "BBa_K592101": (0.9, 0.9, 0.0),
    "BBa_K3427000": (0.9, 0.9, 0.0),
}


# HELPERS 

def _as_list(val):
    """Normalise a participants value (string or list) to a list."""
    if not val:
        return []
    return val if isinstance(val, list) else [val]


def safe_name(s):
    """Convert any display-ID to a safe Python identifier."""
    name = (s.replace("-", "_")
             .replace(".", "_")
             .replace(" ", "_")
             .replace("(", "_")
             .replace(")", "_")
             .lower())
    # Identifiers can't start with a digit — common for real signal names
    # like "3OC6HSL" / "3OC12HSL", which would otherwise produce invalid
    # syntax such as `cell.3oc6hsl_sensed = 0.0`.
    if name and name[0].isdigit():
        name = "_" + name
    return name


def _first(val):
    """Return first element if list, else the value itself."""
    return val[0] if isinstance(val, list) else val


# PARSER

def parse_json(sbol_data, ignore_ids=None):
    """
    Parse SBOL JSON into five structures:

        proteins      — ED entries with type 'Protein'
        modules       — hierarchy entries with constitutive flag
        interactions  — normalised list (same-role participants → list)
        component_map — displayId → role (DNA parts etc.)
        ed_chemicals  — non-protein ED entities (small molecules, complexes …)
    """
    ignore_ids = set(ignore_ids or [])

    component_map = {
        c["displayId"]: c.get("role", "Unknown")
        for c in sbol_data.get("components", [])
        if c["displayId"] not in ignore_ids
    }

    proteins = [
        {
            "display_id": ed.get("name", ed.get("displayId", "")),
            "var_name":   safe_name(ed.get("name", ed.get("displayId", ""))),
        }
        for ed in sbol_data.get("ED", [])
        if ed.get("type", "").lower() == "protein"
    ]

    ed_chemicals = [
        {
            "display_id": ed.get("name", ed.get("displayId", "")),
            "var_name":   safe_name(ed.get("name", ed.get("displayId", ""))),
            "type":       ed.get("type", "Unknown"),
        }
        for ed in sbol_data.get("ED", [])
        if ed.get("type", "").lower() != "protein"
    ]

    modules = []
    for name, data in sbol_data.get("hierarchy", {}).items():
        modules.append({
            "name":         name,
            "components":   data.get("components", []),
            "constitutive": data.get("constitutive", False),
        })

    interactions = []
    for inter in sbol_data.get("interactions", []):
        itype  = inter.get("type", "Unknown")
        parsed = {"type": itype, "participants": {}}
        for p in inter.get("participants", []):
            role        = p.get("role", "Unknown")
            participant = p.get("participant", "")
            existing    = parsed["participants"].get(role)
            if existing is None:
                parsed["participants"][role] = participant
            elif isinstance(existing, list):
                existing.append(participant)
            else:
                parsed["participants"][role] = [existing, participant]
        interactions.append(parsed)

    return proteins, modules, interactions, component_map, ed_chemicals


# SIGNAL TOPOLOGY DETECTION

def detect_signaling_topology(proteins, interactions, component_map, ed_chemicals):
    """
    Auto-detect the three intercellular-signalling facts from SBOL topology.

    Returns
   
    diffusible_signals : set[str]
        display_ids of ED chemicals that appear as Stimulator in any
        Stimulation interaction.  These are the molecules that GridDiffusion
        will diffuse across the colony (e.g. AHL).

    signal_producers : dict[str, list[str]]
        signal_display_id → [protein_display_ids] that synthesise it.
        Detected from Biochemical Reaction / Non-Covalent Binding interactions
        where the product is a diffusible signal (e.g. LuxI → AHL).

    signal_activated_promoters : dict[str, str]
        promoter_displayId → signal_display_id.
        Promoters whose Stimulated target is driven by a chemical Stimulator
        (e.g. BBa_R0062 activated by AHL or LuxR-AHL complex).
    """
    protein_ids  = {p["display_id"] for p in proteins}
    ed_chem_ids  = {c["display_id"] for c in (ed_chemicals or [])}

    # Chemicals that appear as Stimulator
    diffusible_signals: set = set()
    for inter in interactions:
        if inter["type"] in STIMULATION_TYPES:
            for role in ("Stimulator", "Activator"):
                for sid in _as_list(inter["participants"].get(role, "")):
                    if sid in ed_chem_ids:
                        diffusible_signals.add(sid)

    # Proteins that produce a diffusible signal 
    signal_producers: dict = {}
    for inter in interactions:
        if inter["type"] in BIOCHEM_REACTION_TYPES:
            product_id = _first(inter["participants"].get("Product", ""))
            if not product_id or product_id not in diffusible_signals:
                continue
            # Accept Modifier or Reactant role as "producer protein"
            for role in ("Modifier", "Reactant", "Template"):
                for mid in _as_list(inter["participants"].get(role, "")):
                    if mid in protein_ids:
                        signal_producers.setdefault(product_id, []).append(mid)

    # Promoters activated by a diffusible signal
    signal_activated_promoters: dict = {}
    for inter in interactions:
        if inter["type"] in STIMULATION_TYPES:
            stimulated_id = _first(
                inter["participants"].get(
                    "Stimulated", inter["participants"].get("Activated", "")))
            if not stimulated_id or stimulated_id not in component_map:
                continue
            for role in ("Stimulator", "Activator"):
                for sid in _as_list(inter["participants"].get(role, "")):
                    if sid in diffusible_signals:
                        signal_activated_promoters[stimulated_id] = sid

    return diffusible_signals, signal_producers, signal_activated_promoters


# CIRCUIT ANALYSIS

def analyse_circuit(proteins, interactions, component_map, ed_chemicals=None):
    """
    Build complete regulatory lookup tables, including signal topology.

    Returns
    
    protein_regulation : dict
        display_id → {promoter, inhibitors, activators, signal_activator}
    promoter_inhibitors : dict
        promoter_id → [inhibitor protein display_ids]
    direct_inhibitions : list[(inhibitor_id, inhibited_protein_id)]
    diffusible_signals, signal_producers, signal_activated_promoters
        (see detect_signaling_topology)
    """
    protein_names = {p["display_id"]: p["var_name"] for p in proteins}
    protein_ids   = set(protein_names)
    ed_chem_ids   = {c["display_id"] for c in (ed_chemicals or [])}

    # Detect signal topology first
    diffusible_signals, signal_producers, signal_activated_promoters = \
        detect_signaling_topology(proteins, interactions, component_map, ed_chemicals)

    # promoter to protein (from Production interactions)
    protein_promoter: dict = {}
    for inter in interactions:
        if inter["type"] in PRODUCTION_TYPES:
            promoter_id = _first(inter["participants"].get("Promoter", ""))
            product_id  = _first(inter["participants"].get("Product",  ""))
            if promoter_id and product_id:
                protein_promoter[product_id] = promoter_id

    # promoter to inhibitors / activators 
    promoter_inhibitors: dict = {}
    promoter_activators: dict = {}

    for inter in interactions:
        if inter["type"] in INHIBITION_TYPES:
            inhibited_id = _first(inter["participants"].get("Inhibited", ""))
            if inhibited_id not in component_map:
                continue  # not a DNA part — handle as direct inhibition below
            for inh_id in _as_list(inter["participants"].get("Inhibitor", "")):
                if inh_id in protein_ids:
                    promoter_inhibitors.setdefault(inhibited_id, []).append(inh_id)

        elif inter["type"] in STIMULATION_TYPES:
            activated_id = _first(
                inter["participants"].get(
                    "Activated", inter["participants"].get("Stimulated", "")))
            if activated_id not in component_map:
                continue
            for act_id in _as_list(
                    inter["participants"].get(
                        "Activator", inter["participants"].get("Stimulator", ""))):
                # Only protein activators here; chemical ones are in signal_activated_promoters
                if act_id in protein_ids:
                    promoter_activators.setdefault(activated_id, []).append(act_id)

    # per-protein regulation summary 
    protein_regulation: dict = {}
    for protein in proteins:
        pid      = protein["display_id"]
        promoter = protein_promoter.get(pid, "")
        protein_regulation[pid] = {
            "promoter":          promoter,
            "inhibitors":        promoter_inhibitors.get(promoter, []),
            "activators":        promoter_activators.get(promoter, []),
            # NEW: diffusible signal that activates this protein's promoter (or None)
            "signal_activator":  signal_activated_promoters.get(promoter),
        }

    # direct protein-level inhibitions (e.g. aTc to TetR) 
    direct_inhibitions = []
    for inter in interactions:
        if inter["type"] in INHIBITION_TYPES:
            inhibited_id = _first(inter["participants"].get("Inhibited", ""))
            if inhibited_id in protein_ids:          # target is a protein, not a promoter
                for inh_id in _as_list(inter["participants"].get("Inhibitor", "")):
                    direct_inhibitions.append((inh_id, inhibited_id))

    return (protein_regulation, promoter_inhibitors, direct_inhibitions,
            diffusible_signals, signal_producers, signal_activated_promoters)


# MODULE HELPER

def find_controlling_module(protein_display_id, modules, interactions):
    """Return the module dict that contains the promoter driving this protein."""
    for inter in interactions:
        if inter["type"] not in PRODUCTION_TYPES:
            continue
        if _first(inter["participants"].get("Product", "")) != protein_display_id:
            continue
        promoter_id = _first(inter["participants"].get("Promoter", ""))
        for module in modules:
            if promoter_id in module["components"]:
                return module
    return None


# UPDATE LOGIC GENERATOR

def generate_update_logic(proteins, modules, interactions, component_map, params,
                          ed_chemicals=None):
    """
    Generate the body of update() — all lines are indented for the
    ``for id, cell in cells.items():`` loop (8 spaces).

    Signal-activated proteins read ``cell.<signal>_sensed``, which is written
    by signalRates() each timestep so there is only a 1-step delay.
    """
    kin     = params.get("kinetics", {})
    prod_r  = kin.get("production_rate",     0.05)
    max_r   = kin.get("max_production_rate",  0.2)
    deg_r   = kin.get("degradation_rate",     0.01)
    hill_n  = kin.get("hill_coefficient",     2.0)
    act_thr = kin.get("activation_threshold", 0.5)
    rep_thr = kin.get("repression_threshold", 0.5)

    protein_names = {p["display_id"]: p["var_name"] for p in proteins}
    ed_chem_ids   = {c["display_id"] for c in (ed_chemicals or [])}

    (protein_regulation, promoter_inhibitors, direct_inhibitions,
     diffusible_signals, _, _) = \
        analyse_circuit(proteins, interactions, component_map, ed_chemicals)

    lines = []

    # Pre-compute inhibition Hill factors (one per repressed promoter, reused by every protein sharing that promoter)
    repressed_promoters = sorted(promoter_inhibitors)
    if repressed_promoters:
        lines.append("        # — inhibition factors —")
        for prom in repressed_promoters:
            inhibitors = promoter_inhibitors[prom]
            fvar = f"_inh_{safe_name(prom)}"
            if len(inhibitors) == 1:
                iv = protein_names.get(inhibitors[0], safe_name(inhibitors[0]))
                lines.append(
                    f"        {fvar} = {rep_thr}**{hill_n} / "
                    f"({rep_thr}**{hill_n} + cell.{iv}**{hill_n})"
                    f"  # {inhibitors[0]} represses {prom}"
                )
            else:
                parts = [
                    f"({rep_thr}**{hill_n} / ({rep_thr}**{hill_n} + "
                    f"cell.{protein_names.get(iid, safe_name(iid))}**{hill_n}))"
                    for iid in inhibitors
                ]
                lines.append(
                    f"        {fvar} = {' * '.join(parts)}"
                    f"  # {', '.join(inhibitors)} repress {prom}"
                )
        lines.append("")

    # Per-protein production + degradation
    lines.append("        # — protein production and degradation —")
    for protein in proteins:
        pid         = protein["display_id"]
        var         = protein["var_name"]
        reg         = protein_regulation.get(pid, {})
        promoter    = reg.get("promoter", "")
        inhibitors  = reg.get("inhibitors", [])
        activators  = reg.get("activators", [])
        signal_act  = reg.get("signal_activator")   # e.g. "AHL", or None

        module          = find_controlling_module(pid, modules, interactions)
        is_constitutive = module["constitutive"] if module else False

        if inhibitors and signal_act:
            # Combined: signal activates AND protein represses
            fvar      = f"_inh_{safe_name(promoter)}"
            sensed    = f"cell.{safe_name(signal_act)}_sensed"
            act_local = f"_act_{safe_name(promoter)}"
            lines.append(
                f"        # {pid}: {signal_act} activates, "
                f"{', '.join(inhibitors)} represses via {promoter}"
            )
            lines.append(
                f"        {act_local} = {sensed}**{hill_n} / "
                f"({act_thr}**{hill_n} + {sensed}**{hill_n})"
            )
            lines.append(
                f"        cell.{var} += {max_r} * {act_local} * {fvar}"
                f" - {deg_r} * cell.{var}"
            )

        elif inhibitors:
            fvar = f"_inh_{safe_name(promoter)}"
            lines.append(
                f"        # {pid}: repressed by {', '.join(inhibitors)} via {promoter}"
            )
            lines.append(
                f"        cell.{var} += {max_r} * {fvar} - {deg_r} * cell.{var}"
            )

        elif signal_act:
            # Activated by diffusible signal — reads cell.<sig>_sensed from signalRates()
            sensed    = f"cell.{safe_name(signal_act)}_sensed"
            act_local = f"_act_{safe_name(promoter) or safe_name(pid)}"
            lines.append(
                f"        # {pid}: activated by diffusible signal {signal_act} via {promoter}"
            )
            lines.append(
                f"        {act_local} = {sensed}**{hill_n} / "
                f"({act_thr}**{hill_n} + {sensed}**{hill_n})"
            )
            lines.append(
                f"        cell.{var} += {max_r} * {act_local} - {deg_r} * cell.{var}"
            )

        elif activators:
            av        = protein_names.get(activators[0], safe_name(activators[0]))
            act_local = f"_act_{safe_name(promoter) or safe_name(pid)}"
            lines.append(
                f"        # {pid}: activated by {', '.join(activators)} via {promoter}"
            )
            lines.append(
                f"        {act_local} = cell.{av}**{hill_n} / "
                f"({act_thr}**{hill_n} + cell.{av}**{hill_n})"
            )
            lines.append(
                f"        cell.{var} += {max_r} * {act_local} - {deg_r} * cell.{var}"
            )

        elif is_constitutive:
            lines.append(f"        # {pid}: constitutive from {promoter}")
            lines.append(f"        cell.{var} += {prod_r} - {deg_r} * cell.{var}")

        else:
            lines.append(f"        # {pid}: constitutive fallback (promoter: {promoter or 'unresolved'})")
            lines.append(f"        cell.{var} += {prod_r} - {deg_r} * cell.{var}")

        lines.append(f"        cell.{var} = max(0.0, cell.{var})")
        lines.append("")

    # Direct protein-level inhibitions (e.g. external small molecule sequesters a TF)
    for inhibitor_id, inhibited_id in direct_inhibitions:
        tgt = protein_names.get(inhibited_id, safe_name(inhibited_id))
        if inhibitor_id in ed_chem_ids:
            cname = f"{safe_name(inhibitor_id).upper()}_CONC"
            lines.append(
                f"        # Direct inhibition: {inhibitor_id} (external) → {inhibited_id}"
            )
            lines.append(
                f"        cell.{tgt} *= {rep_thr}**{hill_n} / "
                f"({rep_thr}**{hill_n} + {cname}**{hill_n})"
            )
        else:
            iv = protein_names.get(inhibitor_id, safe_name(inhibitor_id))
            lines.append(
                f"        # Direct inhibition: {inhibitor_id} → {inhibited_id}"
            )
            lines.append(
                f"        cell.{tgt} *= {rep_thr}**{hill_n} / "
                f"({rep_thr}**{hill_n} + cell.{iv}**{hill_n})"
            )
        lines.append("")

    return "\n".join(lines) if lines else "        pass"


# COLOR UPDATE HELPER 

def generate_color_update(proteins, params):
    """
    Return 8-space-indented lines to colour cells by FP expression level,
    or None if no fluorescent protein is detected.
    Interpolates from dark-grey (no expression) to full FP colour (saturation).
    """
    rep_thr = params.get("kinetics", {}).get("repression_threshold", 0.5)

    for protein in proteins:
        pid = protein["display_id"]
        var = protein["var_name"]

        # Exact BioBrick ID match
        if pid in _FP_BIOBRICK:
            r, g, b = _FP_BIOBRICK[pid]
        else:
            # Keyword substring match (case-insensitive)
            pid_lower = pid.lower()
            matched = next(
                (col for kw, col in _FP_KEYWORDS.items() if kw in pid_lower), None
            )
            if matched is None:
                continue
            r, g, b = matched

        return (
            f"        # colour by {pid} expression\n"
            f"        _fp_norm = min(1.0, cell.{var} / max({rep_thr}, 1e-9))\n"
            f"        cell.color = [\n"
            f"            {r} * _fp_norm + 0.1 * (1.0 - _fp_norm),\n"
            f"            {g} * _fp_norm + 0.1 * (1.0 - _fp_norm),\n"
            f"            {b} * _fp_norm + 0.1 * (1.0 - _fp_norm),\n"
            f"        ]"
        )
    return None


# SCRIPT ASSEMBLER

def generate_script(proteins, modules, interactions, component_map, params,
                    ed_chemicals=None):
    """
    Assemble a complete, runnable CellModeller Python script.

    Signalling topology (which signal diffuses, who emits it, what it activates)
    is auto-detected from SBOL.  All numerical rates are taken from *params*.
    """
    sim_p      = params.get("simulation",  {})
    cell_types = params.get("cell_types",  [{}])
    sig_p      = params.get("signaling",   {})
    kin        = params.get("kinetics",    {})

    max_cells    = sim_p.get("max_cells",    10000)
    jitter_z     = sim_p.get("jitter_z",     False)
    gamma        = sim_p.get("gamma",        100.0)
    pickle_steps = sim_p.get("pickle_steps", 50)
    random_seed  = sim_p.get("random_seed",  None)

    sig_enabled   = sig_p.get("enabled", False)
    # grid_len  = number of grid cells per axis
    # grid_size = physical size of each grid cell (µm); domain = grid_len × grid_size µm
    grid_len      = sig_p.get("grid_len",   100)
    grid_size     = sig_p.get("grid_size",  4.0)
    boundary_cond = sig_p.get("boundary_condition", "fixed")
    boundary_type = 1 if boundary_cond == "fixed" else 0

    signal_prod_rate = kin.get("signal_production_rate", 0.1)

    # Resolve signals
    (_, _, _,
     diffusible_signals, signal_producers, signal_activated_promoters) = \
        analyse_circuit(proteins, interactions, component_map, ed_chemicals)

    # Merge SBOL-detected signals with any explicitly listed in params.
    # Params-listed signals override defaults for diffusion / degradation rates.
    param_signal_map = {s["name"]: s for s in sig_p.get("signals", [])}

    all_signals = []
    for sid in sorted(diffusible_signals):
        entry = param_signal_map.get(sid, {})
        all_signals.append({
            "name":             sid,
            "diffusion_rate":   entry.get("diffusion_rate",  0.1),
            "degradation_rate": entry.get("degradation_rate", 0.01),
        })
    # Any param-listed signals not found in SBOL topology are included manually
    for name, s in param_signal_map.items():
        if name not in diffusible_signals:
            all_signals.append(s)

    signal_index = {s["name"]: i for i, s in enumerate(all_signals)}
    n_signals    = len(all_signals)
    use_signals  = sig_enabled and n_signals > 0

    # Cell type lookup dicts 
    color_lines, len_lines, growth_lines, noise_lines = [], [], [], []
    for i, ct in enumerate(cell_types):
        c = ct.get("color", [1.0, 0.3, 0.3])
        label = ct.get("display_name", f"Strain {i}")
        color_lines.append( f"    {i}: [{c[0]}, {c[1]}, {c[2]}],  # {label}")
        len_lines.append(   f"    {i}: {ct.get('division_length', 3.5)},")
        growth_lines.append(f"    {i}: {ct.get('growth_rate', 1.0)},")
        noise_lines.append( f"    {i}: {ct.get('division_noise', 0.5)},")

    colors_dict = "{\n" + "\n".join(color_lines) + "\n}"
    lens_dict   = "{\n" + "\n".join(len_lines)   + "\n}"
    growth_dict = "{\n" + "\n".join(growth_lines) + "\n}"
    noise_dict  = "{\n" + "\n".join(noise_lines)  + "\n}"

    # addCell calls 
    add_cell_lines = []
    for i, ct in enumerate(cell_types):
        pos = ct.get("initial_pos", [0.0, float(i) * 3.0, 0.0])
        d   = ct.get("initial_dir", [1.0, 0.0, 0.0])
        add_cell_lines.append(
            f"    sim.addCell(cellType={i}, "
            f"pos=({pos[0]}, {pos[1]}, {pos[2]}), "
            f"dir=({d[0]}, {d[1]}, {d[2]}))"
        )
    add_cells_str = "\n".join(add_cell_lines)

    # Signalling setup block
    if use_signals:
        diff_rates = [str(s["diffusion_rate"])  for s in all_signals]
        deg_rates  = [str(s["degradation_rate"]) for s in all_signals]
        sig_import1 = "from CellModeller.Signalling.GridDiffusion import GridDiffusion"
        sig_import2 = "from CellModeller.Integration.CLCrankNicIntegrator import CLCrankNicIntegrator"
        sig_init = (
            f"    nSignals = {n_signals}  # {[s['name'] for s in all_signals]}\n"
            f"    sig   = GridDiffusion(sim, nSignals,\n"
            f"                [gridLen, gridLen, 1],\n"
            f"                [gridSize, gridSize, 1.0],\n"
            f"                [{', '.join(diff_rates)}],   # diffusion coefficients\n"
            f"                [{', '.join(deg_rates)}])  # degradation rates\n"
            f"    integ = CLCrankNicIntegrator(sim, nSignals, nSignals, maxCells,\n"
            f"                sig, boundaryType={boundary_type})  # {boundary_cond} boundary"
        )
        sim_init_call = "    sim.init(biophys, regul, sig, integ)"
    else:
        sig_import1   = "# from CellModeller.Signalling.GridDiffusion import GridDiffusion"
        sig_import2   = "# from CellModeller.Integration.CLCrankNicIntegrator import CLCrankNicIntegrator"
        sig_init      = "    # Signalling disabled — set signaling.enabled=true in params to activate"
        sim_init_call = "    sim.init(biophys, regul, None, None)"

    # random seed 
    random_seed_line = (
        f"    random.seed({random_seed})"
        if random_seed is not None
        else "    # tip: set simulation.random_seed in params for reproducibility"
    )

    # init() protein attributes
    if proteins:
        init_proteins = "    # proteins\n" + "\n".join(
            f"    cell.{p['var_name']} = 0.0  # {p['display_id']}"
            for p in proteins
        )
        divide_proteins = "\n".join(
            f"    d1.{p['var_name']} = parent.{p['var_name']} / 2.0\n"
            f"    d2.{p['var_name']} = parent.{p['var_name']} / 2.0"
            for p in proteins
        )
    else:
        init_proteins   = "    # no proteins in SBOL ED"
        divide_proteins = ""

    # init() signal-sensing attributes
    # cell.<signal>_sensed bridges signalRates() → update()
    if use_signals and diffusible_signals:
        sig_sensed_attrs = sorted(
            f"{safe_name(s)}_sensed" for s in diffusible_signals
        )
        init_signals = "\n    # signal sensing (written by signalRates, read by update)\n" + \
            "\n".join(f"    cell.{v} = 0.0" for v in sig_sensed_attrs)
        divide_signals = "\n" + "\n".join(
            f"    d1.{v} = 0.0\n    d2.{v} = 0.0" for v in sig_sensed_attrs
        )
    else:
        init_signals   = ""
        divide_signals = ""

    # External (non-diffusible) chemical constants
    external_chems = [
        c for c in (ed_chemicals or [])
        if c["display_id"] not in diffusible_signals
    ]
    if external_chems:
        chem_consts = (
            "\n# EXTERNAL CHEMICAL CONCENTRATIONS"
            "  (non-diffusible inducers — set before running)\n"
            + "\n".join(
                f"{safe_name(c['display_id']).upper()}_CONC = 0.0"
                f"  # {c['display_id']} ({c['type']})"
                for c in external_chems
            ) + "\n"
        )
    else:
        chem_consts = ""

    # Topology comment block
    if diffusible_signals or signal_producers or signal_activated_promoters:
        topo_lines = ["# AUTO-DETECTED SIGNAL TOPOLOGY"]
        if diffusible_signals:
            topo_lines.append(f"#   Diffusible signals      : {sorted(diffusible_signals)}")
        for sig_id, prods in signal_producers.items():
            topo_lines.append(f"#   {sig_id:<22} produced by : {prods}")
        for prom, sig_id in signal_activated_promoters.items():
            topo_lines.append(f"#   {prom:<22} activated by: {sig_id}")
        topo_comment = "\n".join(topo_lines) + "\n"
    else:
        topo_comment = ""

    # update() body
    update_logic = generate_update_logic(
        proteins, modules, interactions, component_map, params, ed_chemicals
    )
    color_update_raw = generate_color_update(proteins, params)
    color_update_str = ("\n" + color_update_raw + "\n") if color_update_raw else ""

    # signalRates() body
    if use_signals:
        sr_lines = []
        for sig_obj in all_signals:
            sid      = sig_obj["name"]
            idx      = signal_index[sid]
            sensed   = f"{safe_name(sid)}_sensed"
            prods    = signal_producers.get(sid, [])

            if prods:
                prod_expr = " + ".join(f"cell.{safe_name(p)}" for p in prods)
                sr_lines.append(
                    f"        # {sid} (index {idx}) — emitted proportional to "
                    f"{', '.join(prods)}"
                )
                sr_lines.append(
                    f"        signals[id, {idx}] = {signal_prod_rate} * ({prod_expr})"
                )
            else:
                sr_lines.append(
                    f"        # {sid} (index {idx}) — no producer auto-detected; "
                    f"set manually if needed"
                )
                sr_lines.append(f"        signals[id, {idx}] = 0.0")

            # Write current local concentration into cell for update() to read
            sr_lines.append(
                f"        cell.{sensed} = signalLevel[id, {idx}]"
            )
            sr_lines.append("")

        sr_body       = "\n".join(sr_lines)
        signal_rates_fn = (
            "\n\n# SIGNAL RATES\n"
            "# signals[id, i]      — emission rate of signal i by cell id (set here)\n"
            "# signalLevel[id, i]  — current local concentration at cell id's position\n\n"
            "def signalRates(cells, signals, pos, signalLevel):\n"
            "    for id, cell in cells.items():\n"
            f"{sr_body}"
        )
    else:
        signal_rates_fn = ""

    # Assemble
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    script = f"""\
# CellModeller simulation script
# Generated : {now}
# Converter : json_to_cellmodeller.py
# Topology (who produces what signal, what it activates) is auto-detected from
# the SBOL JSON.  All numerical rates are defined in the params file / dict.


from CellModeller.Regulation.ModuleRegulator import ModuleRegulator
from CellModeller.Biophysics.BacterialModels.CLBacterium import CLBacterium
{sig_import1}
{sig_import2}
import numpy as np
import random
{chem_consts}
{topo_comment}
# simulation constants
maxCells = {max_cells}
gridLen  = {grid_len}   # grid cells per axis
gridSize = {grid_size}  # µm per grid cell  →  domain = {grid_len * grid_size:.0f} × {grid_len * grid_size:.0f} µm


# cell type lookup tables

cell_colors          = {colors_dict}
cell_lens            = {lens_dict}
cell_growth_rates    = {growth_dict}
cell_division_noise  = {noise_dict}


# SETUP

def setup(sim):
{random_seed_line}
    biophys = CLBacterium(
        sim,
        jitter_z={jitter_z},
        max_cells=maxCells,
        gamma={gamma}
    )
    regul = ModuleRegulator(sim, sim.moduleName)

{sig_init}

{sim_init_call}
    sim.pickleSteps = {pickle_steps}

{add_cells_str}


# INIT

def init(cell):
    cell.targetVol  = (cell_lens[cell.cellType]
                       + random.uniform(0.0, cell_division_noise[cell.cellType]))
    cell.growthRate = cell_growth_rates[cell.cellType]
    cell.color      = cell_colors[cell.cellType]
{init_proteins}{init_signals}


# UPDATE

def update(cells):
    for id, cell in cells.items():
        if cell.volume > cell.targetVol:
            cell.divideFlag = True
{color_update_str}
{update_logic}


# DIVIDE

def divide(parent, d1, d2):
    d1.cellType = parent.cellType
    d2.cellType = parent.cellType
    d1.color    = parent.color
    d2.color    = parent.color
{divide_proteins}{divide_signals}
{signal_rates_fn}
"""
    return script


# CLI ENTRY POINT

def main():
    parser = argparse.ArgumentParser(
        description="Convert an SBOL JSON circuit to a CellModeller simulation script."
    )
    parser.add_argument("--sbol",   required=True, help="Path to SBOL JSON")
    parser.add_argument("--params", required=True, help="Path to parameters JSON")
    parser.add_argument("--output", default="cellmodeller_output.py",
                        help="Output .py file (default: cellmodeller_output.py)")
    args = parser.parse_args()

    print(f"Loading SBOL   : {args.sbol}")
    with open(args.sbol)   as f: sbol_data = json.load(f)
    print(f"Loading params : {args.params}")
    with open(args.params) as f: params    = json.load(f)

    ignore_ids = params.get("sbol_mapping", {}).get("ignore_component_ids", [])
    proteins, modules, interactions, component_map, ed_chemicals = \
        parse_json(sbol_data, ignore_ids=ignore_ids)

    (_, _, _, diffusible_signals, signal_producers, signal_activated_promoters) = \
        analyse_circuit(proteins, interactions, component_map, ed_chemicals)

    print(f"  Proteins               : {[p['display_id'] for p in proteins]}")
    print(f"  Chemicals (ED)         : {[c['display_id'] for c in ed_chemicals]}")
    print(f"  Diffusible signals     : {sorted(diffusible_signals)}")
    print(f"  Signal producers       : {signal_producers}")
    print(f"  Signal-active promoters: {signal_activated_promoters}")

    script = generate_script(
        proteins, modules, interactions, component_map, params, ed_chemicals
    )

    with open(args.output, "w") as f:
        f.write(script)
    print(f"Written: {args.output}")


if __name__ == "__main__":
    main()
