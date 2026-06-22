import json
import argparse
from datetime import datetime


# INTERACTION TYPES

PRODUCTION_TYPES  = {"Genetic Production", "Production"}
INHIBITION_TYPES  = {"Inhibition", "Repression", "Genetic Inhibition"}
STIMULATION_TYPES = {"Stimulation", "Activation", "Genetic Activation"}
DEGRADATION_TYPES = {"Degradation"}


# HELPERS

def _as_list(val):
    """Normalise a participants value (string or list) to a list, dropping empty strings."""
    if not val:
        return []
    return val if isinstance(val, list) else [val]


def safe_name(s):
    """Convert a display ID or name into a safe Python variable name."""
    return s.replace("-", "_").replace(".", "_").replace(" ", "_").lower()


# PARSER

def parse_json(sbol_data, ignore_ids=None):
    """
    Parse the SBOL JSON into five clean structures:
        proteins      — ED entries with type 'Protein'
        modules       — hierarchy entries
        interactions  — list of parsed interaction dicts
        component_map — displayId → role lookup
        ed_chemicals  — ED entries that are NOT proteins (small molecules, etc.)

    ignore_ids: optional set/list of component displayIds to exclude from
                component_map (populated from sbol_mapping.ignore_component_ids
                in the parameters file).

    Participants with the same role in one interaction (e.g. two "Inhibitor"
    entries) are stored as a list so no participant is silently overwritten.
    Single-occurrence roles stay as plain strings.
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

    # Non-protein ED entities (small molecules, inducers, etc.)
    # Treated as external concentration constants, not cell-state variables.
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
    for module_name, module_data in sbol_data.get("hierarchy", {}).items():
        modules.append({
            "name":         module_name,
            "components":   module_data.get("components", []),
            "constitutive": module_data.get("constitutive", False),
        })

    interactions = []
    for inter in sbol_data.get("interactions", []):
        itype = inter.get("type", "Unknown")
        parsed = {"type": itype, "participants": {}}
        for p in inter.get("participants", []):
            role        = p.get("role", "Unknown")
            participant = p.get("participant", "")
            existing    = parsed["participants"].get(role)
            if existing is None:
                # First occurrence — store as a plain string
                parsed["participants"][role] = participant
            elif isinstance(existing, list):
                # Already a list — append
                existing.append(participant)
            else:
                # Second occurrence of same role — promote to list
                parsed["participants"][role] = [existing, participant]
        interactions.append(parsed)

    return proteins, modules, interactions, component_map, ed_chemicals



# CIRCUIT ANALYSIS


def analyse_circuit(proteins, interactions, component_map):
    """
    Build two lookup tables that together describe the full regulatory chain:

        promoter_inhibitors: promoter_id → [inhibitor protein display_ids]
            e.g. "BBa_R0011" → ["LacI"]

        protein_regulation: protein_display_id → {
            "promoter":  promoter_id,
            "inhibitors": [protein display_ids that repress this protein's promoter],
            "activators": [protein display_ids that activate this protein's promoter],
        }

    
    Participants may be strings or lists (when multiple share the same role);
    _as_list() normalises both cases.
    """
    protein_ids = {p["display_id"] for p in proteins}

    # Step 1: map each promoter to the proteins it produces
    promoter_produces = {}  # promoter_id → [protein display_ids]
    protein_promoter  = {}  # protein display_id → promoter_id
    for inter in interactions:
        if inter["type"] in PRODUCTION_TYPES:
            promoter_id = inter["participants"].get("Promoter", "")
            product_id  = inter["participants"].get("Product", "")
            # Promoter and Product should always be single values in well-formed SBOL
            if isinstance(promoter_id, list): promoter_id = promoter_id[0]
            if isinstance(product_id,  list): product_id  = product_id[0]
            if promoter_id and product_id:
                promoter_produces.setdefault(promoter_id, []).append(product_id)
                protein_promoter[product_id] = promoter_id

    # Step 2: map each promoter to its inhibitors and activators
    promoter_inhibitors = {}  # promoter_id → [inhibitor protein display_ids]
    promoter_activators = {}  # promoter_id → [activator protein display_ids]
    for inter in interactions:
        if inter["type"] in INHIBITION_TYPES:
            inhibitor_ids = _as_list(inter["participants"].get("Inhibitor", ""))
            inhibited_id  = inter["participants"].get("Inhibited", "")
            if isinstance(inhibited_id, list): inhibited_id = inhibited_id[0]
            # Only resolve if the inhibited target is a DNA component (promoter),
            # not a protein — protein-level inhibition is handled separately.
            if inhibited_id in component_map:
                for inhibitor_id in inhibitor_ids:
                    if inhibitor_id in protein_ids:
                        promoter_inhibitors.setdefault(inhibited_id, []).append(inhibitor_id)

        elif inter["type"] in STIMULATION_TYPES:
            activator_ids = _as_list(inter["participants"].get("Activator",
                            inter["participants"].get("Stimulator", "")))
            activated_id  = inter["participants"].get("Activated",
                            inter["participants"].get("Stimulated", ""))
            if isinstance(activated_id, list): activated_id = activated_id[0]
            if activated_id in component_map:
                for activator_id in activator_ids:
                    if activator_id in protein_ids:
                        promoter_activators.setdefault(activated_id, []).append(activator_id)

    # Step 3: build per-protein regulation summary
    protein_regulation = {}
    for protein in proteins:
        pid        = protein["display_id"]
        promoter   = protein_promoter.get(pid, "")
        inhibitors = promoter_inhibitors.get(promoter, [])
        activators = promoter_activators.get(promoter, [])
        protein_regulation[pid] = {
            "promoter":   promoter,
            "inhibitors": inhibitors,
            "activators": activators,
        }

    # Step 4: identify direct protein-level inhibitions (e.g. aTc → TetR)
    direct_inhibitions = []
    for inter in interactions:
        if inter["type"] in INHIBITION_TYPES:
            inhibitor_ids = _as_list(inter["participants"].get("Inhibitor", ""))
            inhibited_id  = inter["participants"].get("Inhibited", "")
            if isinstance(inhibited_id, list): inhibited_id = inhibited_id[0]
            if inhibited_id in protein_ids:  # target IS a protein, not a promoter
                for inhibitor_id in inhibitor_ids:
                    direct_inhibitions.append((inhibitor_id, inhibited_id))

    return protein_regulation, promoter_inhibitors, direct_inhibitions


def find_controlling_module(protein_display_id, modules, interactions):
    """Return the module that contains the promoter driving this protein."""
    for inter in interactions:
        if inter["type"] not in PRODUCTION_TYPES:
            continue
        if inter["participants"].get("Product", "") != protein_display_id:
            continue
        promoter_id = inter["participants"].get("Promoter", "")
        for module in modules:
            if promoter_id in module["components"]:
                return module
    return None


# LOGIC GENERATOR

def generate_update_logic(proteins, modules, interactions, component_map, params,
                          ed_chemicals=None):
    """
    Generate the body of update() using the resolved regulatory chain.

    Non-protein ED entities (small molecules such as aTc) are referenced as
    module-level concentration constants (e.g. ATC_CONC) rather than cell
    attributes, because their concentration is set externally by the experimenter.
    """
    kin     = params.get("kinetics", {})
    prod_r  = kin.get("production_rate",     0.05)
    max_r   = kin.get("max_production_rate",  0.2)
    deg_r   = kin.get("degradation_rate",     0.01)
    hill_n  = kin.get("hill_coefficient",     2.0)
    act_thr = kin.get("activation_threshold", 0.5)
    rep_thr = kin.get("repression_threshold", 0.5)

    protein_names    = {p["display_id"]: p["var_name"] for p in proteins}
    ed_chemical_ids  = {c["display_id"] for c in (ed_chemicals or [])}

    protein_regulation, promoter_inhibitors, direct_inhibitions = \
        analyse_circuit(proteins, interactions, component_map)

    lines = []

    # Pre-compute one inhibition factor per repressed promoter
    # These are reused by every protein that shares the same promoter.
    repressed_promoters = sorted(promoter_inhibitors.keys())
    if repressed_promoters:
        lines.append("        # --- Inhibition factors (one per repressed promoter) ---")
        for promoter_id in repressed_promoters:
            inhibitors = promoter_inhibitors[promoter_id]
            factor_var = f"inh_{safe_name(promoter_id)}"
            if len(inhibitors) == 1:
                inh_var = protein_names.get(inhibitors[0], safe_name(inhibitors[0]))
                lines.append(
                    f"        {factor_var} = {rep_thr}**{hill_n} / "
                    f"({rep_thr}**{hill_n} + cell.{inh_var}**{hill_n})"
                    f"  # {inhibitors[0]} represses {promoter_id}"
                )
            else:
                # Multiple inhibitors: multiply their individual factors
                factor_parts = []
                for inh_id in inhibitors:
                    inh_var = protein_names.get(inh_id, safe_name(inh_id))
                    factor_parts.append(
                        f"({rep_thr}**{hill_n} / ({rep_thr}**{hill_n} + cell.{inh_var}**{hill_n}))"
                    )
                lines.append(
                    f"        {factor_var} = {' * '.join(factor_parts)}"
                    f"  # {', '.join(inhibitors)} repress {promoter_id}"
                )
        lines.append("")

    # Production + degradation for each protein
    lines.append("        # Protein production and degradation ")
    for protein in proteins:
        pid      = protein["display_id"]
        var      = protein["var_name"]
        reg      = protein_regulation.get(pid, {})
        promoter = reg.get("promoter", "")
        inhibitors = reg.get("inhibitors", [])
        activators = reg.get("activators", [])

        module         = find_controlling_module(pid, modules, interactions)
        is_constitutive = module["constitutive"] if module else False

        if inhibitors:
            factor_var = f"inh_{safe_name(promoter)}"
            lines.append(
                f"        # {pid}: {', '.join(inhibitors)} → represses {promoter} → reduces {pid}"
            )
            lines.append(
                f"        cell.{var} += {max_r} * {factor_var} - {deg_r} * cell.{var}"
            )
        elif activators:
            factor_var = f"act_{safe_name(promoter)}"
            act_var    = protein_names.get(activators[0], safe_name(activators[0]))
            lines.append(
                f"        # {pid}: {', '.join(activators)} → activates {promoter} → increases {pid}"
            )
            lines.append(
                f"        act_{safe_name(promoter)} = cell.{act_var}**{hill_n} / "
                f"({act_thr}**{hill_n} + cell.{act_var}**{hill_n})"
            )
            lines.append(
                f"        cell.{var} += {max_r} * {factor_var} - {deg_r} * cell.{var}"
            )
        elif is_constitutive:
            lines.append(f"        # {pid}: constitutive from {promoter}")
            lines.append(f"        cell.{var} += {prod_r} - {deg_r} * cell.{var}")
        else:
            lines.append(
                f"        # {pid}: from {promoter} — regulation not resolved, "
                f"check SBOL manually"
            )
            lines.append(f"        cell.{var} += {prod_r} - {deg_r} * cell.{var}")

        lines.append(f"        cell.{var} = max(0.0, cell.{var})")
        lines.append("")

    # Direct protein-level inhibitions
    # Inhibitor may be a tracked protein (cell.x) or an external chemical (MODULE_CONC).
    for inhibitor_id, inhibited_id in direct_inhibitions:
        tgt_var = protein_names.get(inhibited_id, safe_name(inhibited_id))
        if inhibitor_id in ed_chemical_ids:
            # External small molecule — reference the module-level constant
            const_name = f"{safe_name(inhibitor_id).upper()}_CONC"
            lines.append(
                f"        # Direct inhibition: {inhibitor_id} (external chemical) "
                f"inhibits {inhibited_id} — adjust {const_name} at module level"
            )
            lines.append(
                f"        cell.{tgt_var} *= {rep_thr}**{hill_n} / "
                f"({rep_thr}**{hill_n} + {const_name}**{hill_n})"
            )
        else:
            inh_var = protein_names.get(inhibitor_id, safe_name(inhibitor_id))
            lines.append(f"        # Direct inhibition: {inhibitor_id} directly inhibits {inhibited_id}")
            lines.append(
                f"        cell.{tgt_var} *= {rep_thr}**{hill_n} / "
                f"({rep_thr}**{hill_n} + cell.{inh_var}**{hill_n})"
            )
        lines.append("")

    if not lines:
        lines.append("        pass")

    return "\n".join(lines)



# COLOUR HELPER


def generate_color_update(proteins, params):
    """
    If GFP is in the protein list, colour cells by GFP level.
    Otherwise use the static cell type color.
    Returns a string to insert into update(), or None.
    """
    gfp = next((p for p in proteins if "gfp" in p["display_id"].lower()), None)
    if not gfp:
        return None
    var = gfp["var_name"]
    rep_thr = params.get("kinetics", {}).get("repression_threshold", 0.5)
    return (
        f"        # Colour cells by GFP level (bright green = high, dark = low)\n"
        f"        gfp_norm = min(1.0, cell.{var} / {rep_thr})\n"
        f"        cell.color = [1.0 - gfp_norm, gfp_norm, 0.2]"
    )



# SCRIPT ASSEMBLER


def generate_script(proteins, modules, interactions, component_map, params,
                    ed_chemicals=None):
    sim_p        = params.get("simulation", {})
    cell_types   = params.get("cell_types", [{}])
    sig          = params.get("signaling", {})

    max_cells    = sim_p.get("max_cells",    10000)
    jitter_z     = sim_p.get("jitter_z",     False)
    gamma        = sim_p.get("gamma",        100.0)
    pickle_steps = sim_p.get("pickle_steps", 50)
    random_seed  = sim_p.get("random_seed",  None)

    sig_enabled   = sig.get("enabled", False)
    grid_size     = sig.get("grid_size", 100)
    signals       = sig.get("signals", [])
    boundary_cond = sig.get("boundary_condition", "fixed")

    # Cell type dicts
    color_lines, len_lines, growth_lines, noise_lines = [], [], [], []
    for i, ct in enumerate(cell_types):
        c = ct.get("color", [1.0, 0.3, 0.3])
        color_lines.append(  f"    {i}: [{c[0]}, {c[1]}, {c[2]}],  # {ct.get('display_name', f'Strain {i}')}")
        len_lines.append(    f"    {i}: {ct.get('division_length', 3.0)},")
        growth_lines.append( f"    {i}: {ct.get('growth_rate', 1.0)},")
        noise_lines.append(  f"    {i}: {ct.get('division_noise', 0.5)},")

    colors_dict = "{\n" + "\n".join(color_lines) + "\n}"
    lens_dict   = "{\n" + "\n".join(len_lines)   + "\n}"
    growth_dict = "{\n" + "\n".join(growth_lines) + "\n}"
    noise_dict  = "{\n" + "\n".join(noise_lines)  + "\n}"

    # addCell calls
    add_cell_lines = []
    for i, ct in enumerate(cell_types):
        pos = ct.get("initial_pos", [0.0, i * 3.0, 0.0])
        d   = ct.get("initial_dir", [1.0, 0.0, 0.0])
        add_cell_lines.append(
            f"    sim.addCell(cellType={i}, "
            f"pos=({pos[0]}, {pos[1]}, {pos[2]}), "
            f"dir=({d[0]}, {d[1]}, {d[2]}))"
        )
    add_cells_str = "\n".join(add_cell_lines)

    # Chemics
    if sig_enabled and signals:
        diff_rates = [str(s.get("diffusion_rate",  0.1))  for s in signals]
        deg_rates  = [str(s.get("degradation_rate", 0.01)) for s in signals]
        chem_import1   = "from CellModeller.Signalling.GridDiffusion import GridDiffusion"
        chem_import2   = "from CellModeller.Integration.CLCrankNicIntegrator import CLCrankNicIntegrator"
        chem_init     = (
            f"    chem = Chemics(sim,\n"
            f"        nSignals={len(signals)},\n"
            f"        diffusion=[{', '.join(diff_rates)}],\n"
            f"        degradation=[{', '.join(deg_rates)}],\n"
            f"        gridSize={grid_size})  # boundary_condition='{boundary_cond}'\n"
        )
        sim_init_call = "    sim.init(biophys, regul, chem, None)"
    else:
        chem_import1   = "# from CellModeller.Signalling.GridDiffusion import GridDiffusion"
        chem_import2   = "# from CellModeller.Integration.CLCrankNicIntegrator import CLCrankNicIntegrator"
        chem_init     = "    # Signaling disabled — set 'enabled': true in parameters to activate\n"
        sim_init_call = "    sim.init(biophys, regul, None, None)"

    # random_seed line for setup()
    if random_seed is not None:
        random_seed_line = f"    random.seed({random_seed})  # set in simulation.random_seed"
    else:
        random_seed_line = (
            "    # random_seed not set — add 'random_seed' to simulation params for reproducibility"
        )

    # Protein init / divide
    if proteins:
        protein_init = "\n".join(
            f"    cell.{p['var_name']} = 0.0  # {p['display_id']}" for p in proteins
        )
        protein_divide = "\n".join(
            f"    d1.{p['var_name']} = parent.{p['var_name']} / 2.0\n"
            f"    d2.{p['var_name']} = parent.{p['var_name']} / 2.0"
            for p in proteins
        )
    else:
        protein_init   = "    # No proteins found in ED section"
        protein_divide = "    pass"

    # External chemical concentration constants
    if ed_chemicals:
        chem_consts = "\n".join(
            f"{safe_name(c['display_id']).upper()}_CONC = 0.0"
            f"  # {c['display_id']} ({c['type']}) — set to simulate induction"
            for c in ed_chemicals
        )
        chem_consts_section = (
            "\n# EXTERNAL CHEMICAL CONCENTRATIONS  (set before running)\n"
            "# These are non-protein ED entities from the SBOL (e.g. small molecule inducers).\n"
            + chem_consts + "\n"
        )
    else:
        chem_consts_section = ""

    # Update body
    update_logic = generate_update_logic(
        proteins, modules, interactions, component_map, params,
        ed_chemicals=ed_chemicals
    )

    color_update = generate_color_update(proteins, params)
    color_update_str = f"\n{color_update}\n" if color_update else ""

    
    script = f'''# CellModeller Script

from CellModeller.Regulation.ModuleRegulator import ModuleRegulator
from CellModeller.Biophysics.BacterialModels.CLBacterium import CLBacterium
{chem_import1}
{chem_import2}
import numpy as np
import random
{chem_consts_section}

# CELL TYPE PROPERTIES  (edit in cellmodeller_parameters_template.json)

cell_colors = {colors_dict}

cell_lens = {lens_dict}

cell_growth_rates = {growth_dict}

cell_division_noise = {noise_dict}



# SETUP

def setup(sim):
{random_seed_line}
    biophys = CLBacterium(
        sim,
        jitter_z={jitter_z},
        max_cells={max_cells},
        gamma={gamma}
    )
    regul = ModuleRegulator(sim, sim.moduleName)

{chem_init}
{sim_init_call}

    sim.pickleSteps = {pickle_steps}

{add_cells_str}



# INIT

def init(cell):
    cell.targetVol  = cell_lens[cell.cellType] + random.uniform(0.0, cell_division_noise[cell.cellType])
    cell.growthRate = cell_growth_rates[cell.cellType]
    cell.color      = cell_colors[cell.cellType]

    # Proteins tracked from SBOL ED section
{protein_init}


# UPDATE

def update(cells):
    for (id, cell) in cells.items():

        # Division
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

{protein_divide}
'''

    if sig_enabled and signals:
        sig_lines = "\n".join(
            f"        signals[cell_idx, {i}] = 0.1  # {s.get('name', f'signal_{i}')}"
            for i, s in enumerate(signals)
        )
        script += f'''


# SIGNAL RATES

def signalRates(cells, signals, pos, signalLevel):
    for cell_idx, (id, cell) in enumerate(cells.items()):
{sig_lines}
'''
    return script



# MAIN

def main():
    parser = argparse.ArgumentParser(
        description="Convert a pre-processed SBOL JSON to a CellModeller script."
    )
    parser.add_argument("--sbol",   required=True, help="Path to SBOL JSON file")
    parser.add_argument("--params", required=True, help="Path to parameters JSON file")
    parser.add_argument("--output", default="cellmodeller_output.py",
                        help="Output .py filename (default: cellmodeller_output.py)")
    args = parser.parse_args()

    print(f"Loading SBOL JSON:  {args.sbol}")
    with open(args.sbol)   as f: sbol_data = json.load(f)
    print(f"Loading parameters: {args.params}")
    with open(args.params) as f: params    = json.load(f)

    sbol_mapping = params.get("sbol_mapping", {})
    ignore_ids   = sbol_mapping.get("ignore_component_ids", [])
    if ignore_ids:
        print(f"  Ignoring components: {ignore_ids}")

    proteins, modules, interactions, component_map, ed_chemicals = \
        parse_json(sbol_data, ignore_ids=ignore_ids)

    print(f"  Proteins  (ED):   {[p['display_id'] for p in proteins]}")
    print(f"  Chemicals (ED):   {[c['display_id'] for c in ed_chemicals]}")
    print(f"  Modules:          {[m['name'] for m in modules]}")
    print(f"  Interactions:     {[i['type'] for i in interactions]}")

    script = generate_script(
        proteins, modules, interactions, component_map, params,
        ed_chemicals=ed_chemicals
    )

    with open(args.output, "w") as f:
        f.write(script)
    print(f"Script written to:  {args.output}")


if __name__ == "__main__":
    main()
