"""
sbml_to_cellmodeller.py

Direct SBML -> CellModeller script converter.

Why direct, and not through the SBOL-JSON bridge (cellmodeller_converter.py)?
SBOL only records *qualitative* interaction types ("X inhibits Y",
"X stimulates Y"), so cellmodeller_converter.py has to reconstruct rate
equations from scratch using generic Hill-function heuristics driven by a
handful of global "kinetics" parameters. SBML reactions already carry an
exact, per-reaction kinetic law (a MathML rate expression), so forcing an
SBML model through the same qualitative bridge would throw that precision
away and require re-fitting Hill parameters to match it. Instead, this
module:

  1. Parses the SBML file directly (compartments, species, parameters,
     reactions + <kineticLaw> MathML) with the standard library only.
  2. Translates every kinetic law's MathML into both a Python expression
     (for the plain-Python `update()` path) and an OpenCL C expression
     (for the `specRateCL()` / `sigRateCL()` GPU kernels used once
     inter-cell diffusion is enabled), using the *exact* rate law from the
     model.
  3. Sums each species' reaction stoichiometries into a per-species ODE,
     exactly as SBML defines dS/dt = sum_r (stoichiometry_r * rate_r).
  4. Classifies species as either "local" (tracked per cell only) or
     "diffusible" (also exchanged with a CellModeller GridDiffusion field)
     based on which SBML compartment they live in, and assembles a full
     CellModeller simulation script in the same style/shape as
     cellmodeller_converter.generate_script().

CLI usage:
    python sbml_to_cellmodeller.py --sbml circuit.xml --params params.json --output sim.py
"""

import re
import json
import math
import argparse
import xml.etree.ElementTree as ET
from datetime import datetime

# Reuse naming / colour-detection / OpenCL-boundary helpers from the
# existing SBOL-JSON converter instead of duplicating them.
from cellmodeller_converter import safe_name, _FP_KEYWORDS, _FP_BIOBRICK, _BOUNDCOND_MAP


MATHML_NS = "http://www.w3.org/1998/Math/MathML"

# Compartment-name keywords that mark a compartment as "outside the cell"
# (i.e. its species should be exchanged with a CellModeller GridDiffusion
# field rather than only tracked per-cell). Override via
# params["sbml_mapping"]["diffusible_compartments"] / ["local_compartments"].
_EXTRACELLULAR_KEYWORDS = (
    "extracellular", "medium", "media", "environment", "external",
    "supernatant", "exterior", "outside", "out",
)


# XML / MathML HELPERS

def _local(tag):
    """Strip a `{namespace}tag` down to `tag`."""
    return tag.split("}", 1)[1] if "}" in tag else tag


def _findall_local(elem, name):
    return [c for c in elem if _local(c.tag) == name]


def _find_local(elem, name):
    for c in elem:
        if _local(c.tag) == name:
            return c
    return None


def _fmt_num(value, lang):
    """Format a literal float for Python or OpenCL C source."""
    s = repr(float(value))
    return f"{s}f" if lang == "c" else s


def _zero(lang):
    return "0.0f" if lang == "c" else "0.0"


# MATHML -> PYTHON / C EXPRESSION TRANSLATION

class MathTranslator:
    """
    Translates a MathML <math> subtree into either a Python or an OpenCL C
    expression string, resolving <ci> identifiers against a per-reaction
    symbol table (species -> cell.var / species[idx], parameters and
    compartment sizes -> inlined numeric literals).
    """

    def __init__(self, species_var, species_index, symbol_values, warnings, context):
        self.species_var    = species_var     # id -> python var name (for lang='py')
        self.species_index  = species_index   # id -> species[] slot   (for lang='c')
        self.symbol_values  = symbol_values    # id -> float (params + compartment sizes)
        self.warnings       = warnings
        self.context        = context

    def translate(self, math_el, lang):
        return self._walk(math_el, lang)

    # -- identifier resolution --
    def _resolve_ci(self, name, lang):
        if name in self.species_var:
            if lang == "py":
                return f"cell.{self.species_var[name]}"
            return f"species[{self.species_index[name]}]"
        if name in self.symbol_values:
            return _fmt_num(self.symbol_values[name], lang)
        self.warnings.append(
            f"{self.context}: unresolved identifier '{name}', substituted 0"
        )
        return _zero(lang)

    # -- cn numeric literal --
    @staticmethod
    def _cn_float(el):
        ctype = el.attrib.get("type", "real")
        if ctype == "integer":
            return float((el.text or "0").strip())
        if ctype == "e-notation":
            sep = _find_local(el, "sep")
            mantissa = float((el.text or "0").strip())
            exponent = float((sep.tail or "0").strip()) if sep is not None else 0.0
            return mantissa * (10.0 ** exponent)
        if ctype == "rational":
            sep = _find_local(el, "sep")
            numer = float((el.text or "0").strip())
            denom = float((sep.tail or "1").strip()) if sep is not None else 1.0
            return numer / denom if denom else 0.0
        return float((el.text or "0").strip())

    def _walk(self, el, lang):
        tag = _local(el.tag)

        if tag == "math":
            children = list(el)
            return self._walk(children[0], lang) if children else _zero(lang)

        if tag == "cn":
            return _fmt_num(self._cn_float(el), lang)

        if tag == "ci":
            return self._resolve_ci((el.text or "").strip(), lang)

        if tag == "csymbol":
            defn = el.attrib.get("definitionURL", "")
            if "time" in defn:
                self.warnings.append(
                    f"{self.context}: explicit time-dependence is not supported "
                    f"by CellModeller's per-cell update loop; substituted 0"
                )
            else:
                self.warnings.append(
                    f"{self.context}: unsupported csymbol '{defn}', substituted 0"
                )
            return _zero(lang)

        if tag == "true":
            return "True" if lang == "py" else "1"
        if tag == "false":
            return "False" if lang == "py" else "0"
        if tag == "pi":
            return _fmt_num(math.pi, lang)
        if tag == "exponentiale":
            return _fmt_num(math.e, lang)
        if tag == "infinity":
            return _fmt_num(1e30, lang)
        if tag == "notanumber":
            return _fmt_num(0.0, lang)

        if tag == "piecewise":
            return self._walk_piecewise(el, lang)

        if tag == "apply":
            kids = list(el)
            if not kids:
                return _zero(lang)
            op = _local(kids[0].tag)
            return self._walk_apply(op, kids[1:], lang)

        self.warnings.append(f"{self.context}: unsupported MathML tag <{tag}>, substituted 0")
        return _zero(lang)

    def _walk_piecewise(self, el, lang):
        pieces, otherwise = [], None
        for child in el:
            ctag = _local(child.tag)
            if ctag == "piece":
                sub = list(child)
                if len(sub) >= 2:
                    pieces.append((self._walk(sub[1], lang), self._walk(sub[0], lang)))
            elif ctag == "otherwise":
                sub = list(child)
                if sub:
                    otherwise = self._walk(sub[0], lang)
        expr = otherwise if otherwise is not None else _zero(lang)
        for cond, val in reversed(pieces):
            expr = f"({val} if {cond} else {expr})" if lang == "py" else f"({cond} ? {val} : {expr})"
        return expr

    def _walk_apply(self, op, arg_els, lang):
        def W(i):
            return self._walk(arg_els[i], lang)

        if op in ("plus", "times"):
            vals = [self._walk(a, lang) for a in arg_els]
            if not vals:
                return _fmt_num(0.0 if op == "plus" else 1.0, lang)
            joiner = " + " if op == "plus" else " * "
            return "(" + joiner.join(vals) + ")"

        if op == "minus":
            if len(arg_els) == 1:
                return f"(-({W(0)}))"
            return f"({W(0)} - {W(1)})"

        if op == "divide":
            return f"({W(0)} / ({W(1)}))"

        if op == "power":
            base, exp = W(0), W(1)
            return f"(({base}) ** ({exp}))" if lang == "py" else f"pow((float)({base}), (float)({exp}))"

        if op == "root":
            if len(arg_els) == 2 and _local(arg_els[0].tag) == "degree":
                deg_kids = list(arg_els[0])
                degree = self._walk(deg_kids[0], lang) if deg_kids else _fmt_num(2.0, lang)
                base = W(1)
            else:
                degree = _fmt_num(2.0, lang)
                base = W(0)
            if lang == "py":
                return f"(({base}) ** (1.0/({degree})))"
            return f"pow((float)({base}), 1.0f/({degree}))"

        if op == "abs":
            return f"abs({W(0)})" if lang == "py" else f"fabs({W(0)})"
        if op == "exp":
            return f"math.exp({W(0)})" if lang == "py" else f"exp({W(0)})"
        if op == "ln":
            return f"math.log({W(0)})" if lang == "py" else f"log({W(0)})"
        if op == "log":
            if arg_els and _local(arg_els[0].tag) == "logbase":
                base_kids = list(arg_els[0])
                base = self._walk(base_kids[0], lang) if base_kids else _fmt_num(10.0, lang)
                x = W(1)
                return f"(math.log({x}) / math.log({base}))" if lang == "py" else f"(log({x}) / log({base}))"
            x = W(0)
            return f"math.log10({x})" if lang == "py" else f"log10({x})"
        if op == "floor":
            return f"math.floor({W(0)})" if lang == "py" else f"floor({W(0)})"
        if op == "ceiling":
            return f"math.ceil({W(0)})" if lang == "py" else f"ceil({W(0)})"
        if op == "factorial":
            return f"math.factorial(int({W(0)}))" if lang == "py" else W(0)

        for trig, cname in (("sin", "sin"), ("cos", "cos"), ("tan", "tan"),
                            ("sinh", "sinh"), ("cosh", "cosh"), ("tanh", "tanh")):
            if op == trig:
                return f"math.{trig}({W(0)})" if lang == "py" else f"{cname}({W(0)})"

        if op in ("min", "max"):
            vals = [self._walk(a, lang) for a in arg_els]
            if not vals:
                return _zero(lang)
            if lang == "py":
                return f"{op}([{', '.join(vals)}])"
            expr = vals[0]
            for v in vals[1:]:
                expr = f"{op}({expr}, {v})"
            return expr

        if op in ("lt", "gt", "leq", "geq", "eq", "neq"):
            sym = {"lt": "<", "gt": ">", "leq": "<=", "geq": ">=", "eq": "==", "neq": "!="}[op]
            vals = [self._walk(a, lang) for a in arg_els]
            parts = [f"({vals[i]} {sym} {vals[i + 1]})" for i in range(len(vals) - 1)]
            joiner = " and " if lang == "py" else " && "
            return "(" + joiner.join(parts) + ")" if parts else ("True" if lang == "py" else "1")

        if op == "and":
            vals = [self._walk(a, lang) for a in arg_els]
            joiner = " and " if lang == "py" else " && "
            return "(" + joiner.join(vals) + ")"
        if op == "or":
            vals = [self._walk(a, lang) for a in arg_els]
            joiner = " or " if lang == "py" else " || "
            return "(" + joiner.join(vals) + ")"
        if op == "not":
            return f"(not ({W(0)}))" if lang == "py" else f"(!({W(0)}))"

        self.warnings.append(f"{self.context}: unsupported MathML operator '{op}', substituted 0")
        return _zero(lang)


# SBML PARSING

class SBMLModel:
    def __init__(self):
        self.compartments = {}   # id -> {"size": float, "name": str}
        self.species = {}        # id -> {"name","compartment","initial","boundary","constant"}
        self.parameters = {}     # id -> float  (global)
        self.reactions = []      # list of dicts (see parse_reactions)
        self.warnings = []


def parse_sbml(path_or_string, is_path=True):
    """Parse an SBML L2/L3 file (any level/version namespace) into an SBMLModel."""
    if is_path:
        tree = ET.parse(path_or_string)
        root = tree.getroot()
    else:
        root = ET.fromstring(path_or_string)

    model_el = _find_local(root, "model")
    if model_el is None:
        raise ValueError("No <model> element found — is this a valid SBML file?")

    m = SBMLModel()

    loc = _find_local(model_el, "listOfCompartments")
    if loc is not None:
        for c in _findall_local(loc, "compartment"):
            cid = c.attrib.get("id", "")
            size = c.attrib.get("size", c.attrib.get("volume", "1.0"))
            try:
                size_f = float(size)
            except ValueError:
                size_f = 1.0
            m.compartments[cid] = {"size": size_f, "name": c.attrib.get("name", cid)}

    los = _find_local(model_el, "listOfSpecies")
    if los is not None:
        for s in _findall_local(los, "species"):
            sid = s.attrib.get("id", "")
            init = s.attrib.get("initialConcentration", s.attrib.get("initialAmount", "0.0"))
            try:
                init_f = float(init)
            except ValueError:
                init_f = 0.0
            m.species[sid] = {
                "name":        s.attrib.get("name", sid),
                "compartment": s.attrib.get("compartment", ""),
                "initial":     init_f,
                "boundary":    s.attrib.get("boundaryCondition", "false").lower() == "true",
                "constant":    s.attrib.get("constant", "false").lower() == "true",
                "sbo":         s.attrib.get("sboTerm", ""),
            }

    lop = _find_local(model_el, "listOfParameters")
    if lop is not None:
        for p in _findall_local(lop, "parameter"):
            pid = p.attrib.get("id", "")
            val = p.attrib.get("value")
            if val is not None:
                try:
                    m.parameters[pid] = float(val)
                except ValueError:
                    pass

    lor = _find_local(model_el, "listOfReactions")
    if lor is not None:
        for r in _findall_local(lor, "reaction"):
            m.reactions.append(_parse_reaction(r, m))

    return m


def _species_refs(list_el):
    refs = []
    if list_el is None:
        return refs
    for sr in list_el:
        tag = _local(sr.tag)
        if tag in ("speciesReference", "modifierSpeciesReference"):
            sp = sr.attrib.get("species", "")
            stoich = sr.attrib.get("stoichiometry", "1")
            try:
                stoich_f = float(stoich)
            except ValueError:
                stoich_f = 1.0
            refs.append({"species": sp, "stoichiometry": stoich_f})
    return refs


def _parse_reaction(r, model):
    rid = r.attrib.get("id", "unnamed_reaction")
    reactants = _species_refs(_find_local(r, "listOfReactants"))
    products  = _species_refs(_find_local(r, "listOfProducts"))
    modifiers = [ref["species"] for ref in _species_refs(_find_local(r, "listOfModifiers"))]

    local_params = {}
    kinetic_math = None
    kl = _find_local(r, "kineticLaw")
    if kl is not None:
        # SBML L2 <listOfParameters>, L3 <listOfLocalParameters>
        for list_tag in ("listOfLocalParameters", "listOfParameters"):
            lp = _find_local(kl, list_tag)
            if lp is not None:
                for p in _findall_local(lp, "parameter") + _findall_local(lp, "localParameter"):
                    pid = p.attrib.get("id", "")
                    val = p.attrib.get("value")
                    if val is not None:
                        try:
                            local_params[pid] = float(val)
                        except ValueError:
                            pass
        kinetic_math = _find_local(kl, "math")

    return {
        "id": rid,
        "reactants": reactants,
        "products": products,
        "modifiers": modifiers,
        "local_params": local_params,
        "math": kinetic_math,
        "reversible": r.attrib.get("reversible", "true").lower() == "true",
    }


# TOPOLOGY / CLASSIFICATION

def classify_compartments(model, diffusible_override=None, local_override=None):
    """
    Decide which compartments are "diffusible" (exchanged with a
    GridDiffusion field) vs "local" (tracked per cell only), by name
    keyword unless explicitly overridden.
    """
    diffusible_override = set(diffusible_override or [])
    local_override = set(local_override or [])
    diffusible = set()
    for cid, c in model.compartments.items():
        if cid in local_override:
            continue
        if cid in diffusible_override:
            diffusible.add(cid)
            continue
        name = (c["name"] or cid).lower()
        if any(kw in name or kw in cid.lower() for kw in _EXTRACELLULAR_KEYWORDS):
            diffusible.add(cid)
    return diffusible


def classify_species(model, diffusible_compartments, ignore_ids=None):
    """
    Split species into:
      - constant_species: boundaryCondition/constant True -> emitted as
        overridable numeric constants (like the SBOL converter's *_CONC).
      - local_species:    dynamic, tracked per cell (cell.<var> / species[]).
      - diffusible_species: dynamic, tracked per cell AND exchanged with a
        GridDiffusion grid via a configurable membrane rate.
    Only compartments actually present in the model can make a species
    diffusible; a model with no extracellular-like compartment naturally
    produces no diffusible species, mirroring the SBOL converter's
    behaviour of only turning signalling on when something is genuinely
    diffusible.
    """
    ignore_ids = set(ignore_ids or [])
    constant_species, local_species, diffusible_species = [], [], []

    for sid, s in model.species.items():
        if sid in ignore_ids:
            continue
        if s["boundary"] or s["constant"]:
            constant_species.append(sid)
        elif s["compartment"] in diffusible_compartments:
            diffusible_species.append(sid)
        else:
            local_species.append(sid)

    return constant_species, local_species, diffusible_species


def _is_fluorescent(sid, name):
    if sid in _FP_BIOBRICK:
        return _FP_BIOBRICK[sid]
    lowered = (name or sid).lower()
    for kw, rgb in _FP_KEYWORDS.items():
        if kw in lowered:
            return rgb
    return None


# KINETIC LAW -> PER-SPECIES ODE ASSEMBLY

def build_species_odes(model, species_var, species_index, tracked_ids):
    """
    For every tracked (non-constant) species, sum stoichiometry * reaction
    rate across every reaction that references it, in both Python and C
    syntax. Returns:
        odes_py: id -> python expression string (may be "" if unused)
        odes_c:  id -> C expression string
        warnings: list[str]
    """
    warnings = []

    all_symbol_values = dict(model.parameters)
    for cid, c in model.compartments.items():
        all_symbol_values.setdefault(cid, c["size"])

    contributions_py = {sid: [] for sid in tracked_ids}
    contributions_c  = {sid: [] for sid in tracked_ids}

    for rxn in model.reactions:
        if rxn["math"] is None:
            continue

        symbol_values = dict(all_symbol_values)
        symbol_values.update(rxn["local_params"])
        context = f"reaction '{rxn['id']}'"

        translator = MathTranslator(species_var, species_index, symbol_values, warnings, context)
        rate_py = translator.translate(rxn["math"], "py")
        rate_c  = translator.translate(rxn["math"], "c")

        stoich_by_species = {}
        for ref in rxn["reactants"]:
            stoich_by_species[ref["species"]] = stoich_by_species.get(ref["species"], 0.0) - ref["stoichiometry"]
        for ref in rxn["products"]:
            stoich_by_species[ref["species"]] = stoich_by_species.get(ref["species"], 0.0) + ref["stoichiometry"]

        for sid, coeff in stoich_by_species.items():
            if sid not in tracked_ids or coeff == 0.0:
                continue
            if coeff == 1.0:
                contributions_py[sid].append(f"({rate_py})")
                contributions_c[sid].append(f"({rate_c})")
            elif coeff == -1.0:
                contributions_py[sid].append(f"(-({rate_py}))")
                contributions_c[sid].append(f"(-({rate_c}))")
            else:
                coeff_py = _fmt_num(coeff, "py")
                coeff_c = _fmt_num(coeff, "c")
                contributions_py[sid].append(f"({coeff_py} * ({rate_py}))")
                contributions_c[sid].append(f"({coeff_c} * ({rate_c}))")

    odes_py = {sid: " + ".join(terms) for sid, terms in contributions_py.items()}
    odes_c  = {sid: " + ".join(terms) for sid, terms in contributions_c.items()}
    return odes_py, odes_c, warnings


# SCRIPT ASSEMBLY

def generate_script(model, params):
    """
    Build a full CellModeller simulation script directly from a parsed
    SBMLModel, using `params` for everything SBML doesn't specify itself
    (simulation constants, strain/cell-type visuals, grid geometry,
    per-species diffusion/membrane-exchange rates, constant-species
    concentrations). Mirrors the structure of
    cellmodeller_converter.generate_script() so the two families of
    generated scripts look and behave the same way.
    """
    sim_p   = params.get("simulation", {})
    cell_types = params.get("cell_types", [{}])
    sig_p   = params.get("signaling", {})
    mapping = params.get("sbml_mapping", {})

    max_cells    = sim_p.get("max_cells", 10000)
    jitter_z     = sim_p.get("jitter_z", False)
    gamma        = sim_p.get("gamma", 100.0)
    pickle_steps = sim_p.get("pickle_steps", 10)
    random_seed  = sim_p.get("random_seed", None)

    walls = params.get("walls", [])
    wall_lines = "\n".join(
        f"    biophys.addPlane(({w['point'][0]}, {w['point'][1]}, {w['point'][2]}), "
        f"({w['normal'][0]}, {w['normal'][1]}, {w['normal'][2]}), {w.get('coeff', 1.0)})"
        for w in walls
    )

    ignore_ids = set(mapping.get("ignore_species_ids", []))
    diffusible_compartments = classify_compartments(
        model,
        diffusible_override=mapping.get("diffusible_compartments"),
        local_override=mapping.get("local_compartments"),
    )
    constant_ids, local_ids, diffusible_ids = classify_species(
        model, diffusible_compartments, ignore_ids=ignore_ids
    )

    sig_enabled  = sig_p.get("enabled", bool(diffusible_ids)) and bool(diffusible_ids)
    tracked_diffusible_ids = diffusible_ids if sig_enabled else []
    # If signalling ends up disabled, any would-be-diffusible species falls
    # back to being tracked as an ordinary per-cell local species instead
    # of being silently dropped.
    fallback_local_ids = [] if sig_enabled else list(diffusible_ids)
    all_local_ids = local_ids + fallback_local_ids

    species_var = {sid: safe_name(model.species[sid]["name"] or sid)
                   for sid in all_local_ids + tracked_diffusible_ids}
    # de-duplicate any name collisions after safe_name() normalisation
    seen_names = {}
    for sid in all_local_ids + tracked_diffusible_ids:
        base = species_var[sid]
        n = seen_names.get(base, 0)
        if n:
            species_var[sid] = f"{base}_{n}"
        seen_names[base] = n + 1

    species_index = {sid: i for i, sid in enumerate(all_local_ids)}
    next_idx = len(all_local_ids)
    for sid in tracked_diffusible_ids:
        species_index[sid] = next_idx
        next_idx += 1

    tracked_ids = all_local_ids + tracked_diffusible_ids
    odes_py, odes_c, warnings = build_species_odes(model, species_var, species_index, tracked_ids)
    warnings.extend(model.warnings)

    # -- constant species -> overridable numeric constants --
    species_overrides = params.get("species_overrides", {})
    const_lines = []
    for sid in sorted(constant_ids):
        s = model.species[sid]
        cname = safe_name(s["name"] or sid).upper() + "_CONC"
        value = species_overrides.get(sid, s["initial"])
        const_lines.append(f"{cname} = {value}  # {sid} ({s['name']}) -- fixed/external species")
    const_name_by_id = {
        sid: safe_name(model.species[sid]["name"] or sid).upper() + "_CONC"
        for sid in constant_ids
    }
    # constants need to exist inside the OpenCL kernel text too (it is
    # compiled as a standalone C source string, so Python globals aren't
    # visible to it) -- emit #define lines at the top of each kernel.
    cl_defines = "\n".join(
        f"    #define {const_name_by_id[sid]} {species_overrides.get(sid, model.species[sid]['initial'])}f"
        for sid in sorted(constant_ids)
    ) if constant_ids else ""

    # Re-resolve constant-species references inside kinetic laws: they were
    # translated as unresolved (-> 0) above because they aren't in
    # species_var/symbol_values. Do a second pass now that we know their
    # generated constant names, via a light regex-free approach: rebuild
    # odes with constants included in the symbol table up front instead.
    if constant_ids:
        const_symbol_values = {sid: species_overrides.get(sid, model.species[sid]["initial"])
                                for sid in constant_ids}
        odes_py, odes_c, warnings2 = _rebuild_odes_with_constants(
            model, species_var, species_index, tracked_ids,
            const_symbol_values, const_name_by_id
        )
        warnings = warnings2 + model.warnings

    if const_lines:
        chem_consts = (
            "\n# FIXED / EXTERNAL SPECIES CONCENTRATIONS  "
            "(boundaryCondition or constant species -- set before running)\n"
            + "\n".join(const_lines) + "\n"
        )
    else:
        chem_consts = ""

    # -- grid geometry (only meaningful if signalling is on) --
    grid_len     = max(3, sig_p.get("grid_len", 100))
    grid_z_cells = max(3, sig_p.get("grid_z_cells", 3))
    grid_size    = sig_p.get("grid_size", 4.0)
    grid_origin  = sig_p.get("grid_origin", [
        -grid_len * grid_size / 2.0,
        -grid_len * grid_size / 2.0,
        -grid_z_cells * grid_size / 2.0,
    ])
    boundary_cond_in = sig_p.get("boundary_condition", "reflect")
    boundcond = _BOUNDCOND_MAP.get(boundary_cond_in, "reflect")

    param_signal_map = {s["name"]: s for s in sig_p.get("signals", [])}
    all_signals = []
    for sid in sorted(tracked_diffusible_ids):
        sname = model.species[sid]["name"] or sid
        entry = param_signal_map.get(sid, param_signal_map.get(sname, {}))
        diff_rate = entry.get("diffusion_rate", 0.1)
        all_signals.append({
            "species_id":             sid,
            "name":                   sname,
            "diffusion_rate":         diff_rate,
            "membrane_exchange_rate": entry.get("membrane_exchange_rate", diff_rate),
            "initial_concentration":  entry.get("initial_concentration", 0.0),
        })
    signal_index   = {s["species_id"]: i for i, s in enumerate(all_signals)}
    membrane_rates = {s["species_id"]: s["membrane_exchange_rate"] for s in all_signals}
    n_signals      = len(all_signals)
    use_signals    = sig_enabled and n_signals > 0 and len(tracked_ids) > 0

    if not use_signals:
        # signalling ended up empty after all (e.g. no reactions reference
        # any diffusible-compartment species) -- fully disable it.
        tracked_diffusible_ids = []

    # -- cell type / strain tables (identical shape to the SBOL converter) --
    color_lines, len_lines, growth_lines, noise_lines, conc_lines = [], [], [], [], []
    for i, ct in enumerate(cell_types):
        c = ct.get("color", [1.0, 0.3, 0.3])
        label = ct.get("display_name", f"Strain {i}")
        color_lines.append(f"    {i}: [{c[0]}, {c[1]}, {c[2]}],  # {label}")
        len_lines.append(f"    {i}: {ct.get('division_length', 3.5)},")
        growth_lines.append(f"    {i}: {ct.get('growth_rate', 1.0)},")
        noise_lines.append(f"    {i}: {ct.get('division_noise', 0.005)},")
        conc = ct.get("initial_concentrations", {})
        entries = ", ".join(
            f"'{species_var[sid]}': {conc.get(sid, model.species[sid]['initial'])}"
            for sid in all_local_ids
        )
        conc_lines.append(f"    {i}: {{{entries}}},  # {label}")

    colors_dict = "{\n" + "\n".join(color_lines) + "\n}"
    lens_dict   = "{\n" + "\n".join(len_lines) + "\n}"
    growth_dict = "{\n" + "\n".join(growth_lines) + "\n}"
    noise_dict  = "{\n" + "\n".join(noise_lines) + "\n}"
    conc_dict   = ("{\n" + "\n".join(conc_lines) + "\n}") if all_local_ids else "{}"

    add_cell_lines = []
    for i, ct in enumerate(cell_types):
        pos = ct.get("initial_pos", [0.0, float(i) * 3.0, 0.0])
        d = ct.get("initial_dir", [1.0, 0.0, 0.0])
        add_cell_lines.append(
            f"    sim.addCell(cellType={i}, "
            f"pos=({pos[0]}, {pos[1]}, {pos[2]}), "
            f"dir=({d[0]}, {d[1]}, {d[2]}))"
        )
    add_cells_str = "\n".join(add_cell_lines)

    if use_signals:
        n_species = len(species_index)
        diff_rates = [str(s["diffusion_rate"]) for s in all_signals]
        sig_import1 = "from CellModeller.Signalling.GridDiffusion import GridDiffusion"
        sig_import2 = "from CellModeller.Integration.CLCrankNicIntegrator import CLCrankNicIntegrator"
        sig_init = (
            f"    nSignals = {n_signals}  # {[s['name'] for s in all_signals]}\n"
            f"    nSpecies = {n_species}  # tracked per-cell species (local + diffusible pools)\n"
            f"    sig   = GridDiffusion(sim, nSignals,\n"
            f"                (gridLen, gridLen, gridZCells),\n"
            f"                (gridSize, gridSize, gridSize),   # must be isotropic\n"
            f"                (gridOrigX, gridOrigY, gridOrigZ),\n"
            f"                [{', '.join(diff_rates)}])   # diffusion coefficients\n"
            f"    integ = CLCrankNicIntegrator(sim, nSignals, nSpecies, maxCells,\n"
            f"                sig, boundcond='{boundcond}')  # {boundary_cond_in} boundary"
        )
        sim_init_call = "    sim.init(biophys, regul, sig, integ)"
        grid_constants_block = (
            f"gridLen     = {grid_len}       # grid cells along x and y\n"
            f"gridZCells  = {grid_z_cells}   # grid cells along z (small = thin/2D)\n"
            f"gridSize    = {grid_size}      # micrometers per grid cell (must be isotropic)\n"
            f"gridOrigX, gridOrigY, gridOrigZ = {grid_origin[0]}, {grid_origin[1]}, {grid_origin[2]}"
        )
    else:
        sig_import1 = "# from CellModeller.Signalling.GridDiffusion import GridDiffusion"
        sig_import2 = "# from CellModeller.Integration.CLCrankNicIntegrator import CLCrankNicIntegrator"
        sig_init = "    # Signalling disabled — no diffusible-compartment species reference a reaction"
        sim_init_call = "    sim.init(biophys, regul, None, None)"
        grid_constants_block = (
            f"gridLen  = {grid_len}   # grid cells per axis\n"
            f"gridSize = {grid_size}  # µm per grid cell"
        )

    random_seed_line = (
        f"    random.seed({random_seed})"
        if random_seed is not None
        else "    # tip: set simulation.random_seed in params for reproducibility"
    )

    # -- topology summary comment --
    topo_lines = ["# AUTO-DETECTED SBML TOPOLOGY"]
    topo_lines.append(f"#   Compartments             : {list(model.compartments.keys())}")
    topo_lines.append(f"#   Diffusible compartments  : {sorted(diffusible_compartments)}")
    topo_lines.append(f"#   Local species            : {all_local_ids}")
    topo_lines.append(f"#   Diffusible species       : {tracked_diffusible_ids}")
    topo_lines.append(f"#   Fixed/external species   : {sorted(constant_ids)}")
    if warnings:
        topo_lines.append("#   Conversion warnings:")
        for w in warnings:
            topo_lines.append(f"#     - {w}")
    topo_comment = "\n".join(topo_lines) + "\n"

    # -- fluorescent-protein colour-by-expression --
    proteins_for_color = [
        {"display_id": sid, "var_name": species_var[sid]}
        for sid in tracked_ids
        if _is_fluorescent(sid, model.species[sid]["name"])
    ]
    color_update_raw = _generate_color_update(
        proteins_for_color, species_index=species_index if use_signals else None
    )
    color_update_str = ("\n" + color_update_raw + "\n") if color_update_raw else ""

    if use_signals:
        init_lines = [f"        cell_initial_concentrations[cell.cellType]['{species_var[sid]}'],"
                      f"  # {sid} -> species[{species_index[sid]}]" for sid in all_local_ids]
        init_lines += [
            f"        {next(s for s in all_signals if s['species_id'] == sid)['initial_concentration']},"
            f"  # {sid} -> species[{species_index[sid]}]"
            for sid in tracked_diffusible_ids
        ]
        init_species = (
            "    # local + diffusible species, all tracked via cell.species[]\n"
            "    cell.species[:] = [\n" + "\n".join(init_lines) + "\n    ]"
        )
        divide_species = (
            "    d1.species[:] = parent.species[:] / 2.0\n"
            "    d2.species[:] = parent.species[:] / 2.0"
        )
        init_signals = (
            f"\n    # locally-sensed extracellular signal levels\n"
            f"    cell.signals[:] = [0.0] * {n_signals}"
        )
        divide_signals = ""

        specratecl_body = _generate_specratecl(
            tracked_ids, all_local_ids, tracked_diffusible_ids, species_index,
            signal_index, odes_c, membrane_rates, cl_defines,
        )
        sigratecl_body = _generate_sigratecl(tracked_diffusible_ids, species_index, signal_index, membrane_rates)

        cl_functions = f'''

# SPEC / SIGNAL RATES (OpenCL C — required by CLCrankNicIntegrator)

def specRateCL():
    return \'\'\'
{specratecl_body}
    \'\'\'

def sigRateCL():
    return \'\'\'
{sigratecl_body}
    \'\'\'
'''
        update_logic = ""
    else:
        if all_local_ids:
            init_species = "    # species\n" + "\n".join(
                f"    cell.{species_var[sid]} = "
                f"cell_initial_concentrations[cell.cellType]['{species_var[sid]}']"
                f"  # {sid}"
                for sid in all_local_ids
            )
            divide_species = "\n".join(
                f"    d1.{species_var[sid]} = parent.{species_var[sid]} / 2.0\n"
                f"    d2.{species_var[sid]} = parent.{species_var[sid]} / 2.0"
                for sid in all_local_ids
            )
        else:
            init_species = "    # no dynamic species in this SBML model"
            divide_species = ""
        init_signals = ""
        divide_signals = ""
        cl_functions = ""
        update_logic = _generate_update_logic(all_local_ids, species_var, odes_py)

    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    update_block = f"\n{update_logic}\n" if update_logic else ""

    script = f"""\
# CellModeller simulation script
# Generated : {now}
# Converter : sbml_to_cellmodeller.py  (direct SBML -> CellModeller, no SBOL/JSON bridge)
{"# NOTE: signalling is ON — species kinetics live in specRateCL(), not update()." if use_signals else ""}


from CellModeller.Regulation.ModuleRegulator import ModuleRegulator
from CellModeller.Biophysics.BacterialModels.CLBacterium import CLBacterium
{sig_import1}
{sig_import2}
import numpy as np
import random
import math
{chem_consts}
{topo_comment}
# simulation constants
maxCells = {max_cells}
{grid_constants_block}


# cell type lookup tables

cell_colors          = {colors_dict}
cell_lens            = {lens_dict}
cell_growth_rates    = {growth_dict}
cell_division_noise  = {noise_dict}
cell_initial_concentrations = {conc_dict}


# SETUP

def setup(sim):
{random_seed_line}
    biophys = CLBacterium(
        sim,
        jitter_z={jitter_z},
        max_cells=maxCells,
        max_planes={max(len(walls), 1)},
        gamma={gamma}
    )
{wall_lines}
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
{init_species}{init_signals}


# UPDATE

def update(cells):
    for id, cell in cells.items():
        if cell.volume > cell.targetVol:
            cell.divideFlag = True
{color_update_str}{update_block}

# DIVIDE

def divide(parent, d1, d2):
    d1.cellType = parent.cellType
    d2.cellType = parent.cellType
    d1.color    = parent.color
    d2.color    = parent.color
{divide_species}{divide_signals}
{cl_functions}"""
    return script, warnings


def _rebuild_odes_with_constants(model, species_var, species_index, tracked_ids,
                                  const_values, const_name_by_id):
    """Same as build_species_odes, but constant/boundary species resolve to
    their named _CONC constant (py: bare name; c: #define'd macro) instead
    of an 'unresolved identifier' warning."""
    warnings = []
    all_symbol_values = dict(model.parameters)
    for cid, c in model.compartments.items():
        all_symbol_values.setdefault(cid, c["size"])

    contributions_py = {sid: [] for sid in tracked_ids}
    contributions_c  = {sid: [] for sid in tracked_ids}

    for rxn in model.reactions:
        if rxn["math"] is None:
            continue
        symbol_values = dict(all_symbol_values)
        symbol_values.update(rxn["local_params"])
        context = f"reaction '{rxn['id']}'"

        translator = MathTranslator(species_var, species_index, symbol_values, warnings, context)

        # monkey-patch constant-species resolution for this translator
        orig_resolve = translator._resolve_ci

        def resolve_with_const(name, lang, _orig=orig_resolve):
            if name in const_name_by_id:
                return const_name_by_id[name] if lang == "c" else const_name_by_id[name]
            return _orig(name, lang)

        translator._resolve_ci = resolve_with_const

        rate_py = translator.translate(rxn["math"], "py")
        rate_c  = translator.translate(rxn["math"], "c")

        stoich_by_species = {}
        for ref in rxn["reactants"]:
            stoich_by_species[ref["species"]] = stoich_by_species.get(ref["species"], 0.0) - ref["stoichiometry"]
        for ref in rxn["products"]:
            stoich_by_species[ref["species"]] = stoich_by_species.get(ref["species"], 0.0) + ref["stoichiometry"]

        for sid, coeff in stoich_by_species.items():
            if sid not in tracked_ids or coeff == 0.0:
                continue
            if coeff == 1.0:
                contributions_py[sid].append(f"({rate_py})")
                contributions_c[sid].append(f"({rate_c})")
            elif coeff == -1.0:
                contributions_py[sid].append(f"(-({rate_py}))")
                contributions_c[sid].append(f"(-({rate_c}))")
            else:
                contributions_py[sid].append(f"({_fmt_num(coeff, 'py')} * ({rate_py}))")
                contributions_c[sid].append(f"({_fmt_num(coeff, 'c')} * ({rate_c}))")

    odes_py = {sid: " + ".join(terms) for sid, terms in contributions_py.items()}
    odes_c  = {sid: " + ".join(terms) for sid, terms in contributions_c.items()}
    return odes_py, odes_c, warnings


def _generate_update_logic(local_ids, species_var, odes_py):
    if not local_ids:
        return "        pass"
    lines = ["        # — species production/consumption, summed from every SBML reaction —"]
    for sid in local_ids:
        var = species_var[sid]
        expr = odes_py.get(sid, "")
        if expr:
            lines.append(f"        cell.{var} += {expr}  # {sid}")
        else:
            lines.append(f"        # {sid}: not referenced by any reaction (held at initial value)")
        lines.append(f"        cell.{var} = max(0.0, cell.{var})")
        lines.append("")
    return "\n".join(lines)


def _generate_color_update(proteins, species_index=None):
    """Same shading logic as cellmodeller_converter.generate_color_update,
    kept local so this module has no dependency on a `params`-style
    kinetics dict (SBML species don't share one global repression
    threshold — a fixed normalisation constant is used instead)."""
    NORM = 2.0
    for protein in proteins:
        pid, var = protein["display_id"], protein["var_name"]
        rgb = _is_fluorescent(pid, pid)
        if rgb is None:
            continue
        r, g, b = rgb
        expr = f"cell.species[{species_index[pid]}]" if species_index else f"cell.{var}"
        return (
            f"        # colour by {pid} expression\n"
            f"        _fp_norm = min(1.0, {expr} / {NORM})\n"
            f"        cell.color = [\n"
            f"            {r} * _fp_norm + 0.1 * (1.0 - _fp_norm),\n"
            f"            {g} * _fp_norm + 0.1 * (1.0 - _fp_norm),\n"
            f"            {b} * _fp_norm + 0.1 * (1.0 - _fp_norm),\n"
            f"        ]"
        )
    return None


def _generate_specratecl(tracked_ids, local_ids, diffusible_ids, species_index,
                         signal_index, odes_c, membrane_rates, cl_defines):
    lines = []
    if cl_defines:
        lines.append("    // — fixed/external species constants —")
        lines.append(cl_defines)
        lines.append("")

    lines.append("    // — species production/consumption, summed from every SBML reaction —")
    for sid in local_ids:
        idx = species_index[sid]
        expr = odes_c.get(sid, "")
        if expr:
            lines.append(f"    rates[{idx}] = {expr};  // {sid}")
        else:
            lines.append(f"    rates[{idx}] = 0.0f;  // {sid}: not referenced by any reaction")
    lines.append("")

    if diffusible_ids:
        lines.append("    // — diffusible-compartment species pools "
                      "(reaction kinetics here, grid exchange in sigRateCL) —")
        for sid in diffusible_ids:
            idx = species_index[sid]
            gidx = signal_index[sid]
            mrate = membrane_rates.get(sid, 0.1)
            expr = odes_c.get(sid, "0.0f")
            lines.append(
                f"    rates[{idx}] = ({expr if expr else '0.0f'})"
                f" - {mrate}f * (species[{idx}] - signals[{gidx}]) * area / volume;  // {sid}"
            )
        lines.append("")

    return "\n".join(lines) if lines else "    // nothing to do"


def _generate_sigratecl(diffusible_ids, species_index, signal_index, membrane_rates):
    if not diffusible_ids:
        return "    // no diffusible species to exchange with the grid"
    lines = []
    for sid in diffusible_ids:
        idx = species_index[sid]
        gidx = signal_index[sid]
        mrate = membrane_rates.get(sid, 0.1)
        lines.append(f"    // {sid}: secretion into the grid (species[{idx}] -> signals[{gidx}])")
        lines.append(
            f"    rates[{gidx}] = {mrate}f * (species[{idx}] - signals[{gidx}]) * area / gridVolume;"
        )
    return "\n".join(lines)


# CLI ENTRY POINT

def main():
    parser = argparse.ArgumentParser(
        description="Convert an SBML file directly to a CellModeller simulation script."
    )
    parser.add_argument("--sbml", required=True, help="Path to SBML XML")
    parser.add_argument("--params", required=True, help="Path to parameters JSON")
    parser.add_argument("--output", default="cellmodeller_output.py",
                        help="Output .py file (default: cellmodeller_output.py)")
    args = parser.parse_args()

    print(f"Loading SBML   : {args.sbml}")
    model = parse_sbml(args.sbml)
    print(f"Loading params : {args.params}")
    with open(args.params) as f:
        params = json.load(f)

    print(f"  Compartments : {list(model.compartments.keys())}")
    print(f"  Species      : {list(model.species.keys())}")
    print(f"  Reactions    : {[r['id'] for r in model.reactions]}")

    script, warnings = generate_script(model, params)
    if warnings:
        print("Conversion warnings:")
        for w in warnings:
            print(f"  - {w}")

    with open(args.output, "w") as f:
        f.write(script)
    print(f"Written: {args.output}")


if __name__ == "__main__":
    main()
