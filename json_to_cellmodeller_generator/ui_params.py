"""
ui_params.py

Interactive ipywidgets form for building the `parameters` dict consumed by
`generate_script()` in cellmodeller_converter.py, with first-class support
for signalling (diffusible chemicals such as AHL in quorum sensing).

Usage (in a Jupyter notebook):

    import json
    from ui_params import display_form

    sbol_data = json.load(open("circuit.json"))
    parameters, topology = display_form(sbol_data)

    # ... adjust widgets, click "Generate CellModeller script" ...
    # `parameters` is updated in place every time "Save parameters" or
    # "Generate CellModeller script" is clicked, so the same dict you got
    # back above will reflect the latest edits — e.g. you can later do:

    from cellmodeller_converter import parse_json, generate_script
    proteins, modules, interactions, component_map, ed_chemicals = parse_json(sbol_data)
    script = generate_script(proteins, modules, interactions, component_map,
                              parameters, ed_chemicals)
"""

import ipywidgets as widgets
from IPython.display import display, clear_output

from cellmodeller_converter import parse_json, generate_script
from params import (
    DEFAULT_SIMULATION,
    DEFAULT_CELL_TYPE,
    DEFAULT_KINETICS,
    DEFAULT_SIGNAL_RATES,
    DEFAULT_SIGNALING,
    detect_signals,
    detect_reporter_color,
    COLOR_NAME_TO_RGB,
)


# small layout helpers

def _spacer(width="15px"):
    return widgets.Label(value="", layout=widgets.Layout(width=width))


def _section(title, *children):
    return widgets.VBox(
        [widgets.HTML(value=f"<h3>{title}</h3>")] + list(children),
        layout=widgets.Layout(
            border="1px solid lightgray", margin="10px 0", padding="10px"
        ),
    )


def _rgb_to_hex(rgb):
    r, g, b = (max(0, min(255, int(round(c * 255)))) for c in rgb)
    return f"#{r:02x}{g:02x}{b:02x}"


def _hex_to_rgb(hex_color):
    hex_color = (hex_color or "#ff4d4d").lstrip("#")
    r = int(hex_color[0:2], 16) / 255.0
    g = int(hex_color[2:4], 16) / 255.0
    b = int(hex_color[4:6], 16) / 255.0
    return [round(r, 3), round(g, 3), round(b, 3)]


# main entry point

def display_form(sbol_data, ignore_component_ids=None):
    """
    Build and display an ipywidgets form for the CellModeller simulation
    parameters, pre-populated with the signalling topology auto-detected
    from `sbol_data` (diffusible signals, their producer proteins, and the
    promoters they activate).

    Args:
      sbol_data (dict): SBOL JSON circuit description (keys "ED",
                        "interactions", "hierarchy", "components").
      ignore_component_ids (list[str], optional): DNA-part displayIds to
                        exclude up front (passed straight through to
                        `parse_json`).

    Returns:
      tuple:
        - parameters (dict): the live parameters dict. It starts populated
          with defaults / auto-detected signals, and is updated in place
          every time "Save parameters" or "Generate CellModeller script" is
          clicked — so the reference returned here always reflects the
          latest form state once you've clicked one of those buttons.
        - topology (dict): the auto-detected signalling topology (see
          `params.detect_signals`), for printing a summary or building
          further tooling.
    """
    topology = detect_signals(sbol_data, ignore_component_ids=ignore_component_ids)
    proteins = topology["proteins"]
    diffusible_signals = topology["diffusible_signals"]
    signal_producers = topology["signal_producers"]
    signal_activated_promoters = topology["signal_activated_promoters"]

    parameters = {
        "simulation": dict(DEFAULT_SIMULATION),
        "cell_types": [],
        "signaling": {**DEFAULT_SIGNALING, "signals": []},
        "kinetics": dict(DEFAULT_KINETICS),
        "sbol_mapping": {"ignore_component_ids": list(ignore_component_ids or [])},
    }

    output_area = widgets.Output()

    # Simulation section
    max_cells_w = widgets.IntText(
        value=DEFAULT_SIMULATION["max_cells"], description="Max cells:",
        style={"description_width": "initial"},
    )
    gamma_w = widgets.FloatText(
        value=DEFAULT_SIMULATION["gamma"], description="Gamma (friction):",
        style={"description_width": "initial"},
    )
    pickle_steps_w = widgets.IntText(
        value=DEFAULT_SIMULATION["pickle_steps"], description="Pickle every N steps:",
        style={"description_width": "initial"},
    )
    jitter_z_w = widgets.Checkbox(
        value=DEFAULT_SIMULATION["jitter_z"], description="Jitter Z (3-D wobble)",
        indent=False,
    )
    use_seed_w = widgets.Checkbox(value=False, description="Fix random seed", indent=False)
    seed_w = widgets.IntText(
        value=0, description="Seed:", style={"description_width": "initial"}, disabled=True,
    )

    def _toggle_seed(change):
        seed_w.disabled = not change["new"]

    use_seed_w.observe(_toggle_seed, names="value")

    simulation_box = _section(
        "Simulation",
        widgets.HBox([max_cells_w, _spacer(), gamma_w, _spacer(), pickle_steps_w]),
        widgets.HBox([jitter_z_w, _spacer(), use_seed_w, _spacer(), seed_w]),
    )

    # Cell types / strains section
    cell_type_entries = []
    cell_types_container = widgets.VBox([])

    def _refresh_cell_types():
        cell_types_container.children = tuple(e["box"] for e in cell_type_entries)
        for i, e in enumerate(cell_type_entries):
            e["index_label"].value = f"<b>Strain {i}</b>"

    def _make_cell_type_box(preset=None):
        preset = preset or {}
        index_label = widgets.HTML(value="<b>Strain</b>")
        name_w = widgets.Text(
            value=preset.get("display_name", "Strain"), description="Name:",
            style={"description_width": "initial"},
        )
        color_w = widgets.ColorPicker(
            value=_rgb_to_hex(preset.get("color", DEFAULT_CELL_TYPE["color"])),
            description="Colour:", style={"description_width": "initial"},
        )
        div_len_w = widgets.FloatText(
            value=preset.get("division_length", DEFAULT_CELL_TYPE["division_length"]),
            description="Division length (µm):", style={"description_width": "initial"},
        )
        growth_w = widgets.FloatText(
            value=preset.get("growth_rate", DEFAULT_CELL_TYPE["growth_rate"]),
            description="Growth rate:", style={"description_width": "initial"},
        )
        noise_w = widgets.FloatText(
            value=preset.get("division_noise", DEFAULT_CELL_TYPE["division_noise"]),
            description="Division noise:", style={"description_width": "initial"},
        )
        init_pos = preset.get("initial_pos", DEFAULT_CELL_TYPE["initial_pos"])
        pos_x_w = widgets.FloatText(value=init_pos[0], description="Start X:",
                                     style={"description_width": "initial"},
                                     layout=widgets.Layout(width="120px"))
        pos_y_w = widgets.FloatText(value=init_pos[1], description="Y:",
                                     style={"description_width": "initial"},
                                     layout=widgets.Layout(width="100px"))
        pos_z_w = widgets.FloatText(value=init_pos[2], description="Z:",
                                     style={"description_width": "initial"},
                                     layout=widgets.Layout(width="100px"))
        remove_btn = widgets.Button(
            description="Remove strain", button_style="danger",
            layout=widgets.Layout(width="130px"),
        )
        
        conc_widgets = {}
        conc_row_items = [widgets.Label(value="Initial protein amounts:")]
        preset_conc = preset.get("initial_concentrations", {})
        for p in proteins:
            w = widgets.FloatText(
                value=preset_conc.get(p["display_id"], 0.0),
                description=f"{p['display_id']}:",
                style={"description_width": "initial"},
                layout=widgets.Layout(width="150px"),
            )
            conc_widgets[p["display_id"]] = w
            conc_row_items.append(w)
        conc_row = widgets.HBox(conc_row_items)
        
        row1 = widgets.HBox([index_label, _spacer(), name_w, _spacer(), color_w, _spacer(), remove_btn])
        row2 = widgets.HBox([div_len_w, _spacer(), growth_w, _spacer(), noise_w])
        row3 = widgets.HBox([widgets.Label(value="Initial position:"), pos_x_w, pos_y_w, pos_z_w])
        rows = [row1, row2, row3]
        if proteins:
            rows.append(conc_row)
        box = widgets.VBox(
            rows,
            layout=widgets.Layout(border="1px dashed lightgray", margin="5px 0", padding="8px"),
        )

        entry = {
            "box": box,
            "index_label": index_label,
            "widgets": {
                "name": name_w, "color": color_w, "div_len": div_len_w,
                "growth": growth_w, "noise": noise_w,
                "pos_x": pos_x_w, "pos_y": pos_y_w, "pos_z": pos_z_w,
                "conc": conc_widgets,
            },
        }

        def _on_remove(_btn):
            if len(cell_type_entries) <= 1:
                with output_area:
                    clear_output(wait=True)
                    print("At least one strain is required.")
                return
            if entry in cell_type_entries:
                cell_type_entries.remove(entry)
                _refresh_cell_types()

        remove_btn.on_click(_on_remove)
        return entry

    def _add_cell_type(_btn=None, preset=None):
        cell_type_entries.append(_make_cell_type_box(preset=preset))
        _refresh_cell_types()

    add_strain_btn = widgets.Button(description="Add strain", layout=widgets.Layout(width="120px"))
    add_strain_btn.on_click(_add_cell_type)

    # Seed the first strain, suggesting a colour from a detected reporter
    # protein (e.g. GFP to green) so the default already looks right.
    first_preset = dict(DEFAULT_CELL_TYPE)
    first_preset["display_name"] = "Strain 0"
    for protein in proteins:
        color_name = detect_reporter_color(protein["display_id"])
        if color_name and color_name in COLOR_NAME_TO_RGB:
            first_preset["color"] = COLOR_NAME_TO_RGB[color_name]
            first_preset["display_name"] = f"Strain 0 ({protein['display_id']})"
            break
    _add_cell_type(preset=first_preset)

    cell_types_box = _section("Cell types / strains", add_strain_btn, cell_types_container)

    # Kinetics section
    prod_rate_w = widgets.FloatText(
        value=DEFAULT_KINETICS["production_rate"],
        description="Constitutive production rate:", style={"description_width": "initial"},
    )
    max_rate_w = widgets.FloatText(
        value=DEFAULT_KINETICS["max_production_rate"],
        description="Max regulated production rate:", style={"description_width": "initial"},
    )
    deg_rate_w = widgets.FloatText(
        value=DEFAULT_KINETICS["degradation_rate"],
        description="Protein degradation rate:", style={"description_width": "initial"},
    )
    hill_w = widgets.FloatText(
        value=DEFAULT_KINETICS["hill_coefficient"],
        description="Hill coefficient:", style={"description_width": "initial"},
    )
    act_thr_w = widgets.FloatText(
        value=DEFAULT_KINETICS["activation_threshold"],
        description="Activation threshold:", style={"description_width": "initial"},
    )
    rep_thr_w = widgets.FloatText(
        value=DEFAULT_KINETICS["repression_threshold"],
        description="Repression threshold:", style={"description_width": "initial"},
    )
    sig_prod_w = widgets.FloatText(
        value=DEFAULT_KINETICS["signal_production_rate"],
        description="Signal emission rate (per unit producer protein):",
        style={"description_width": "initial"},
    )

    kinetics_box = _section(
        "Reaction kinetics  (shared Hill-function parameters for every gene)",
        widgets.HBox([prod_rate_w, _spacer(), max_rate_w, _spacer(), deg_rate_w]),
        widgets.HBox([hill_w, _spacer(), act_thr_w, _spacer(), rep_thr_w]),
        widgets.HBox([sig_prod_w]),
    )

    # Signalling section
    sig_enabled_w = widgets.Checkbox(
        value=bool(diffusible_signals),
        description="Enable inter-cell signalling (GridDiffusion)", indent=False,
    )
    grid_len_w = widgets.IntText(
        value=DEFAULT_SIGNALING["grid_len"], description="Grid cells (x, y):",
        style={"description_width": "initial"},
        layout=widgets.Layout(width="220px"),
    )
    grid_z_cells_w = widgets.BoundedIntText(
        value=DEFAULT_SIGNALING["grid_z_cells"], min=2, max=1000,
        description="Grid cells (z):",
        style={"description_width": "initial"},
        layout=widgets.Layout(width="220px"),
    )
    grid_size_w = widgets.FloatText(
        value=DEFAULT_SIGNALING["grid_size"], description="µm per grid cell:",
        style={"description_width": "initial"},
        layout=widgets.Layout(width="220px"),
    )
    boundary_w = widgets.Dropdown(
        options=["reflect", "fixed", "periodic"], value=DEFAULT_SIGNALING["boundary_condition"],
        description="Boundary:", style={"description_width": "initial"},
        layout=widgets.Layout(width="220px"),
    )
    grid_help_w = widgets.HTML(
        value=(
            "<i>Keep z small (minimum 2, e.g. 2-4) for a 2D/monolayer simulation — "
            "µm per grid cell must be the same in x, y and z.</i>"
        )
    )

    signal_rate_widgets = {}  # signal display_id: {"diffusion": w, "membrane": w}

    if diffusible_signals:
        signal_rows = []
        for sid in sorted(diffusible_signals):
            producers = signal_producers.get(sid, [])
            activates = sorted(
                prom for prom, sig in signal_activated_promoters.items() if sig == sid
            )
            info_text = f"<b>{sid}</b>"
            info_text += (
                f" — produced by: {', '.join(producers)}"
                if producers else " — no producer auto-detected; set rates manually if needed"
            )
            if activates:
                info_text += f"; activates: {', '.join(activates)}"
            info_html = widgets.HTML(value=info_text)

            diff_w = widgets.FloatText(
                value=DEFAULT_SIGNAL_RATES["diffusion_rate"], description="Diffusion rate:",
                style={"description_width": "initial"},
            )
            mem_w = widgets.FloatText(
                value=DEFAULT_SIGNAL_RATES["membrane_exchange_rate"],
                description="Membrane exchange rate:",
                style={"description_width": "initial"},
            )
            row = widgets.VBox(
                [info_html, widgets.HBox([diff_w, _spacer(), mem_w])],
                layout=widgets.Layout(border="1px dashed lightgray", margin="5px 0", padding="8px"),
            )
            signal_rows.append(row)
            signal_rate_widgets[sid] = {"diffusion": diff_w, "membrane": mem_w}
        detected_signals_box = widgets.VBox(signal_rows)
    else:
        detected_signals_box = widgets.HTML(
            value=(
                "<i>No diffusible signal was auto-detected in this circuit "
                "(no chemical acts as a Stimulator in a Stimulation interaction). "
                "Signalling will be skipped unless you add a custom signal below "
                "or your SBOL is updated to include one.</i>"
            )
        )

    # Manual/custom signal support, for signals the SBOL topology doesn't
    # capture but you still want GridDiffusion to simulate.
    custom_signal_entries = []
    custom_signals_container = widgets.VBox([])

    def _refresh_custom_signals():
        custom_signals_container.children = tuple(e["box"] for e in custom_signal_entries)

    def _add_custom_signal(_btn=None):
        name_w = widgets.Text(value="", description="Signal name:", style={"description_width": "initial"})
        diff_w = widgets.FloatText(
            value=DEFAULT_SIGNAL_RATES["diffusion_rate"], description="Diffusion rate:",
            style={"description_width": "initial"},
        )
        mem_w = widgets.FloatText(
            value=DEFAULT_SIGNAL_RATES["membrane_exchange_rate"], description="Membrane exchange rate:",
            style={"description_width": "initial"},
        )
        remove_btn = widgets.Button(description="Remove", button_style="danger",
                                     layout=widgets.Layout(width="90px"))
        row = widgets.HBox([name_w, _spacer(), diff_w, _spacer(), mem_w, _spacer(), remove_btn])
        box = widgets.VBox([row], layout=widgets.Layout(border="1px dashed lightgray", margin="5px 0", padding="8px"))
        entry = {"box": box, "widgets": {"name": name_w, "diffusion": diff_w, "membrane": mem_w}}

        def _on_remove(_b):
            if entry in custom_signal_entries:
                custom_signal_entries.remove(entry)
                _refresh_custom_signals()

        remove_btn.on_click(_on_remove)
        custom_signal_entries.append(entry)
        _refresh_custom_signals()

    add_custom_signal_btn = widgets.Button(description="Add custom signal", layout=widgets.Layout(width="160px"))
    add_custom_signal_btn.on_click(_add_custom_signal)

    signaling_box = _section(
        "Signalling",
        widgets.HBox([sig_enabled_w]),
        widgets.HBox([grid_len_w, _spacer(), grid_z_cells_w, _spacer(), grid_size_w, _spacer(), boundary_w]),
        grid_help_w,
        widgets.HTML(value="<b>Auto-detected diffusible signals</b>"),
        detected_signals_box,
        widgets.HTML(value="<b>Custom signals</b> (not detected from SBOL — optional)"),
        add_custom_signal_btn,
        custom_signals_container,
    )

    # sbol_mapping
    ignore_ids_w = widgets.Text(
        value=", ".join(ignore_component_ids or []),
        description="Ignore component IDs (comma-separated):",
        style={"description_width": "initial"},
        layout=widgets.Layout(width="600px"),
    )
    advanced_box = _section(
        "Advanced",
        ignore_ids_w,
        widgets.HTML(
            value=(
                "<i>Changing this re-filters which DNA parts are considered "
                "when re-parsing the circuit on Save/Generate. If you change "
                "it and it affects which signals are detected, re-run "
                "display_form() to refresh the signalling widgets above.</i>"
            )
        ),
    )

    # Actions: Save / Generate
    script_preview_w = widgets.Textarea(
        value="", layout=widgets.Layout(width="100%", height="350px"), disabled=True,
    )
    output_path_w = widgets.Text(
        value="cellmodeller_output.py", description="Save script to:",
        style={"description_width": "initial"},
    )

    def _current_ignore_ids():
        return [s.strip() for s in ignore_ids_w.value.split(",") if s.strip()]

    def _collect_parameters():
        parameters["simulation"] = {
            "max_cells": max_cells_w.value,
            "jitter_z": jitter_z_w.value,
            "gamma": gamma_w.value,
            "pickle_steps": pickle_steps_w.value,
            "random_seed": seed_w.value if use_seed_w.value else None,
        }

        cell_types = []
        for entry in cell_type_entries:
            w = entry["widgets"]
            cell_types.append({
                "display_name": w["name"].value,
                "color": _hex_to_rgb(w["color"].value),
                "division_length": w["div_len"].value,
                "growth_rate": w["growth"].value,
                "division_noise": w["noise"].value,
                "initial_pos": [w["pos_x"].value, w["pos_y"].value, w["pos_z"].value],
                "initial_dir": [1.0, 0.0, 0.0],
                "initial_concentrations": {pid: cw.value for pid, cw in w["conc"].items()},
            })
        parameters["cell_types"] = cell_types

        parameters["kinetics"] = {
            "production_rate": prod_rate_w.value,
            "max_production_rate": max_rate_w.value,
            "degradation_rate": deg_rate_w.value,
            "hill_coefficient": hill_w.value,
            "activation_threshold": act_thr_w.value,
            "repression_threshold": rep_thr_w.value,
            "signal_production_rate": sig_prod_w.value,
        }

        signals = []
        for sid, rw in signal_rate_widgets.items():
            signals.append({
                "name": sid,
                "diffusion_rate": rw["diffusion"].value,
                "membrane_exchange_rate": rw["membrane"].value,
            })
        for entry in custom_signal_entries:
            w = entry["widgets"]
            if w["name"].value.strip():
                signals.append({
                    "name": w["name"].value.strip(),
                    "diffusion_rate": w["diffusion"].value,
                    "membrane_exchange_rate": w["membrane"].value,
                })

        parameters["signaling"] = {
            "enabled": sig_enabled_w.value,
            "grid_len": grid_len_w.value,
            "grid_z_cells": grid_z_cells_w.value,
            "grid_size": grid_size_w.value,
            "boundary_condition": boundary_w.value,
            "signals": signals,
        }

        parameters["sbol_mapping"] = {"ignore_component_ids": _current_ignore_ids()}
        return parameters

    def _on_save(_btn):
        _collect_parameters()
        with output_area:
            clear_output(wait=True)
            print("Parameters saved.")

    save_btn = widgets.Button(description="Save parameters", button_style="info")
    save_btn.on_click(_on_save)

    def _on_generate(_btn):
        _collect_parameters()
        with output_area:
            clear_output(wait=True)
            try:
                (proteins_, modules_, interactions_,
                 component_map_, ed_chemicals_) = parse_json(
                    sbol_data, ignore_ids=parameters["sbol_mapping"]["ignore_component_ids"]
                )
                script = generate_script(
                    proteins_, modules_, interactions_, component_map_,
                    parameters, ed_chemicals_,
                )
                script_preview_w.value = script
                with open(output_path_w.value, "w") as f:
                    f.write(script)
                print(f"Script generated and written to {output_path_w.value}")
                if parameters["signaling"]["enabled"] and not parameters["signaling"]["signals"]:
                    print(
                        "Note: signalling is enabled but no signals are defined — "
                        "GridDiffusion will be set up with 0 signals."
                    )
                elif not parameters["signaling"]["enabled"] and diffusible_signals:
                    print(
                        "Note: SBOL contains diffusible signal(s) "
                        f"{sorted(diffusible_signals)} but signalling is disabled — "
                        "tick 'Enable inter-cell signalling' to simulate them."
                    )
            except Exception as exc:
                script_preview_w.value = ""
                print(f"Error generating script: {exc}")

    generate_btn = widgets.Button(description="Generate CellModeller script", button_style="success")
    generate_btn.on_click(_on_generate)

    actions_box = _section(
        "Generate",
        widgets.HBox([save_btn, _spacer(), generate_btn]),
        output_path_w,
        output_area,
        widgets.HTML(value="<b>Script preview</b>"),
        script_preview_w,
    )

    # Populate `parameters` once up front so it's usable even if the form
    # is never touched before being read.
    _collect_parameters()

    form = widgets.VBox([
        simulation_box, cell_types_box, kinetics_box, signaling_box,
        advanced_box, actions_box,
    ])
    display(form)

    return parameters, topology
