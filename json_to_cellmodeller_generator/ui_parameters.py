import ipywidgets as widgets
from IPython.display import display, clear_output
from json_to_gro_generator.gro_file_compiler import extract_genes_and_qs_actions, extract_ncb_production_genes_and_actions, extract_biochemical_reactions, normalize_signal_name
import numpy as np

def display_form(sbol_data):
  """
  Generates and displays an interactive form using ipywidgets to gather simulation parameters
  based on SBOL data.

  Args:
    sbol_data (dict): A dictionary containing SBOL data, expected to have keys like
                      "interactions", "hierarchy", "components", and "ED".

  Returns:
    tuple: A tuple containing:
      - parameters (dict): The collected parameters from the form.
      - all_genes_for_gro_generation (list): A list of all genes (detected and auxiliary) for GRO generation.
      - qs_actions_map (dict): A map of quorum sensing actions.
      - biochemical_reactions_data (list): Data on biochemical reactions.
  """
  genes_detected, qs_actions_map = extract_genes_and_qs_actions(
    sbol_data.get("interactions", []),
    sbol_data.get("hierarchy", {}),
    sbol_data.get("components", []),
    sbol_data.get("ED", [])
  )

  # Create a modifiable copy of qs_actions_map to be potentially updated by other functions
  qs_actions_map_modifiable = qs_actions_map.copy()

  auxiliary_genes_ncb, ncb_emission_actions_info = extract_ncb_production_genes_and_actions(
    sbol_data.get("interactions", []),
    sbol_data.get("ED", []),
    sbol_data.get("components", []),
    sbol_data.get("hierarchy", {}),
    qs_actions_map_modifiable,
  )

  biochemical_reactions_data = extract_biochemical_reactions(
    sbol_data.get("interactions", []),
    sbol_data.get("ED", []),
    sbol_data.get("hierarchy", {}),
    sbol_data.get("components", [])
  )
  # Update qs_actions_map with any modifications made to the copy
  qs_actions_map = qs_actions_map_modifiable
  formatted_auxiliary_gene_names_ncb = [g_aux["name"] for g_aux in auxiliary_genes_ncb]
  available_gene_names_ui = [gene["name"].strip('"') for gene in genes_detected]
  ed_list = sbol_data.get("ED", [])

  # Prepare data for S1->S2 enzymatic reactions for UI and signal generation logic
  s1_s2_enzymatic_reactions_for_ui = []
  if biochemical_reactions_data:
    s1_s2_enzymatic_reactions_for_ui = [
      r for r in biochemical_reactions_data
      if r.get("type") == "enzymatic_conversion_S1_to_S2"
    ]

  signals_produced_by_biochem_conversion = set()
  for reaction_data_item in s1_s2_enzymatic_reactions_for_ui:
    if "product_signal" in reaction_data_item and "original_name" in reaction_data_item["product_signal"]:
      signals_produced_by_biochem_conversion.add(reaction_data_item["product_signal"]["original_name"])

  signals_emitted_by_ncb = set()
  if ncb_emission_actions_info:
    signals_emitted_by_ncb = {
      action_info["emitted_signal_original_name"]
      for action_info in ncb_emission_actions_info
      if action_info.get("type") == "ncb_emission"
    }

  content_widget_layout_small = widgets.Layout(width='auto', flex='1 1 100px', min_width='80px')
  content_widget_layout_medium = widgets.Layout(width='auto', flex='1 1 150px', min_width='120px')
  content_widget_layout_large = widgets.Layout(width='auto', flex='1 1 200px', min_width='180px')
  spacer_layout = widgets.Layout(flex='0.5 1 15px', min_width='10px')

  def create_spacer():
    """
    Creates a spacer widget.

    Returns:
      widgets.Label: An empty label widget used as a spacer.
    """
    return widgets.Label(value="", layout=spacer_layout)

  # --- Global Parameters Section ---
  dt_minutes_widget = widgets.FloatText(description="Timestep (dt, minutes):", value=0.1,style={'description_width': 'initial'}, layout=content_widget_layout_large)
  population_max_widget = widgets.IntText(description="Population max:", value=20000, style={'description_width': 'initial'}, layout=content_widget_layout_large)
  doubling_time_widget = widgets.FloatText(description="Cell Doubling Time (minutes):", value=20.0,style={'description_width': 'initial'}, layout=content_widget_layout_large)
  
  global_params_hbox = widgets.HBox(
    [dt_minutes_widget, create_spacer(), doubling_time_widget, create_spacer(), population_max_widget],
    layout=widgets.Layout(width='100%', justify_content='flex-start', align_items='center')
  )
  global_box = widgets.VBox([
    widgets.HTML(value="<h2>Global Simulation Parameters</h2>"),
    global_params_hbox
  ], layout=widgets.Layout(border='1px solid lightgray', margin='10px 0', padding='10px'))

  # Gene Parameters Section 
  gene_widgets = {}
  all_genes_boxes = []
  for gene in genes_detected:
    gene_name_key = gene["name"].strip('"')
    proteins_for_gene = sorted(list({p.strip('"') for p in gene["proteins"]}))

    gene_box_content = []
    gene_info_html = f"<b>Gene: {gene_name_key}</b> (Produces: {', '.join(proteins_for_gene)})"
    gene_box_content.append(widgets.HTML(value=gene_info_html))
    protein_params_widgets = {}
    if len(proteins_for_gene) == 1:
      protein_id = proteins_for_gene[0]
      gene_box_content.append(widgets.Label("Timing (Activation/Degradation, minutes):"))
      
      act_time_w = widgets.FloatText(description="Act. Time:", value=15.0, layout=content_widget_layout_medium, style={'description_width': 'initial'})
      act_var_w = widgets.FloatText(description="Act. Var.:", value=1.0, layout=content_widget_layout_small, style={'description_width': 'initial'})
      deg_time_w = widgets.FloatText(description="Deg. Time:", value=20.0, layout=content_widget_layout_medium, style={'description_width': 'initial'})
      deg_var_w = widgets.FloatText(description="Deg. Var.:", value=1.0, layout=content_widget_layout_small, style={'description_width': 'initial'})

      protein_params_widgets[protein_id] = {
        "act_time": act_time_w, "act_var": act_var_w,
        "deg_time": deg_time_w, "deg_var": deg_var_w
      }
      protein_hbox = widgets.HBox([act_time_w, create_spacer(), act_var_w, create_spacer(), deg_time_w, create_spacer(), deg_var_w], layout=widgets.Layout(width='100%', justify_content='flex-start', align_items='center'))
      gene_box_content.append(protein_hbox)
    else:
      gene_box_content.append(widgets.Label("Timing (Activation/Degradation, minutes):"))
      for protein_id in proteins_for_gene:
        protein_header = widgets.HTML(value=f"<b>Protein: {protein_id}</b>", layout=widgets.Layout(margin='5px 0 2px 0'))
        gene_box_content.append(protein_header)
        
        act_time_w = widgets.FloatText(description="Act. Time:", value=15.0, layout=content_widget_layout_medium, style={'description_width': 'initial'})
        act_var_w = widgets.FloatText(description="Act. Var.:", value=1.0, layout=content_widget_layout_small, style={'description_width': 'initial'})
        deg_time_w = widgets.FloatText(description="Deg. Time:", value=20.0, layout=content_widget_layout_medium, style={'description_width': 'initial'})
        deg_var_w = widgets.FloatText(description="Deg. Var.:", value=1.0, layout=content_widget_layout_small, style={'description_width': 'initial'})

        protein_params_widgets[protein_id] = {
          "act_time": act_time_w, "act_var": act_var_w,
          "deg_time": deg_time_w, "deg_var": deg_var_w
        }
        protein_hbox = widgets.HBox([act_time_w, create_spacer(), act_var_w, create_spacer(), deg_time_w, create_spacer(), deg_var_w], layout=widgets.Layout(width='100%', justify_content='flex-start', align_items='center'))
        gene_box_content.append(protein_hbox)

    gene_box_content.append(widgets.Label("Noise Parameters:"))
    noise_widgets = {
      "toOn": widgets.FloatSlider(description="P(noise ON):", min=0.0, max=1.0, step=0.01, value=0.0, layout=content_widget_layout_medium, style={'description_width': 'initial'}),
      "toOff": widgets.FloatSlider(description="P(noise Off):", min=0.0, max=1.0, step=0.01, value=0.0, layout=content_widget_layout_medium, style={'description_width': 'initial'}),
      "noise_time": widgets.FloatText(description="Noise Act. Time:", value=100.0, layout=content_widget_layout_medium, style={'description_width': 'initial'})
    }
    noise_params_hbox = widgets.HBox([
      noise_widgets["noise_time"], create_spacer(),
      noise_widgets["toOn"], create_spacer(),
      noise_widgets["toOff"]
    ], layout=widgets.Layout(width='100%', justify_content='flex-start', align_items='center'))
    gene_box_content.append(noise_params_hbox)

    gene_widgets[gene_name_key] = {
      "protein_params": protein_params_widgets,
      "proteins_order": proteins_for_gene,
      "noise_params": noise_widgets
    }
    
    gene_box = widgets.VBox(gene_box_content, layout=widgets.Layout(border='1px dashed lightgray', margin='5px 0', padding='5px', width='100%'))
    all_genes_boxes.append(gene_box)

  genes_section_box = widgets.VBox([widgets.HTML(value="<h2>Gene Parameters</h2>")] + all_genes_boxes,
                                   layout=widgets.Layout(border='1px solid lightgray', margin='10px 0', padding='10px', width='100%'))
  
  # Plasmid Definition Section 
  plasmids_ui_list = []
  plasmid_container = widgets.VBox([])
  ecoli_population_ui_list = []
  conjugation_plasmids_checkbox_vbox = widgets.VBox([])
  conjugation_checkboxes_dict = {}
  dynamic_conjugation_rate_widgets = {}
  conjugation_rates_container = widgets.VBox([])

  def get_defined_plasmid_names():
    """
    Retrieves the names of plasmids defined in the UI.

    Returns:
      list: A list of defined plasmid names (strings).
    """
    defined_names = []
    for plasmid_entry in plasmids_ui_list:
      plasmid_w = plasmid_entry["widgets"]
      try:
        name_widget = plasmid_w["name_input"]
        if isinstance(name_widget, widgets.Text) and name_widget.value.strip():
          defined_names.append(name_widget.value.strip())
      except Exception:
        continue
    return defined_names

  def update_conjugation_rates_ui(change_event):
    """
    Updates the UI for conjugation rate inputs based on selected plasmids.
    """
    nonlocal dynamic_conjugation_rate_widgets
    selected_plasmids_for_conj = [name for name, cb in conjugation_checkboxes_dict.items() if cb.value]
    new_rate_widgets_for_conj = {}
    widgets_to_display_for_conj = []
    
    for plasmid_name_conj_str in selected_plasmids_for_conj:
      if plasmid_name_conj_str in dynamic_conjugation_rate_widgets:
        rate_widget_conj_item = dynamic_conjugation_rate_widgets[plasmid_name_conj_str]
      else:
        rate_widget_conj_item = widgets.FloatText(
          value=1.0, 
          description=f"Rate (events/doubling time) for {plasmid_name_conj_str}:",
          style={'description_width': 'initial'},
          layout=widgets.Layout(width='auto', min_width='250px', max_width='100%')
        )
      new_rate_widgets_for_conj[plasmid_name_conj_str] = rate_widget_conj_item
      widgets_to_display_for_conj.append(rate_widget_conj_item)

    dynamic_conjugation_rate_widgets = new_rate_widgets_for_conj
    conjugation_rates_container.children = tuple(widgets_to_display_for_conj)
    
  def update_conjugation_plasmid_options():
    """
    Updates the conjugation plasmid checkbox options based on currently defined plasmids.
    """
    nonlocal conjugation_checkboxes_dict, dynamic_conjugation_rate_widgets
    defined_plasmid_names = get_defined_plasmid_names()
    current_selection = {name for name, cb in conjugation_checkboxes_dict.items() if cb.value}
    
    new_checkboxes_list = []
    new_checkbox_dict = {}

    for name in defined_plasmid_names:
      cb = widgets.Checkbox(description=name, value=(name in current_selection), indent=False, layout=widgets.Layout(width='auto'))
      cb.observe(update_conjugation_rates_ui, names='value')
      new_checkboxes_list.append(cb)
      new_checkbox_dict[name] = cb

    conjugation_checkboxes_dict = new_checkbox_dict
    conjugation_plasmids_checkbox_vbox.children = tuple(new_checkboxes_list)
    update_conjugation_rates_ui(None) 

  def update_ecoli_plasmid_options():
    """
    Updates plasmid checkbox options for all defined E. coli populations.
    """
    plasmid_names = get_defined_plasmid_names()
    for population_entry in ecoli_population_ui_list:
      pop_widgets = population_entry["widgets"]
      vbox_widget = pop_widgets.get("plasmids_vbox")
      current_cb_dict = pop_widgets.get("plasmids_checkboxes", {})
      
      if vbox_widget:
        current_selection = {name for name, cb in current_cb_dict.items() if cb.value}
        new_cbs = []
        new_dict = {}
        for name in plasmid_names:
          cb = widgets.Checkbox(description=name, value=(name in current_selection), indent=False, layout=widgets.Layout(width='auto'))
          new_cbs.append(cb)
          new_dict[name] = cb
        
        vbox_widget.children = tuple(new_cbs)
        pop_widgets["plasmids_checkboxes"] = new_dict

  def on_plasmid_name_change(change): 
    """
    Callback function triggered when a plasmid name input changes.
    Updates options for E.coli plasmid selectors and conjugation plasmid selector.

    Args:
      change (dict): The change event object from the widget.
    """
    update_ecoli_plasmid_options()
    update_conjugation_plasmid_options()

  def add_plasmid_ui(button_event=None):
    """
    Adds a new set of UI elements for defining a plasmid.

    Args:
      button_event (widgets.Button, optional): The button click event. Defaults to None.
    """
    nonlocal plasmids_ui_list
    plasmid_idx = len(plasmids_ui_list) + 1
    plasmid_name_input_widget = widgets.Text(description=f"Plasmid {plasmid_idx} Name:", value=f"p{plasmid_idx}", style={'description_width': 'initial'}, layout=widgets.Layout(flex='1 1 auto', min_width='150px'))
    plasmid_name_input_widget.observe(on_plasmid_name_change, names='value')
    remove_plasmid_button_widget = widgets.Button(description="X", button_style='danger', layout=widgets.Layout(width='35px', height='35px', margin='0 0 0 5px', display='flex', align_items='center', justify_content='center'))
    name_and_remove_hbox_widget = widgets.HBox([plasmid_name_input_widget, remove_plasmid_button_widget], layout=widgets.Layout(width='100%', align_items='center'))

    gene_selection_label_widget = widgets.Label("Assign Genes:")
    gene_checkbox_vbox = widgets.VBox([])
    gene_checkboxes_dict = {}
    for gene_name in available_gene_names_ui:
        cb = widgets.Checkbox(description=gene_name, value=False, indent=False)
        gene_checkbox_vbox.children += (cb,)
        gene_checkboxes_dict[gene_name] = cb

    plasmid_card_vbox_widget = widgets.VBox([name_and_remove_hbox_widget, gene_selection_label_widget, gene_checkbox_vbox], layout=widgets.Layout(border='1px dashed lightgray', margin='5px', padding='10px', width='auto'))
    current_plasmid_widgets_dict = {"name_input": plasmid_name_input_widget, "gene_checkboxes": gene_checkboxes_dict}
    plasmid_entry_dict = {"ui_box": plasmid_card_vbox_widget, "widgets": current_plasmid_widgets_dict}
    plasmids_ui_list.append(plasmid_entry_dict)

    def on_remove_plasmid_clicked(button_inner_event, entry_to_remove=plasmid_entry_dict):
      if entry_to_remove in plasmids_ui_list:
        plasmids_ui_list.remove(entry_to_remove)
        update_plasmid_display()
    remove_plasmid_button_widget.on_click(lambda btn_event_remove, entry=plasmid_entry_dict: on_remove_plasmid_clicked(btn_event_remove, entry_to_remove=entry))
    update_plasmid_display()

  def update_plasmid_display():
    """
    Updates the displayed list of plasmid UI elements and refreshes related selectors.
    """
    plasmid_container.children = tuple([entry["ui_box"] for entry in plasmids_ui_list])
    for idx, entry in enumerate(plasmids_ui_list):
      try:
        entry["widgets"]["name_input"].description = f"Plasmid {idx + 1} Name:"
      except Exception: 
        pass 
    update_ecoli_plasmid_options()
    update_conjugation_plasmid_options()

  add_plasmid_button = widgets.Button(description="Add Additional Plasmid", layout={'width': 'auto'})
  add_plasmid_button.on_click(add_plasmid_ui)
  plasmid_definition_instruction = widgets.Label("Define plasmids and assign detected genes. If none defined, specific defaults apply during save.")
  plasmids_section_box = widgets.VBox([
    widgets.HTML(value="<h2>Plasmid Definition</h2>"),
    plasmid_definition_instruction, add_plasmid_button, plasmid_container
  ], layout=widgets.Layout(border='1px solid lightgray', margin='10px 0', padding='10px'))

  if not plasmids_ui_list: # Add one plasmid UI by default if list is empty
    add_plasmid_ui(None)

  # --- E.coli Populations Section ---
  ecoli_populations_container = widgets.VBox([])

  def add_ecoli_population_ui(button_event=None):
    """
    Adds a new set of UI elements for defining an E. coli population group.

    Args:
      button_event (widgets.Button, optional): The button click event. Defaults to None.
    """
    nonlocal ecoli_population_ui_list
    population_idx = len(ecoli_population_ui_list) + 1
    title_widget_html = widgets.HTML(value=f"<h4>Population Group {population_idx}</h4>", layout=widgets.Layout(flex='1 1 auto', margin='0 0 5px 0'))
    remove_population_button = widgets.Button(description="X", button_style='danger', layout=widgets.Layout(width='35px', height='35px', display='flex', align_items='center', justify_content='center'))
    title_and_remove_hbox_layout = widgets.HBox([title_widget_html, remove_population_button], layout=widgets.Layout(justify_content='space-between', width='100%', align_items='center', margin='0 0 5px 0'))

    num_ecolis_widget_int = widgets.IntText(description="Number of Cells:", value=100, style={'description_width': 'initial'}, layout=content_widget_layout_medium)
    center_x_widget_float = widgets.FloatText(description="Center X:", value=0.0, style={'description_width': 'initial'}, layout=content_widget_layout_small)
    center_y_widget_float = widgets.FloatText(description="Center Y:", value=0.0, style={'description_width': 'initial'}, layout=content_widget_layout_small)
    radius_widget_float = widgets.FloatText(description="Radius:", value=100.0, style={'description_width': 'initial'}, layout=content_widget_layout_small)
    population_params_hbox = widgets.HBox([num_ecolis_widget_int, create_spacer(), center_x_widget_float, create_spacer(), center_y_widget_float, create_spacer(), radius_widget_float], layout=widgets.Layout(width='100%', justify_content='flex-start', align_items='center'))

    plasmid_selector_label_ecoli = widgets.Label("Select Plasmids:")
    plasmid_checkbox_vbox_ecoli = widgets.VBox([], layout={'width': '100%', 'min_width': '230px', 'margin':'5px 0 0 0'})
    plasmid_checkboxes_dict_ecoli = {} # Se llenará en update_ecoli_plasmid_checkboxes

    population_card_vbox_layout = widgets.VBox([title_and_remove_hbox_layout, population_params_hbox, plasmid_selector_label_ecoli, plasmid_checkbox_vbox_ecoli], layout=widgets.Layout(border='1px dashed lightgray', margin='10px 0', padding='5px 10px 10px 10px', align_items='stretch', width='auto'))
    current_population_widgets_dict = {
      "num": num_ecolis_widget_int, "x": center_x_widget_float, "y": center_y_widget_float, "radius": radius_widget_float, 
      "plasmids_vbox": plasmid_checkbox_vbox_ecoli, 
      "plasmids_checkboxes": plasmid_checkboxes_dict_ecoli, 
      "title_widget": title_widget_html
    }
    population_entry_item = {"ui_box": population_card_vbox_layout, "widgets": current_population_widgets_dict}
    ecoli_population_ui_list.append(population_entry_item)

    def on_remove_population_clicked(btn_inner_event, entry_to_remove=population_entry_item):
      if entry_to_remove in ecoli_population_ui_list:
        ecoli_population_ui_list.remove(entry_to_remove)
        update_ecoli_populations_display()
    remove_population_button.on_click(lambda btn_event_r, entry=population_entry_item: on_remove_population_clicked(btn_event_r, entry_to_remove=entry))

    update_ecoli_populations_display()
    update_ecoli_plasmid_options() 

  def update_ecoli_populations_display():
    """
    Updates the displayed list of E. coli population UI elements.
    """
    ecoli_populations_container.children = tuple([entry["ui_box"] for entry in ecoli_population_ui_list])
    for idx, entry in enumerate(ecoli_population_ui_list):
      try:
        entry["widgets"]["title_widget"].value = f"<h4>Population Group {idx + 1}</h4>"
      except Exception:
        pass 

  add_ecoli_population_button = widgets.Button(description="Add Additional Population", layout={'width': 'auto'})
  add_ecoli_population_button.on_click(add_ecoli_population_ui)
  initial_ecoli_instruction_html = widgets.HTML(value="Define one or more groups of E. coli cells to initialize the simulation.")
  initial_ecoli_section_box = widgets.VBox([widgets.HTML(value="<h2>Initial Cell Populations</h2>"), initial_ecoli_instruction_html, add_ecoli_population_button, ecoli_populations_container], layout=widgets.Layout(border='1px solid lightgray', margin='10px 0', padding='10px'))

  if not ecoli_population_ui_list: # Add one population UI by default
    add_ecoli_population_ui(None)

  # Signal Parameters Section 
  signal_widgets_dict = {}
  signal_boxes_list = []
  chemical_names_in_ed = {ed_item.get("name") for ed_item in ed_list if ed_item.get("type") == "Simple chemical"}

  def create_add_point_handler(point_list_reference, points_container_reference):
    """
    Creates a handler function for adding new signal point UI elements.

    Args:
      point_list_reference (list): A reference to the list holding the point UI HBox widgets.
      points_container_reference (widgets.VBox): The VBox container to display the point HBoxes.

    Returns:
      function: The actual handler function to be called on button click.
    """
    def actual_add_handler(button_ignored): 
      """Handles adding a new signal point UI."""
      text_layout_cfg = widgets.Layout(width='auto', flex='1 1 60px', min_width='50px')
      point_x_widget = widgets.FloatText(description="X:", value=0.0, layout=text_layout_cfg, style={'description_width': 'auto'})
      point_y_widget = widgets.FloatText(description="Y:", value=0.0, layout=text_layout_cfg, style={'description_width': 'auto'})
      point_conc_widget = widgets.FloatText(description="Conc:", value=15.0, step=1.0, layout=text_layout_cfg, style={'description_width': 'auto'})
      input_fields_hbox_widget = widgets.HBox([point_x_widget, create_spacer(), point_y_widget, create_spacer(), point_conc_widget], layout=widgets.Layout(gap='10px', flex='3 1 auto'))
      point_constant_emission_cb = widgets.Checkbox(description="Constant Emission", value=False, indent=False, tooltip="Emit this concentration from this point at every timestep", layout=widgets.Layout(width='auto', margin='0 10px 0 15px', flex='0 0 auto'))
      remove_button = widgets.Button(description="X", button_style='danger', layout={'width': '40px', 'margin': '0 0 0 5px', 'flex': '0 0 auto'})
      current_point_hbox_widget = widgets.HBox([input_fields_hbox_widget, create_spacer(), point_constant_emission_cb, create_spacer(), remove_button], layout=widgets.Layout(border='1px solid lightgray', padding='5px', width='100%', align_items='center', justify_content='space-between', margin='2px 0'))

      point_list_reference.append(current_point_hbox_widget)
      points_container_reference.children = tuple(point_list_reference)

      def remove_point_ui_handler(btn_remove, box_to_remove=current_point_hbox_widget):
        """Removes a signal point UI."""
        if box_to_remove in point_list_reference:
          point_list_reference.remove(box_to_remove)
          points_container_reference.children = tuple(point_list_reference)
      remove_button.on_click(remove_point_ui_handler)
    return actual_add_handler

  if chemical_names_in_ed:
    signal_section_title_html = widgets.HTML(value="<h2>Signal Parameters</h2>")
    signal_boxes_list.append(signal_section_title_html)

    for chem_name_original_str in sorted(list(chemical_names_in_ed)):
      normalized_chem_name_for_qs = normalize_signal_name(chem_name_original_str)
      is_sensed_as_qs_effector = normalized_chem_name_for_qs in qs_actions_map
      is_emitted_via_ncb = chem_name_original_str in signals_emitted_by_ncb
      is_product_of_biochem_conv = chem_name_original_str in signals_produced_by_biochem_conversion
      is_internally_generated = is_emitted_via_ncb or is_product_of_biochem_conv

      kdiff_slider = widgets.FloatSlider(description="Diffusion Rate (kdiff):", layout=content_widget_layout_large, style={'description_width': 'initial'}, min=0.0, max=1.0, step=0.001, value=0.1, readout_format='.3f')
      kdeg_slider = widgets.FloatSlider(description="Degradation Rate (kdeg):", layout=content_widget_layout_large, style={'description_width': 'initial'}, min=0.0, max=0.1, step=0.0001, value=0.001, readout_format='.4f')

      qs_symbol_dropdown_widget = widgets.Dropdown(description="QS Operator:", layout=content_widget_layout_medium, options=[">", "<"], value=">", style={'description_width': 'initial'}, disabled=not is_sensed_as_qs_effector)
      qs_threshold_input_widget = widgets.FloatText(description="QS Threshold:", layout=content_widget_layout_medium, value=0.3, style={'description_width': 'initial'}, disabled=not is_sensed_as_qs_effector)

      current_signal_initial_points_container = widgets.VBox([])
      current_signal_point_widgets_list = []
      add_initial_point_button_widget = None

      add_point_handler_fn = create_add_point_handler(current_signal_point_widgets_list, current_signal_initial_points_container)

      if not is_internally_generated:
        add_initial_point_button_widget = widgets.Button(description="Add Additional Initial Signal Point", layout={'width': 'auto'})
        add_initial_point_button_widget.on_click(add_point_handler_fn)
        add_point_handler_fn(None) # Add one point by default for externally added signals
        current_signal_initial_points_container.layout.margin = '10px 0 0 0'

      signal_widgets_dict[chem_name_original_str] = {
        "kdiff": kdiff_slider, "kdeg": kdeg_slider,
        "Symbol_getQS": qs_symbol_dropdown_widget, "Threshold_getQS": qs_threshold_input_widget,
        "initial_points_container": current_signal_initial_points_container,
        "point_widgets_list": current_signal_point_widgets_list,
        "add_point_button": add_initial_point_button_widget
      }

      diffusion_degradation_hbox = widgets.HBox(
        [kdiff_slider, create_spacer(), kdeg_slider],
        layout=widgets.Layout(width='100%', justify_content='flex-start', align_items='center')
      )
      
      title_suffix = "" # Sufijo para el título
      
      if is_emitted_via_ncb:
        ncb_action_detail_item = next((act for act in ncb_emission_actions_info if act["emitted_signal_original_name"] == chem_name_original_str), None)
        if ncb_action_detail_item:
          reactants_str_list = ", ".join(ncb_action_detail_item.get("reactants_original_names", []))
          title_suffix = f" <i>(Produced by NCB of: <b>[{reactants_str_list}]</b>)</i>"
      elif is_product_of_biochem_conv:
        reaction_details_for_product = next((r for r in s1_s2_enzymatic_reactions_for_ui if r["product_signal"]["original_name"] == chem_name_original_str), None)
        if reaction_details_for_product:
          substrate_name_str = reaction_details_for_product["reactant_signal"]["original_name"]
          enzyme_name_str = reaction_details_for_product["enzyme"]["original_name"]
          title_suffix = f" <i>(Generated by conversion of <b>{substrate_name_str}</b> via <b>{enzyme_name_str}</b>)</i>"
        else:
           title_suffix = " <i>(Product of biochemical conversion)</i>"
      
      current_signal_box_content = [
        widgets.HTML(value=f"<b>Signal: {chem_name_original_str}</b>{title_suffix}"), 
        diffusion_degradation_hbox,
      ]

      if is_sensed_as_qs_effector:
        qs_protein_name_display = qs_actions_map.get(normalized_chem_name_for_qs, ("N/A",))[0]
        qs_hbox_layout = widgets.HBox(
          [qs_symbol_dropdown_widget, create_spacer(), qs_threshold_input_widget],
          layout=widgets.Layout(width='100%', justify_content='flex-start', align_items='center')
        )
        current_signal_box_content.extend([
          widgets.HTML(value=f"<i>This signal is sensed/is an effector, controlling protein <b>{qs_protein_name_display}</b>:</i>"),
          qs_hbox_layout
        ])
      if not is_internally_generated:
        current_signal_box_content.append(widgets.HTML(value="<b>Initial Signal Placement:</b>"))
        if add_initial_point_button_widget:
          current_signal_box_content.append(add_initial_point_button_widget)
        current_signal_box_content.append(current_signal_initial_points_container)

      individual_signal_vbox_layout = widgets.VBox(current_signal_box_content, layout=widgets.Layout(border='1px dashed lightgray', margin='5px 0', padding='10px', align_items='stretch'))
      signal_boxes_list.append(individual_signal_vbox_layout)

  signals_section_box_widget = widgets.VBox() 
  if signal_boxes_list: 
    signals_layout_cfg = widgets.Layout(border='1px solid lightgray', margin='10px 0', padding='10px')
    signals_section_box_widget = widgets.VBox(signal_boxes_list, layout=signals_layout_cfg)

  # Biochemical Reaction Parameters Section
  biochem_conversion_widgets_dict = {}
  biochem_conversion_boxes_ui_list = []

  if s1_s2_enzymatic_reactions_for_ui:
    biochem_conversion_section_title_html = widgets.HTML(
      value="<h2>Biochemical Conversion Parameters</h2>"
    )
    biochem_conversion_boxes_ui_list.append(biochem_conversion_section_title_html)

    for idx, reaction_data_conv in enumerate(s1_s2_enzymatic_reactions_for_ui):
      enzyme_display_name_str = reaction_data_conv["enzyme"]["original_name"]
      enzyme_gro_id_display_str = reaction_data_conv["enzyme"]["gro_id"] # GRO Action
      s1_display_name_str = reaction_data_conv["reactant_signal"]["original_name"]
      s2_display_name_str = reaction_data_conv["product_signal"]["original_name"]
      reaction_key_ui_str = reaction_data_conv.get("sbol_interaction_id", f"{s1_display_name_str}_to_{s2_display_name_str}_by_{enzyme_display_name_str}")
      ui_label_html_str = (f"<b>Reaction {idx+1}: {s1_display_name_str} → {s2_display_name_str}</b> <i> Catalyzed by: {enzyme_display_name_str} (GRO ID: {enzyme_gro_id_display_str})</i>")

      conversion_rate_widget = widgets.FloatText(description="Conversion Rate (S1 Abs & S2 Em, units/dt):", value=1.0, step=0.1, tooltip="Units of S1 absorbed and S2 emitted per cell per dt", style={'description_width': 'initial'}, layout=content_widget_layout_large)
      absorption_type_widget = widgets.Dropdown(description="S1 Abs. Type:", options=["exact", "area", "average"], value="area", style={'description_width': 'initial'}, layout=content_widget_layout_medium)
      emission_type_widget_conv = widgets.Dropdown(description="S2 Em. Type:", options=["exact", "area", "average"], value="area", style={'description_width': 'initial'}, layout=content_widget_layout_medium)

      biochem_conversion_widgets_dict[reaction_key_ui_str] = {"rate": conversion_rate_widget, "absorption_type": absorption_type_widget, "emission_type": emission_type_widget_conv}
      all_conversion_params_hbox_layout = widgets.HBox(
        [conversion_rate_widget, create_spacer(), absorption_type_widget, create_spacer(), emission_type_widget_conv],
        layout=widgets.Layout(width='100%', justify_content='flex-start', align_items='center')
      )
      single_reaction_box_ui_layout = widgets.VBox([widgets.HTML(value=ui_label_html_str), all_conversion_params_hbox_layout], layout=widgets.Layout(border='1px dashed lightgray', margin='10px 0', padding='10px'))
      biochem_conversion_boxes_ui_list.append(single_reaction_box_ui_layout)

  biochem_conversion_section_box_ui_widget = widgets.VBox(
    biochem_conversion_boxes_ui_list,
    layout=widgets.Layout(border='1px solid lightgray', margin='10px 0', padding='10px', width='auto',
                          # Hide if no reactions to show
                          display='block' if s1_s2_enzymatic_reactions_for_ui else 'none')
  )

  # Conjugation Parameters Section
  conjugation_checkbox = widgets.Checkbox(description="Enable Bacterial Conjugation", value=False, indent=False)
  conjugation_plasmids_label = widgets.Label("Select Plasmids for Conjugation:")

  conjugation_details_box = widgets.VBox(
    [conjugation_plasmids_label, conjugation_plasmids_checkbox_vbox, conjugation_rates_container],
    layout=widgets.Layout(display='none', border='1px dashed lightgray', padding='5px', margin='5px', align_items='stretch', width='auto')
  )

  def toggle_conjugation_details(change_event):
    """
    Shows or hides the conjugation details UI based on the checkbox state.

    Args:
      change_event (dict): The change event object from the checkbox.
    """
    nonlocal dynamic_conjugation_rate_widgets
    is_enabled = change_event["new"]
    conjugation_details_box.layout.display = 'block' if is_enabled else 'none'
    for cb in conjugation_checkboxes_dict.values():
      cb.disabled = not is_enabled
      if not is_enabled:
        cb.value = False
    
    if is_enabled:
      update_conjugation_plasmid_options() 
    else:
      dynamic_conjugation_rate_widgets.clear()
      conjugation_rates_container.children = []

  conjugation_checkbox.observe(toggle_conjugation_details, names='value')
  conjugation_section_box = widgets.VBox([widgets.HTML(value="<h2>Bacterial Conjugation Parameters</h2>"), conjugation_checkbox, conjugation_details_box], layout=widgets.Layout(border='1px solid lightgray', margin='10px 0', padding='10px'))

  # NCB Emission Parameters Section
  ncb_emission_widgets_dict = {}
  ncb_emission_boxes_list = []

  # Collect unique NCB emitted signals and their conditioning proteins for UI generation
  unique_ncb_emitted_signals_info = {}
  for action_info_item_ui in ncb_emission_actions_info:
    if action_info_item_ui["type"] == "ncb_emission":
      original_signal_name = action_info_item_ui["emitted_signal_original_name"]
      if original_signal_name not in unique_ncb_emitted_signals_info:
        unique_ncb_emitted_signals_info[original_signal_name] = {"condition_protein": action_info_item_ui["condition_protein"]}

  if unique_ncb_emitted_signals_info:
    ncb_emission_section_title_html = widgets.HTML(value="<h2>Non-Covalent Binding Induced Emission Parameters</h2>")
    ncb_emission_boxes_list.append(ncb_emission_section_title_html)
    for emitted_signal_name, ncb_info_data in unique_ncb_emitted_signals_info.items():
      condition_protein_name = ncb_info_data["condition_protein"]
      ui_label_text_ncb = f"<b>Emitted Signal: {emitted_signal_name}</b> (Controls: {condition_protein_name})"
      concentration_widget_ncb = widgets.FloatText(description="Emission Conc:", value=3.0, step=0.1, style={'description_width': 'initial'}, layout=content_widget_layout_medium)
      emission_type_widget_ncb = widgets.Dropdown(description="Emission Type:", options=["exact", "area", "average"], value="exact", style={'description_width': 'initial'}, layout=content_widget_layout_medium)
      ncb_emission_widgets_dict[emitted_signal_name] = {"concentration": concentration_widget_ncb, "emission_type": emission_type_widget_ncb}
      ncb_params_hbox_layout = widgets.HBox(
        [concentration_widget_ncb, create_spacer(), emission_type_widget_ncb],
        layout=widgets.Layout(width='100%', justify_content='flex-start', align_items='center')
      )
      signal_emission_box_layout = widgets.VBox([widgets.HTML(value=ui_label_text_ncb), ncb_params_hbox_layout], layout=widgets.Layout(border='1px dashed lightgray', margin='5px 0', padding='5px'))
      ncb_emission_boxes_list.append(signal_emission_box_layout)

  ncb_emission_section_box_widget = widgets.VBox(
      ncb_emission_boxes_list,
      layout=widgets.Layout(border='1px solid lightgray', margin='10px 0', padding='10px', display='block' if unique_ncb_emitted_signals_info else 'none')
  )

  # Save Button and Output Area 
  save_button = widgets.Button(description="Save Parameters", button_style='success', icon='save')
  output_area = widgets.Output()
  parameters_dict = {}

  def save_parameters_callback(button_click_event):
    """
    Callback function for the save button. Collects all parameters from the UI,
    validates some, and stores them in the 'parameters_dict'.

    Args:
      button_click_event (widgets.Button): The button click event.
    """
    nonlocal parameters_dict, genes_detected
    parameters_dict.clear()
    user_feedback_messages = []

    # Process dt and Growth Rate
    dt_minutes_for_gro_val = dt_minutes_widget.value
    if dt_minutes_for_gro_val <= 0:
      user_feedback_messages.append("Warning: Timestep (dt) must be > 0. Using default 0.1 minutes.")
      dt_minutes_for_gro_val = 0.1

    doubling_time_min_input_val = doubling_time_widget.value
    growth_rate_for_gro_val = 0.0
    if doubling_time_min_input_val > 0:
      growth_rate_for_gro_val = np.log(2) / doubling_time_min_input_val
    else:
      default_doubling_time_minutes = 20.0
      growth_rate_for_gro_val = np.log(2) / default_doubling_time_minutes
      user_feedback_messages.append(f"Warning: Doubling time ('{doubling_time_min_input_val}') must be > 0. Using default growth rate corresponding to {default_doubling_time_minutes} min doubling time.")

    # Collect plasmid definitions
    plasmid_definitions_to_save = {}
    for plasmid_box_idx, plasmid_ui_entry in enumerate(plasmids_ui_list):
      plasmid_widgets_set = plasmid_ui_entry["widgets"]
      try:
        plasmid_name_str_save = plasmid_widgets_set["name_input"].value.strip()
        if plasmid_name_str_save:
          gene_checkboxes = plasmid_widgets_set.get("gene_checkboxes", {})
          selected_genes_list_ui = [name for name, cb in gene_checkboxes.items() if cb.value]
          genes_for_plasmid_gro_set = {f'"{gene_ui_name}"' for gene_ui_name in selected_genes_list_ui}
          for aux_gene_name_formatted in formatted_auxiliary_gene_names_ncb:
            genes_for_plasmid_gro_set.add(aux_gene_name_formatted)
          if genes_for_plasmid_gro_set: # Only save if there are genes (real or auxiliary)
            plasmid_definitions_to_save[plasmid_name_str_save] = sorted(list(genes_for_plasmid_gro_set))
          else: # No real genes selected, and no auxiliary genes applicable
            user_feedback_messages.append(f"Info: Plasmid '{plasmid_name_str_save}' was defined in the UI but no genes were assigned (genes from the selector are required). This plasmid will be ignored.")
        else: # Plasmid has no name
            user_feedback_messages.append(f"Info: Plasmid UI at position {plasmid_box_idx+1} has no name and will be ignored.")
      except Exception as e: 
        user_feedback_messages.append(f"Error reading UI for plasmid at UI position {plasmid_box_idx+1}: {e}. This plasmid definition will be skipped.")
        continue

    if not plasmid_definitions_to_save:
      message = (
        "Warning: No valid plasmids were configured. "
        "This could be because all plasmid UI entries were left unnamed, had no genes selected from the 'Assign Genes' selector, "
        "or were removed."
        "This will impact E.coli cell setup (if they require plasmids) and conjugation features. "
        "The final .gro file might not include plasmid specifications."
      )
      if not plasmids_ui_list:
        message = (
        "Warning: No plasmid definition sections were found (all UI cards removed). "
        )
      user_feedback_messages.append(message)

    # Collect E.coli populations
    initial_ecoli_populations_data_list = []
    for population_idx_save, population_ui_entry in enumerate(ecoli_population_ui_list):
      population_widgets_set = population_ui_entry["widgets"]
      try:
        num_cells = population_widgets_set["num"].value
        center_x_val = population_widgets_set["x"].value
        center_y_val = population_widgets_set["y"].value
        radius_val = population_widgets_set["radius"].value
        # Plasmids selector can be empty, which is valid (cells with no plasmids)
        plasmid_checkboxes = population_widgets_set.get("plasmids_checkboxes", {})
        selected_plasmids_for_population = [name for name, cb in plasmid_checkboxes.items() if cb.value]
        if num_cells <= 0:
            user_feedback_messages.append(f"Info: E.coli population group {population_idx_save+1} has a non-positive number of cells ({num_cells}) and will be ignored.")
            continue 

        initial_ecoli_populations_data_list.append({
          "num_ecolis": num_cells,
          "center_x": center_x_val,
          "center_y": center_y_val,
          "radius": radius_val,
          "plasmids": selected_plasmids_for_population 
        })
      except Exception as e: 
        user_feedback_messages.append(f"Error reading UI for E.coli population group {population_idx_save+1}: {e}. This population group will be skipped.")
        continue

    if not initial_ecoli_populations_data_list:
      message = (
        "Warning: No valid E. coli population groups were configured. "
        "This could be because all population UI entries had issues (e.g., zero cells) or were removed. "
        "The simulation setup might rely on default population creation in the .gro file, or no cells will be simulated."
      )
      if not ecoli_population_ui_list:
        message = (
          "Warning: No E. coli population group sections were found (all UI cards removed). "
        )
      user_feedback_messages.append(message)

    parameters_dict.update({
      "dt": dt_minutes_for_gro_val,
      "population_max": population_max_widget.value,
      "growth_rate": growth_rate_for_gro_val,
      "gene_parameters": {},
      "signal_parameters": {},
      "plasmid_configuration": {"defined_plasmids": plasmid_definitions_to_save},
      "initial_ecoli_populations": initial_ecoli_populations_data_list,
      "conjugation_parameters": {},
      "ncb_emission_parameters": {},
      "info_for_ncb_emission_actions": ncb_emission_actions_info, 
      "biochemical_conversion_parameters": {},
      "ed_list": ed_list 
    })

    # Save gene parameters
    for gene_key_str, widgets_collection in gene_widgets.items():
      gene_params = {}
      protein_params_w = widgets_collection.get("protein_params", {})

      act_times = []
      act_vars = []
      deg_times = []
      deg_vars = []

      for protein_id in widgets_collection.get("proteins_order", []):
        if protein_id in protein_params_w:
          params = protein_params_w[protein_id]
          act_times.append(params["act_time"].value)
          act_vars.append(params["act_var"].value)
          deg_times.append(params["deg_time"].value)
          deg_vars.append(params["deg_var"].value)

      gene_params["act_times"] = act_times
      gene_params["act_vars"] = act_vars
      gene_params["deg_times"] = deg_times
      gene_params["deg_vars"] = deg_vars

      noise_params_w = widgets_collection.get("noise_params", {})
      if noise_params_w:
        gene_params["toOn"] = noise_params_w["toOn"].value
        gene_params["toOff"] = noise_params_w["toOff"].value
        gene_params["noise_time"] = noise_params_w["noise_time"].value

      parameters_dict["gene_parameters"][gene_key_str] = gene_params

    # Save signal parameters
    for chem_name_orig, widgets_chem_collection in signal_widgets_dict.items():
      params_for_this_signal = {}
      for key_chem_str, widget_chem_item in widgets_chem_collection.items():
        if key_chem_str not in ["point_widgets_list", "initial_points_container", "add_point_button"]:
          params_for_this_signal[key_chem_str] = widget_chem_item.value

      initial_points_data_for_chem = []
      if "point_widgets_list" in widgets_chem_collection and widgets_chem_collection["point_widgets_list"] is not None:
        for point_outer_hbox_widget in widgets_chem_collection["point_widgets_list"]:
          try:
            input_fields_hbox_item = point_outer_hbox_widget.children[0] 
            constant_emission_cb_item = point_outer_hbox_widget.children[2] 

            x_val_chem_point = input_fields_hbox_item.children[0].value
            y_val_chem_point = input_fields_hbox_item.children[2].value 
            conc_val_chem_point = input_fields_hbox_item.children[4].value
            constant_emission_val = constant_emission_cb_item.value

            if conc_val_chem_point > 0: # Only save points with positive concentration
              initial_points_data_for_chem.append({
                "x": x_val_chem_point,
                "y": y_val_chem_point,
                "conc": conc_val_chem_point,
                "constant_emission": constant_emission_val
              })
          except IndexError:
            user_feedback_messages.append(f"Warning: UI structure for a signal point of '{chem_name_orig}' is not as expected. Skipping one point.")
          except Exception as e: 
            user_feedback_messages.append(f"Warning: Error reading data for a specific point of signal '{chem_name_orig}': {e}. Skipping one point.")

      params_for_this_signal["initial_points"] = initial_points_data_for_chem
      parameters_dict["signal_parameters"][chem_name_orig] = params_for_this_signal

    # Save conjugation parameters
    is_conjugation_enabled_val = conjugation_checkbox.value
    conjugation_settings_to_save = {}
    if is_conjugation_enabled_val:
      valid_plasmid_names_for_conj_save = list(plasmid_definitions_to_save.keys())
      selected_plasmids_in_ui_for_conj = [name for name, cb in conjugation_checkboxes_dict.items() if cb.value]
      for plasmid_name_conj_save_str in selected_plasmids_in_ui_for_conj:
        if plasmid_name_conj_save_str not in valid_plasmid_names_for_conj_save: user_feedback_messages.append(f"Info: Plasmid '{plasmid_name_conj_save_str}' ... Ignoring."); continue
        if plasmid_name_conj_save_str in dynamic_conjugation_rate_widgets:
          events_per_doubling_time = dynamic_conjugation_rate_widgets[plasmid_name_conj_save_str].value
          gro_conjugation_rate_val = max(0.0, events_per_doubling_time)
          if events_per_doubling_time < 0: user_feedback_messages.append(f"Warning: Conjugation rate ... Using 0.")
          conjugation_settings_to_save[plasmid_name_conj_save_str] = gro_conjugation_rate_val
        else: conjugation_settings_to_save[plasmid_name_conj_save_str] = 0.0; user_feedback_messages.append(f"Info: Conjugation rate widget for '{plasmid_name_conj_save_str}' ... Defaulting to 0.")
    parameters_dict["conjugation_parameters"] = {"enabled": is_conjugation_enabled_val, "settings": conjugation_settings_to_save}

    # Save NCB emission parameters (key is the original emitted signal name)
    ncb_emissions_data_to_save = {}
    for signal_original_key_str, emission_widget_collection in ncb_emission_widgets_dict.items():
      ncb_emissions_data_to_save[signal_original_key_str] = {
        "concentration": emission_widget_collection["concentration"].value,
        "emission_type": emission_widget_collection["emission_type"].value
      }
    parameters_dict["ncb_emission_parameters"] = ncb_emissions_data_to_save

    # Save Biochemical Reaction parameters
    if biochem_conversion_widgets_dict: 
      for reaction_key_save_str, reaction_widgets_collection in biochem_conversion_widgets_dict.items():
        parameters_dict["biochemical_conversion_parameters"][reaction_key_save_str] = {
          "rate": reaction_widgets_collection["rate"].value,
          "absorption_type": reaction_widgets_collection["absorption_type"].value,
          "emission_type": reaction_widgets_collection["emission_type"].value,
        }

    # Clear output area and show messages
    with output_area:
      clear_output(wait=True)
      for msg in user_feedback_messages:
        print(msg)
      print("Parameters saved successfully.")

    # Update options in case plasmid names changed and affected other selectors
    update_ecoli_plasmid_options()
    update_conjugation_plasmid_options()

  save_button.on_click(save_parameters_callback)

  # Assemble Form
  form_children_list = [global_box, genes_section_box, plasmids_section_box, initial_ecoli_section_box]
  if chemical_names_in_ed: form_children_list.append(signals_section_box_widget)
  if unique_ncb_emitted_signals_info: form_children_list.append(ncb_emission_section_box_widget)
  if s1_s2_enzymatic_reactions_for_ui: form_children_list.append(biochem_conversion_section_box_ui_widget)
  form_children_list.extend([conjugation_section_box, save_button, output_area])
  main_form_widget = widgets.VBox(form_children_list)

  if not plasmids_ui_list: add_plasmid_ui(None)
  if not ecoli_population_ui_list: add_ecoli_population_ui(None)

  # Initial update of plasmid options for selectors
  update_ecoli_plasmid_options()
  update_conjugation_plasmid_options()
  toggle_conjugation_details({'new': conjugation_checkbox.value})

  display(main_form_widget)
  all_genes_for_gro_generation = genes_detected + auxiliary_genes_ncb

  return parameters_dict, all_genes_for_gro_generation, qs_actions_map, biochemical_reactions_data
