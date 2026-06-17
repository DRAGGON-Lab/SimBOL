# CellModeller Parameters Guide
## What you need beyond the SBOL/JSON

The SBOL/JSON file describes *what* your genetic circuit is — its parts and
interactions. It does not describe *how fast* things happen or how cells
behave physically. The parameters below fill that gap.

---

## 1. Simulation-Level Parameters

These control the overall behaviour of the simulation environment.

| Parameter | Type | Description | Typical value |
|---|---|---|---|
| `max_cells` | int | Maximum number of cells before simulation stops | 10000 |
| `jitter_z` | bool | `false` = 2D colony, `true` = 3D | false |
| `gamma` | float | Frictional drag on cell growth. Higher = cells push harder | 100.0 |
| `pickle_steps` | int | Save state every N steps (lower = more frequent) | 50 |
| `random_seed` | int | Seed for reproducibility | 12345 |

---

## 2. Cell Type Parameters

One block per cell type / strain in your design. If your SBOL describes a
single genetic construct, you likely have one cell type. Multi-strain
designs need one block each.

| Parameter | Type | Description | Typical value |
|---|---|---|---|
| `display_name` | string | Human-readable label | "Strain A" |
| `color` | [R, G, B] | Visualisation color, values 0.0–1.0 | [1.0, 0.3, 0.3] |
| `growth_rate` | float | How fast the cell grows each timestep | 1.0 |
| `division_length` | float | Cell length that triggers division | 3.0 |
| `division_noise` | float | Random variation added to division length | 0.5 |
| `initial_pos` | [x, y, z] | Starting position of the seed cell | [0.0, 0.0, 0.0] |
| `initial_dir` | [x, y, z] | Starting orientation of the seed cell | [1.0, 0.0, 0.0] |

---

## 3. Biochemical Kinetics Parameters

These are the rates that govern gene expression inside each cell. SBOL
describes *what* regulates *what*, but not *how fast*. You must supply these.

| Parameter | Type | Description | Typical value |
|---|---|---|---|
| `production_rate` | float | Basal rate of protein production per timestep | 0.05 |
| `max_production_rate` | float | Maximum rate when fully activated | 0.2 |
| `degradation_rate` | float | Fraction of protein degraded per timestep | 0.01 |
| `hill_coefficient` | float | Cooperativity of activation/repression (n in Hill eq.) | 2.0 |
| `activation_threshold` | float | Signal concentration for half-max activation | 0.5 |
| `repression_threshold` | float | Protein level for half-max repression | 0.5 |

These values apply per interaction found in the SBOL. If you have multiple
interactions with different rates, you can override them individually in the
`interactions` block of the parameters file.

---

## 4. Chemical Signaling Parameters

Only needed if your SBOL design includes diffusible signals between cells.

| Parameter | Type | Description | Typical value |
|---|---|---|---|
| `enabled` | bool | Whether to use the Chemics signaling module | false |
| `grid_size` | int | Resolution of the diffusion grid | 100 |
| `diffusion_rate` | float | How fast each signal spreads (per signal) | 0.1 |
| `signal_degradation_rate` | float | How fast each signal breaks down | 0.01 |
| `boundary_condition` | string | `"periodic"` or `"fixed"` | "fixed" |

---

## 5. SBOL Mapping Overrides (optional)

By default, the converter will try to infer mappings automatically.
Use this section to correct or override those inferences.

| Parameter | Type | Description |
|---|---|---|
| `component_to_celltype` | dict | Map SBOL component display IDs to cell type indices |
| `signal_component_ids` | list | List of SBOL component IDs that are diffusible signals |
| `ignore_component_ids` | list | SBOL components to skip entirely (e.g. backbone, chassis) |

---

## What the SBOL/JSON provides automatically

You do NOT need to supply these — the converter reads them from the SBOL:

- Names and roles of genetic components (promoter, CDS, terminator, etc.)
- Interaction types (inhibition, stimulation, degradation)
- Which component participates in which interaction (activator, inhibited, etc.)
- Number of distinct genetic constructs (informs number of cell types)
- Whether signals are intracellular proteins or extracellular molecules
