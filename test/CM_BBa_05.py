# CellModeller simulation script
# Generated : 2026-07-08 16:11:47
# Converter : json_to_cellmodeller.py
# Topology (who produces what signal, what it activates) is auto-detected from
# the SBOL JSON.  All numerical rates are defined in the params file / dict.



from CellModeller.Regulation.ModuleRegulator import ModuleRegulator
from CellModeller.Biophysics.BacterialModels.CLBacterium import CLBacterium
# from CellModeller.Signalling.GridDiffusion import GridDiffusion
# from CellModeller.Integration.CLCrankNicIntegrator import CLCrankNicIntegrator
import numpy as np
import random

# EXTERNAL CHEMICAL CONCENTRATIONS  (non-diffusible inducers -- set before running)
ATC_CONC = 1.0  # aTc (Simple chemical) -- not referenced by any interaction in this circuit -- unused


# simulation constants
maxCells = 10000
gridLen  = 100   # grid cells per axis
gridSize = 4.0  # µm per grid cell


# cell type lookup tables

cell_colors          = {
    0: [0.0, 1.0, 0.2],  # Strain 0 (GFP)
}
cell_lens            = {
    0: 3.5,
}
cell_growth_rates    = {
    0: 1.0,
}
cell_division_noise  = {
    0: 0.005,
}
cell_initial_concentrations = {
    0: {'tetr': 1.0, 'gfp': 0.0, 'ci': 0.0, 'laci': 0.0},  # Strain 0 (GFP)
}


# SETUP

def setup(sim):
    # tip: set simulation.random_seed in params for reproducibility
    biophys = CLBacterium(
        sim,
        jitter_z=False,
        max_cells=maxCells,
        gamma=100.0
    )
    regul = ModuleRegulator(sim, sim.moduleName)

    # Signalling disabled — set signaling.enabled=true in params to activate

    sim.init(biophys, regul, None, None)
    sim.pickleSteps = 10

    sim.addCell(cellType=0, pos=(0.0, 0.0, 0.0), dir=(1.0, 0.0, 0.0))


# INIT

def init(cell):
    cell.targetVol  = (cell_lens[cell.cellType]
                       + random.uniform(0.0, cell_division_noise[cell.cellType]))
    cell.growthRate = cell_growth_rates[cell.cellType]
    cell.color      = cell_colors[cell.cellType]
    # proteins
    cell.tetr = cell_initial_concentrations[cell.cellType]['tetr']  # TetR
    cell.gfp = cell_initial_concentrations[cell.cellType]['gfp']  # GFP
    cell.ci = cell_initial_concentrations[cell.cellType]['ci']  # cI
    cell.laci = cell_initial_concentrations[cell.cellType]['laci']  # LacI


# UPDATE

def update(cells):
    for id, cell in cells.items():
        if cell.volume > cell.targetVol:
            cell.divideFlag = True

        # colour by GFP expression
        _fp_norm = min(1.0, cell.gfp / max(2.0, 1e-9))
        cell.color = [
            0.0 * _fp_norm + 0.1 * (1.0 - _fp_norm),
            1.0 * _fp_norm + 0.1 * (1.0 - _fp_norm),
            0.2 * _fp_norm + 0.1 * (1.0 - _fp_norm),
        ]

        # — inhibition factors —
        _inh_bba_k909012 = (2.0**4.0 / (2.0**4.0 + cell.ci**4.0)) * (2.0**4.0 / (2.0**4.0 + cell.laci**4.0))  # cI, LacI repress BBa_K909012
        _inh_bba_r0040 = 2.0**4.0 / (2.0**4.0 + cell.tetr**4.0)  # TetR represses BBa_R0040

        # — protein production and degradation —
        # TetR: constitutive from BBa_J23100
        cell.tetr += 2.0 - 0.1 * cell.tetr
        cell.tetr = max(0.0, cell.tetr)

        # GFP: repressed by cI, LacI via BBa_K909012
        cell.gfp += 2.0 * _inh_bba_k909012 - 0.1 * cell.gfp
        cell.gfp = max(0.0, cell.gfp)

        # cI: repressed by TetR via BBa_R0040
        cell.ci += 2.0 * _inh_bba_r0040 - 0.1 * cell.ci
        cell.ci = max(0.0, cell.ci)

        # LacI: repressed by cI, LacI via BBa_K909012
        cell.laci += 2.0 * _inh_bba_k909012 - 0.1 * cell.laci
        cell.laci = max(0.0, cell.laci)

        # Direct inhibition: aTc (external) → TetR
        cell.tetr *= 2.0**4.0 / (2.0**4.0 + ATC_CONC**4.0)



# DIVIDE

def divide(parent, d1, d2):
    d1.cellType = parent.cellType
    d2.cellType = parent.cellType
    d1.color    = parent.color
    d2.color    = parent.color
    d1.tetr = parent.tetr / 2.0
    d2.tetr = parent.tetr / 2.0
    d1.gfp = parent.gfp / 2.0
    d2.gfp = parent.gfp / 2.0
    d1.ci = parent.ci / 2.0
    d2.ci = parent.ci / 2.0
    d1.laci = parent.laci / 2.0
    d2.laci = parent.laci / 2.0
