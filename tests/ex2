# CellModeller simulation script
# Generated : 2026-06-30 17:10:13
# Converter : json_to_cellmodeller.py
# Topology (who produces what signal, what it activates) is auto-detected from
# the SBOL JSON.  All numerical rates are defined in the params file / dict.


from CellModeller.Regulation.ModuleRegulator import ModuleRegulator
from CellModeller.Biophysics.BacterialModels.CLBacterium import CLBacterium
# from CellModeller.Signalling.GridDiffusion import GridDiffusion
# from CellModeller.Integration.CLCrankNicIntegrator import CLCrankNicIntegrator
import numpy as np
import random


# simulation constants
maxCells = 10000
gridLen  = 100   # grid cells per axis
gridSize = 4.0  # µm per grid cell  →  domain = 400 × 400 µm


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
    cell.ci = 20.0 # cI
    cell.gfp = 0.0  # GFP
    cell.laci = 0.0  # LacI
    cell.tetr = 0.0  # TetR


# UPDATE

def update(cells):
    for id, cell in cells.items():
        if cell.volume > cell.targetVol:
            cell.divideFlag = True

        # colour by GFP expression
        _fp_norm = min(1.0, cell.gfp / 20.0)
        cell.color = [
            0.0 * _fp_norm + 0.1 * (1.0 - _fp_norm),
            1.0 * _fp_norm + 0.1 * (1.0 - _fp_norm),
            0.2 * _fp_norm + 0.1 * (1.0 - _fp_norm),
        ]

        # — inhibition factors —
        _inh_bba_r0011 = 2.0**4.0 / (2.0**4.0 + cell.laci**4.0)  # LacI represses BBa_R0011
        _inh_bba_r0040 = 2.0**4.0 / (2.0**4.0 + cell.tetr**4.0)  # TetR represses BBa_R0040
        _inh_bba_r0050 = 2.0**4.0 / (2.0**4.0 + cell.ci**4.0)  # cI represses BBa_R0050

        # — protein production and degradation —
        # cI: repressed by TetR via BBa_R0040
        cell.ci += 2.0 * _inh_bba_r0040 - 0.1 * cell.ci
        cell.ci = max(0.0, cell.ci)

        # GFP: repressed by LacI via BBa_R0011
        cell.gfp += 2.0 * _inh_bba_r0011 - 0.1 * cell.gfp
        cell.gfp = max(0.0, cell.gfp)

        # LacI: repressed by cI via BBa_R0050
        cell.laci += 2.0 * _inh_bba_r0050 - 0.1 * cell.laci
        cell.laci = max(0.0, cell.laci)

        # TetR: repressed by LacI via BBa_R0011
        cell.tetr += 2.0 * _inh_bba_r0011 - 0.1 * cell.tetr
        cell.tetr = max(0.0, cell.tetr)



# DIVIDE

def divide(parent, d1, d2):
    d1.cellType = parent.cellType
    d2.cellType = parent.cellType
    d1.color    = parent.color
    d2.color    = parent.color
    d1.ci = parent.ci / 2.0
    d2.ci = parent.ci / 2.0
    d1.gfp = parent.gfp / 2.0
    d2.gfp = parent.gfp / 2.0
    d1.laci = parent.laci / 2.0
    d2.laci = parent.laci / 2.0
    d1.tetr = parent.tetr / 2.0
    d2.tetr = parent.tetr / 2.0

