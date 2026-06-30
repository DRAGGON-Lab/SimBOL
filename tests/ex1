# CellModeller simulation script
# Generated : 2026-06-24 18:02:09
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
gridLen  = 100    # grid cells per axis
gridSize = 4.0    # µm per grid cell  →  domain = 400 × 400 µm

# repressilator circuit parameters
ALPHA = 2.0    # max protein production per step
BETA  = 0.1    # protein degradation per step
K     = 2.0    # Hill half-saturation constant (protein units)
N     = 4.0    # Hill coefficient
P_MAX = ALPHA / BETA    # 20 — unrepressed steady-state level

# cell type lookup tables
cell_colors = {
    0: [0.0, 0.1, 0.0],  # dim green at t=0; overwritten by update() each step
}
cell_lens = {
    0: 3.5,
}
cell_growth_rates = {
    0: 0.5,
}
cell_division_noise = {
    0: 0.005,   
}


# SETUP

def setup(sim):
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


    cell.ci   = P_MAX   # cI at peak
    cell.gfp  = 0.0  
    cell.laci = 0.0
    cell.tetr = 0.0


# UPDATE

def update(cells):
    for id, cell in cells.items():
        if cell.volume > cell.targetVol:
            cell.divideFlag = True

        
        _fp_norm = min(1.0, cell.gfp / P_MAX)
        cell.color = [
            0.0 * _fp_norm + 0.1 * (1.0 - _fp_norm),   # R: near-zero
            1.0 * _fp_norm + 0.1 * (1.0 - _fp_norm),   # G: bright when GFP high
            0.2 * _fp_norm + 0.1 * (1.0 - _fp_norm),   # B: slight tint
        ]

        # Hill repression factors
        _inh_bba_r0011 = K**N / (K**N + cell.laci**N)  # LacI represses BBa_R0011
        _inh_bba_r0040 = K**N / (K**N + cell.tetr**N)  # TetR represses BBa_R0040
        _inh_bba_r0050 = K**N / (K**N + cell.ci**N)    # cI  represses BBa_R0050

        # protein production and degradation
        # cI: repressed by TetR via BBa_R0040
        cell.ci   += ALPHA * _inh_bba_r0040 - BETA * cell.ci
        cell.ci    = max(0.0, cell.ci)

        # GFP: repressed by LacI via BBa_R0011
        cell.gfp  += ALPHA * _inh_bba_r0011 - BETA * cell.gfp
        cell.gfp   = max(0.0, cell.gfp)

        # LacI: repressed by cI via BBa_R0050
        cell.laci += ALPHA * _inh_bba_r0050 - BETA * cell.laci
        cell.laci  = max(0.0, cell.laci)

        # TetR: repressed by LacI via BBa_R0011
        cell.tetr += ALPHA * _inh_bba_r0011 - BETA * cell.tetr
        cell.tetr  = max(0.0, cell.tetr)


# DIVIDE

def divide(parent, d1, d2):
    d1.cellType = parent.cellType
    d2.cellType = parent.cellType
    d1.color    = parent.color
    d2.color    = parent.color
    d1.ci   = parent.ci   / 2.0
    d2.ci   = parent.ci   / 2.0
    d1.gfp  = parent.gfp  / 2.0
    d2.gfp  = parent.gfp  / 2.0
    d1.laci = parent.laci / 2.0
    d2.laci = parent.laci / 2.0
    d1.tetr = parent.tetr / 2.0
    d2.tetr = parent.tetr / 2.0
