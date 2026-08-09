# CellModeller simulation script
# Generated : 2026-07-13 15:44:03
# Converter : json_to_cellmodeller.py
# Topology (who produces what signal, what it activates) is auto-detected from
# the SBOL JSON.  All numerical rates are defined in the params file / dict.
# NOTE: signalling is ON — protein/signal kinetics live in specRateCL(), not update().


from CellModeller.Regulation.ModuleRegulator import ModuleRegulator
from CellModeller.Biophysics.BacterialModels.CLBacterium import CLBacterium
from CellModeller.Signalling.GridDiffusion import GridDiffusion
from CellModeller.Integration.CLCrankNicIntegrator import CLCrankNicIntegrator
import numpy as np
import random

# EXTERNAL CHEMICAL CONCENTRATIONS  (non-diffusible inducers -- set before running)
PRECURSOR_AHL_CONC = 5.0  # Precursor AHL (Simple chemical) -- substrate -- gates production of Complex LuxR-AHL in specRateCL; 0.0 keeps that pathway off no matter how many cells exist
AHL_CONC = 0.0  # AHL (Simple chemical) -- produced in-model from another species -- this constant has no effect unless you add exogenous-supplementation kinetics yourself

# AUTO-DETECTED SIGNAL TOPOLOGY
#   Diffusible signals      : ['Complex LuxR-AHL']
#   Complex LuxR-AHL       produced by : ['LuxR', 'LuxI']
#   BBa_R0062              activated by: Complex LuxR-AHL

# simulation constants
maxCells = 10000
gridLen     = 400       # grid cells along x and y
gridZCells  = 3   # grid cells along z (small = thin/2D)
gridSize    = 1.0      # micrometers per grid cell (must be isotropic)
gridOrigX, gridOrigY, gridOrigZ = -200.0, -200.0, -1.5


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
    0: {'luxr': 0.0, 'gfp': 0.0, 'luxi': 0.0},  # Strain 0 (GFP)
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

    nSignals = 1  # ['Complex LuxR-AHL']
    nSpecies = 4  # tracked per-cell species (proteins + signal pools)
    sig   = GridDiffusion(sim, nSignals,
                (gridLen, gridLen, gridZCells),
                (gridSize, gridSize, gridSize),   # must be isotropic
                (gridOrigX, gridOrigY, gridOrigZ),
                [0.01])   # diffusion coefficients
    integ = CLCrankNicIntegrator(sim, nSignals, nSpecies, maxCells,
                sig, boundcond='reflect')  # reflect boundary

    sim.init(biophys, regul, sig, integ)
    sim.pickleSteps = 10

    sim.addCell(cellType=0, pos=(0.0, 0.0, 0.0), dir=(1.0, 0.0, 0.0))


# INIT

def init(cell):
    cell.targetVol  = (cell_lens[cell.cellType]
                       + random.uniform(0.0, cell_division_noise[cell.cellType]))
    cell.growthRate = cell_growth_rates[cell.cellType]
    cell.color      = cell_colors[cell.cellType]
    # proteins + signal pools, all tracked via cell.species[]
    cell.species[:] = [
        cell_initial_concentrations[cell.cellType]['luxr'],  # LuxR -> species[0]
        cell_initial_concentrations[cell.cellType]['gfp'],  # GFP -> species[1]
        cell_initial_concentrations[cell.cellType]['luxi'],  # LuxI -> species[2]
        0.0,  # Complex LuxR-AHL -> species[3]
    ]
    # locally-sensed extracellular signal levels
    cell.signals[:] = [0.0] * 1


# UPDATE

def update(cells):
    for id, cell in cells.items():
        if cell.volume > cell.targetVol:
            cell.divideFlag = True

        # colour by GFP expression
        _fp_norm = min(1.0, cell.species[1] / max(2.0, 1e-9))
        cell.color = [
            0.0 * _fp_norm + 0.1 * (1.0 - _fp_norm),
            1.0 * _fp_norm + 0.1 * (1.0 - _fp_norm),
            0.2 * _fp_norm + 0.1 * (1.0 - _fp_norm),
        ]


# DIVIDE

def divide(parent, d1, d2):
    d1.cellType = parent.cellType
    d2.cellType = parent.cellType
    d1.color    = parent.color
    d2.color    = parent.color
    d1.species[:] = parent.species[:] / 2.0
    d2.species[:] = parent.species[:] / 2.0


# SPEC / SIGNAL RATES (OpenCL C — required by CLCrankNicIntegrator)
# specRateCL() sets rates[] for every tracked species (proteins + signal pools).
# sigRateCL()  sets rates[] for the flux of each signal exchanged with the grid.
# Available in both: gridVolume, area, volume, cellType, species[], signals[]
# NOTE: specRateCL() returns an f-string. Any {CONSTANT_NAME} placeholder in
# the body below (e.g. a substrate's *_CONC gate) is evaluated against this
# script's own module-level globals every time specRateCL() is called -- edit
# the constants above before running, not this function.

def specRateCL():
    return f'''
    // — protein production and degradation —
    // LuxR: constitutive from BBa_J23100
    rates[0] = 2.0f - 0.1f * species[0];

    // GFP: activated by diffusible signal Complex LuxR-AHL via BBa_R0062
    float _act_bba_r0062 = pow(signals[0], 4.0f) / (pow(2.0f, 4.0f) + pow(signals[0], 4.0f));
    rates[1] = 2.0f * _act_bba_r0062 - 0.1f * species[1];

    // LuxI: constitutive from BBa_J23119
    rates[2] = 2.0f - 0.1f * species[2];

    // — diffusible-signal pools (produced here, exchanged with grid in sigRateCL) —
    // Complex LuxR-AHL: produced by LuxR, LuxI
    // Complex LuxR-AHL: gated by external substrate(s) Precursor AHL — set the matching _CONC constant(s) above to activate this pathway
    rates[3] = (0.1f * (species[0] + species[2])) * {PRECURSOR_AHL_CONC}f - 0.1f * species[3] - 0.1f * (species[3] - signals[0]) * area / volume;

    '''

def sigRateCL():
    return '''
    // Complex LuxR-AHL: secretion into the grid (species[3] -> signals[0])
    rates[0] = 0.1f * (species[3] - signals[0]) * area / gridVolume;
    '''
