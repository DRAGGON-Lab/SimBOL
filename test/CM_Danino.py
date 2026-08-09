"""
CellModeller model of the synchronized quorum-sensing genetic clock from:
  Danino T, Mondragon-Palomino O, Tsimring L, Hasty J.
  "A synchronized quorum of genetic clocks." Nature 463, 326-330 (2010).
"""

import random
import numpy
from CellModeller.Regulation.ModuleRegulator import ModuleRegulator
from CellModeller.Biophysics.BacterialModels.CLBacterium import CLBacterium
from CellModeller.Signalling.GridDiffusion import GridDiffusion
from CellModeller.Integration.CLCrankNicIntegrator import CLCrankNicIntegrator
from CellModeller.GUI import Renderers

max_cells = 5000

# Module-level definitions
numSpecies = 3   # LuxI, AiiA, GFP
numSignals = 2   # AHL, Nutrient

# Trap + perpendicular flow-channel geometry
# The trap is a straight corridor along x, closed at TRAP_BACK_X and open at
# TRAP_OPEN_X. The flow channel runs perpendicular to it (along y), centered
# on the trap mouth (x = TRAP_OPEN_X) and extending CHANNEL_HALF_LENGTH
# along y on both sides of the opening.
TRAP_HALF_X = 60.0    # trap extends from the opening to the closed back wall
TRAP_HALF_Y = 15.0    # trap corridor half-width
TRAP_HALF_Z = 3.0      # corridor half-height (thin, ~monolayer)
TRAP_OPEN_X = -TRAP_HALF_X   # trap/channel boundary (the trap mouth)
TRAP_BACK_X = TRAP_HALF_X    # closed wall

CHANNEL_WIDTH = 40.0                           # channel extent in x, away from the trap mouth
CHANNEL_FAR_X = TRAP_OPEN_X - CHANNEL_WIDTH    # the channel's solid far wall
CHANNEL_HALF_LENGTH = 120.0                    # channel extends this far along y

#REMOVAL_Y = -CHANNEL_HALF_LENGTH + 20.0        # cells are "washed out" (deleted) once carried this far downstream,
                                                # comfortably inside the grid and short of the channel's own far end

CELL_R = 0.5   # rod/sphere radius used throughout setup() and update()

#FLOW_SPEED = 2.0   # units per unit time, applied along -y (downstream) to cells

#FLOW_BUFFER = 5.0
#FLOW_RAMP = 10.0
#FLOW_MAX_STEP = 0.15 * CELL_R

# AHL/Nutrient diffusion grid
# Covers the trap AND the full width/length of the channel
grid_dim = (64, 72, 4)
grid_size = (4, 4, 4)
grid_orig = (-140, -144, -8)

D_AHL = 40.0
D_NUTRIENT = 20.0

# Nutrient parameters
NUTRIENT_TARGET = 10.0        # steady-state level maintained inside the trap
NUTRIENT_SUPPLY_RATE = 2.0    # how fast it's replenished toward that target
NUTRIENT_DECAY_RATE = 0.5     # how fast it depletes once outside the trap
BASE_GROWTH_RATE = 1.0        # growth rate achieved when nutrient is abundant
NUTRIENT_K = 5.0              # half-saturation constant for growth-rate coupling

# If a cell is still this many times past its own targetVol without having
# managed to divide, treat it as mechanically jammed
JAM_VOL_MULT = 1.5


class GridDiffusionWithFlowAndNutrient(GridDiffusion):
    """
    signals[0] = AHL    
    signals[1] = Nutrient
    """
    def __init__(self, *args, trap_x=TRAP_OPEN_X,
                 ahl_sink_rate=5.0,
                 nutrient_target=NUTRIENT_TARGET,
                 nutrient_supply_rate=NUTRIENT_SUPPLY_RATE,
                 nutrient_decay_rate=NUTRIENT_DECAY_RATE,
                 **kwargs):
        super().__init__(*args, adv=None, **kwargs)
        self.trap_x = trap_x
        self.ahl_sink_rate = ahl_sink_rate
        self.nutrient_target = nutrient_target
        self.nutrient_supply_rate = nutrient_supply_rate
        self.nutrient_decay_rate = nutrient_decay_rate

        nx = self.gridDim[1]
        xs = self.gridOrig[0] + numpy.arange(nx) * self.gridSize[0]
        self.outside_mask = (xs < self.trap_x)    # True = outside trap (in the channel; AHL sink zone)
        self.inside_mask = ~self.outside_mask      # True = inside trap (nutrient supply zone)

    def transportRates(self, signalRates, signalLevels, boundcond='constant', mode='normal'):
        super().transportRates(signalRates, signalLevels, boundcond, mode)

        # CLCrankNicIntegrator probes this method once at startup with
        # mode='greens' to build its implicit solver's operator structure.
        if mode == 'greens':
            return

        # self.gridDim already includes the leading signal-count dimension
        signalRatesView = signalRates.reshape(self.gridDim)
        signalLevelsView = signalLevels.reshape(self.gridDim)

        # AHL: extra loss outside the trap (washed away by the flow channel)
        signalRatesView[0][self.outside_mask, :, :] -= (
            self.ahl_sink_rate * signalLevelsView[0][self.outside_mask, :, :]
        )

        # Nutrient: resupplied inside the trap, decays outside
        signalRatesView[1][self.inside_mask, :, :] += (
            self.nutrient_supply_rate *
            (self.nutrient_target - signalLevelsView[1][self.inside_mask, :, :])
        )
        signalRatesView[1][self.outside_mask, :, :] -= (
            self.nutrient_decay_rate * signalLevelsView[1][self.outside_mask, :, :]
        )


def addWallOfSpheres(biophys, p1, p2, radius, coeff=1.0, spacing_factor=0.5):
    """Place overlapping spheres from p1 to p2 to form a finite wall segment.
    norm=1 -> obstacle: cells are pushed away from the spheres."""
    p1 = numpy.array(p1, dtype=float)
    p2 = numpy.array(p2, dtype=float)
    length = numpy.linalg.norm(p2 - p1)
    spacing = radius * 2 * spacing_factor   # overlap so there are no gaps
    n = max(2, int(numpy.ceil(length / spacing)) + 1)
    for t in numpy.linspace(0, 1, n):
        pt = tuple(p1 + t * (p2 - p1))
        biophys.addSphere(pt, radius, coeff, 1)   # norm=1 -> solid obstacle


def setup(sim):
    # Stash a module-level reference so update() can still reach sim.phys,
    # sim.dt, and sim.delCell.
    global _sim
    _sim = sim

    r = CELL_R
    planeWeight = 1.0


    biophys = CLBacterium(sim, max_cells=max_cells, max_contacts=32,
                           jitter_z=False, max_sqs=64 * 64,
                           max_planes=3, max_spheres=1500,
                           max_substeps=16)

    # Trap side walls: these run the length of the trap itself, from the
    # closed back wall (TRAP_BACK_X) to just short of the mouth.
    WALL_SETBACK = 3.0
    addWallOfSpheres(biophys, (TRAP_OPEN_X + WALL_SETBACK, -TRAP_HALF_Y - r, 0.0), (TRAP_BACK_X, -TRAP_HALF_Y - r, 0.0), r)  # bottom
    addWallOfSpheres(biophys, (TRAP_OPEN_X + WALL_SETBACK, TRAP_HALF_Y + r, 0.0), (TRAP_BACK_X, TRAP_HALF_Y + r, 0.0), r)    # top
    addWallOfSpheres(biophys, (TRAP_BACK_X + r, -TRAP_HALF_Y, 0.0), (TRAP_BACK_X + r, TRAP_HALF_Y, 0.0), r)   # closed back wall

    # Flow channel: perpendicular to the trap, centered on the mouth.
    biophys.addPlane((CHANNEL_FAR_X, 0, 0), (1, 0, 0), planeWeight)

    addWallOfSpheres(biophys, (TRAP_OPEN_X, -CHANNEL_HALF_LENGTH + 3.0, 0.0), (TRAP_OPEN_X, -TRAP_HALF_Y, 0.0), r)
    addWallOfSpheres(biophys, (TRAP_OPEN_X, TRAP_HALF_Y, 0.0), (TRAP_OPEN_X, CHANNEL_HALF_LENGTH - 3.0, 0.0), r)

    # Bridge the corner gap left by WALL_SETBACK, at both mouth corners.
    addWallOfSpheres(biophys, (TRAP_OPEN_X + WALL_SETBACK, -TRAP_HALF_Y - r, 0.0), (TRAP_OPEN_X, -TRAP_HALF_Y - r, 0.0), r)
    addWallOfSpheres(biophys, (TRAP_OPEN_X, -TRAP_HALF_Y - r, 0.0), (TRAP_OPEN_X, -TRAP_HALF_Y, 0.0), r)
    addWallOfSpheres(biophys, (TRAP_OPEN_X + WALL_SETBACK, TRAP_HALF_Y + r, 0.0), (TRAP_OPEN_X, TRAP_HALF_Y + r, 0.0), r)
    addWallOfSpheres(biophys, (TRAP_OPEN_X, TRAP_HALF_Y + r, 0.0), (TRAP_OPEN_X, TRAP_HALF_Y, 0.0), r)

    # z-confinement (applies to trap + channel alike)
    biophys.addPlane((0, 0, TRAP_HALF_Z), (0, 0, -1), planeWeight)    # closed ceiling
    biophys.addPlane((0, 0, -TRAP_HALF_Z), (0, 0, 1), planeWeight)    # closed floor

    # Signalling: AHL (index 0) and Nutrient (index 1)
    sig = GridDiffusionWithFlowAndNutrient(sim, 2, grid_dim, grid_size, grid_orig,
                                            D=[D_AHL, D_NUTRIENT])

    # Regulation: this file also defines the GRN (numSpecies/specRateCL etc.)
    regul = ModuleRegulator(sim, sim.moduleName)

    # Integrator couples intracellular species (LuxI, AiiA, GFP) to the
    # extracellular AHL + Nutrient fields. nSignals bumped to 2 to match sig.
    integ = CLCrankNicIntegrator(sim, 2, 3, max_cells, sig, boundcond='reflect')

    sim.init(biophys, regul, sig, integ)

    # Single founder cell, started near the closed (+x) end of the trap.
    sim.addCell(cellType=0, pos=(TRAP_BACK_X - 5.0, 0, 0))

    if sim.is_gui:
        therenderer = Renderers.GLBacteriumRenderer(sim)
        sim.addRenderer(therenderer)

    sim.pickleSteps = 20


def init(cell):
    cell.targetVol = 3.5 + random.uniform(-0.3, 0.3)
    cell.growthRate = BASE_GROWTH_RATE   # overwritten every step in update() based on local nutrient
    cell.species[:] = [random.uniform(0, 0.2),
                        random.uniform(0, 0.2),
                        0.0]
    cell.color = (0.0, 0.0, 0.0)


def update(cells):
    to_remove = []
    for id, cell in cells.items():
        gfp = cell.species[2]
        if not (gfp == gfp) or gfp < 0:  # NaN check (NaN != NaN) + negative guard
            gfp = 0.0
        # Smooth saturating normalization (Michaelis-Menten style), so you
        # see the burst build up rather than a hard-clamped on/off flash.
        GFP_COLOR_SCALE = 1.0
        g = gfp / (gfp + GFP_COLOR_SCALE)
        cell.color = (0.05, g, 0.05)

        # Growth rate follows the local nutrient level (Monod/Michaelis-Menten):
        # near-max inside the trap where nutrient is resupplied, falling off
        # toward zero once a cell has drifted out into the flow channel.
        nutrient = max(0.0, cell.signals[1])
        cell.growthRate = BASE_GROWTH_RATE * nutrient / (NUTRIENT_K + nutrient)

        if cell.volume > cell.targetVol:
            cell.divideFlag = True

        # Safety valve if a cell is mechanically jammed.
        if cell.volume > JAM_VOL_MULT * cell.targetVol:
            cell.growthRate = 0.0

        # Flow: once a cell has exited the trap into the channel.
        #if cell.pos[0] < TRAP_OPEN_X - FLOW_BUFFER:
            #if cell.pos[1] > REMOVAL_Y:
                #depth = (TRAP_OPEN_X - FLOW_BUFFER) - cell.pos[0]
                #ramp = min(1.0, depth / FLOW_RAMP)
                #step = FLOW_SPEED * ramp * _sim.dt
                #step = min(step, FLOW_MAX_STEP)
                #delta = numpy.array([0.0, -1.0, 0.0]) * step
                #_sim.phys.moveCell(cell, tuple(delta))
            #else:
                #to_remove.append(id)
        #elif cell.pos[0] < TRAP_OPEN_X:
            #pass  # in the decompression buffer just past the mouth -- no flow yet

    #for id in to_remove:
        #_sim.delCell(id)


def divide(parent, d1, d2):
    d1.targetVol = 3.5 + random.uniform(-0.3, 0.3)
    d2.targetVol = 3.5 + random.uniform(-0.3, 0.3)
    # species carried over as concentrations (roughly conserved on division,
    # small noise added so daughters aren't perfectly identical)
    for d in (d1, d2):
        d.species[:] = [max(0.0, s * random.uniform(0.9, 1.1))
                         for s in parent.species]


# Reaction kernels (OpenCL C, inlined by CellModeller)
# These return strings of C code, not Python -- species[], signals[], and
# rates[] are float arrays inside the compiled kernel. signals[] now has 2
# entries: signals[0] = AHL, signals[1] = Nutrient.

def sigRateCL():
    return '''
    const float k_synth = 8.0f;   // AHL production per unit LuxI
    const float k_degr  = 4.0f;   // AHL degradation per unit AiiA

    // clamp to non-negative: concentrations can't physically go negative,
    // and letting them dip below 0 numerically here can otherwise poison
    // the Hill calc in specRateCL (see note there) with NaN forever
    float LuxI = fmax(species[0], 0.0f);
    float AiiA = fmax(species[1], 0.0f);
    float A    = fmax(signals[0], 0.0f);

    // net AHL exchanged with the shared grid by this cell
    rates[0] = k_synth*LuxI - k_degr*AiiA*A;

    // Nutrient is not directly produced/consumed per-cell here
    // exchange term is added.
    rates[1] = 0.0f;
    '''


def specRateCL():
    return '''
    // basal + AHL-activated transcription from the (shared) lux promoter,
    // independent first-order decay/dilution of each protein
    const float alpha0  = 0.02f;   // basal expression
    const float alpha1  = 8.0f;    // max activated expression
    const float K       = 2.0f;    // AHL EC50
    const float n        = 3.0f;   // Hill coefficient (cooperative LuxR-AHL binding)

    const float gamma_LuxI = 1.2f;  // LuxI turnover
    const float gamma_AiiA = 0.3f;  // AiiA turnover (slower -> relaxation osc.)
    const float gamma_GFP  = 0.5f;  // GFP is the readout, degrades slowest

    // IMPORTANT: powr(x,y) returns NaN for any x < 0 (per OpenCL spec, unlike
    // pow). The integrator can occasionally overshoot AHL slightly below 0
    // numerically, and that single negative value will poison every species
    // with NaN from then on. Clamp before the Hill calc to prevent this.
    float A  = fmax(signals[0], 0.0f);
    float An = powr(A, n);
    float hill = An / (powr(K, n) + An);

    float LuxI = fmax(species[0], 0.0f);
    float AiiA = fmax(species[1], 0.0f);
    float GFP  = fmax(species[2], 0.0f);

    rates[0] = alpha0 + alpha1*hill - gamma_LuxI*LuxI; // LuxI
    rates[1] = alpha0 + alpha1*hill - gamma_AiiA*AiiA; // AiiA
    rates[2] = alpha0 + alpha1*hill - gamma_GFP*GFP;   // GFP (reporter)
    '''
