from __future__ import division, print_function
import argparse
from HyresBuilder import HyresFF
# OpenMM Imports
from openmm.unit import *
from openmm.app import *
from openmm import *
import numpy as np

# 0) set variables in the simulation
gpu_id = "0"
top_inp = './force_fields/top_RNA.inp'
param_inp = './force_fields/param_RNA.inp'

# input parameters
parser = argparse.ArgumentParser()
parser.add_argument('-c', "--pdb", default='conf.pdb', help="pdb file")
parser.add_argument('-p', "--psf", default='conf.psf', help="psf file")
parser.add_argument('-t', "--temp", default=303, type=float, help="system temperature")
parser.add_argument('-b', "--box", nargs='+', type=float, help="box dimensions in nanometer, e.g., '50 50 50' ")
parser.add_argument('-s', "--salt", default=150.0, type=float, help="salt concentration in mM")
parser.add_argument('-m', "--Mg", default=0.0, type=float, help="Mg2+ concentration in mM")

args = parser.parse_args()
pdb_file = args.pdb
psf_file = args.psf
T = args.temp
# pbc box length
if len(args.box) == 1:
    lx, ly, lz = args.box[0], args.box[0], args.box[0]
elif len(args.box) == 3:
    lx = args.box[0]*unit.nanometer
    ly = args.box[1]*unit.nanometer
    lz = args.box[2]*unit.nanometer
else:
    print("Error: You must provide either one or three values for box.")
    exit(1)
c_ion = args.salt/1000.0                                   # concentration of ions in M
c_Mg = args.Mg                                           # concentration of Mg in mM

# simulation parameters
dt = 0.0001*unit.picoseconds		                                # time step
total_step = 250000000                                             # total step
equil_step = 10000
temperture = T*unit.kelvin                                             # temperature
log_freq = 1000                                                 # frequency of log file
dcd_freq = 5000                                                 # frequency of dcd file
pdb_freq = 12500000                                             # frequence of dpd_traj file
pressure = 1*unit.atmosphere                                    # pressure in NPT
friction = 0.1/unit.picosecond                                  # friction coefficient in Langevin

# force field parameters
Td = T-273
er_t = 87.74-0.4008*Td+9.398*10**(-4)*Td**2-1.41*10**(-6)*Td**3
print('relative electric constant: ', er_t*20.3/77.6)                        
dh = 0.304/(np.sqrt(c_ion))
print('Debye-Huckel screening length: ', dh)
nMg = 0.526*(c_Mg/0.680)**(0.283)/(1+(c_Mg/0.680)**(0.283)) + 0.0012*(Td-30) 
lmd0 = nMg*0.64+0.536                                                        
print('lmd: ', lmd0)
ffs = {
    'temp': T,                                                  # Temperature
    'lmd': lmd0,                                                # Charge scaling factor of P-
    'dh': dh*unit.nanometer,                                  # Debye Huckel screening length
    'ke': 138.935456,                                           # Coulomb constant, ONE_4PI_EPS0
    'er': er_t*20.3/77.6,                                         # relative dielectric constant
    'eps_hb': 2.0*unit.kilocalorie_per_mole,                    # hydrogen bond strength
    'sigma_hb': 0.29*unit.nanometer,                            # sigma of hydrogen bond
    'eps_base': 2.05*unit.kilocalorie_per_mole,                 # base stacking strength
}

# 1) import coordinates and topology form charmm pdb and psf
print('\n################## load coordinates, topology and parameters ###################')
pdb = PDBFile(pdb_file)
psf = CharmmPsfFile(psf_file)
psf.setBox(lx, lx, lx)
top = psf.topology
params = CharmmParameterSet(top_inp, param_inp)
system = psf.createSystem(params, nonbondedMethod=CutoffPeriodic, constraints=HBonds)
system = HyresFF.SimpleRNASystem(psf, system, ffs)
system.setDefaultPeriodicBoxVectors(a, b, c)

with open('system.xml', 'w') as output:
    output.write(XmlSerializer.serialize(system))

print('\n# Now, the system has:')
for force in system.getForces():
    print('      ', force.getName())

# simulation
print('\n################### prepare simulation system with NVT ####################')
integrator = LangevinMiddleIntegrator(temperture, friction, dt)
plat = Platform.getPlatformByName('CUDA')
prop = {'Precision': 'mixed', 'DeviceIndex': gpu_id}
simulation = Simulation(top, system, integrator, plat, prop)
simulation.context.setPositions(pdb.positions)
simulation.context.setVelocitiesToTemperature(temperture)

print('\n################### Minimization, Equilibriation, Production simulation ####################')
print('\n# minimizeEnergy:')
print(simulation.context.getState(getEnergy=True).getPotentialEnergy())
simulation.minimizeEnergy(maxIterations=50000, tolerance=0.1)
print(simulation.context.getState(getEnergy=True).getPotentialEnergy())

print('\n# Equilibriation:')
simulation.step(equil_step)
simulation.integrator.setStepSize(0.008*unit.picoseconds)
simulation.step(equil_step)

simulation.reporters.append(PDBReporter('system.pdb', pdb_freq))
simulation.reporters.append(DCDReporter('system.dcd', dcd_freq))
simulation.reporters.append(StateDataReporter('system.log', log_freq, progress=True, remainingTime=True, speed=True, totalSteps=total_step, temperature=True))
#simulation.reporters.append(CheckpointReporter('system.chk', dcd_freq*10))

print('\n# NVT simulation:')
simulation.step(total_step)

simulation.saveCheckpoint('system.chk')
print('\n# Finished!')
