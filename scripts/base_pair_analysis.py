#!/usr/bin/python

"""
this script is used for the analysis of general hydrogen bond, meaning donor-acceptor system, based on MDAnalysis
authour: Shanlong Li@UMass
date: 09-04-2024
"""

from __future__ import division
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis

u = mda.Universe('conf.psf', 'system.dcd')
# calculate GC pairs
GCpairs = HydrogenBondAnalysis(
    universe=u,
    donors_sel="resname GUA and name NB",
    hydrogens_sel="resname GUA and name NC",
    acceptors_sel="resname CYT and name NB",
    d_a_cutoff=10.0,
    d_h_cutoff=3.5,
    d_h_a_angle_cutoff=150, 
    update_selections=False
)

GCpairs.run(step=1, verbose=True)

counts = GCpairs.count_by_time()
with open('pair_number_GC.dat', 'w') as f:
    for frame, num in zip(GCpairs.frames, counts):
        print((frame+1)*u.trajectory.dt, num, file=f)

# calculate AU pairs
AUpairs = HydrogenBondAnalysis(
    universe=u,
    donors_sel="resname ADE and name NB",
    hydrogens_sel="resname ADE and name NC",
    acceptors_sel="resname URA and name NB",
    d_a_cutoff=10.0,
    d_h_cutoff=3.5,
    d_h_a_angle_cutoff=150, 
    update_selections=False
)

AUpairs.run(step=1, verbose=True)

counts = AUpairs.count_by_time()
with open('pair_number_AU.dat', 'w') as f:
    for frame, num in zip(AUpairs.frames, counts):
        print((frame+1)*u.trajectory.dt, num, file=f)

print('Finished!')

 
