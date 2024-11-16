import sys
from psfgen import PsfGen
from HyresBuilder import utils


'''
usage: python psfgen_RNA.py input_pdb_file output_psf_file
note: 1. change the directory of top_RNA.inp
      2. change the segid, which is the name of RNA
'''

pdb = sys.argv[1]
psf = sys.argv[2]
top_inp, param_inp = utils.load_ff('RNA')

gen = PsfGen()
gen.read_topology(top_inp)
gen.add_segment(segid='iCon', pdbfile=pdb, auto_angles=False, auto_dihedrals=False)
gen.write_psf(filename=psf)

print('Finished!')
