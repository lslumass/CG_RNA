import sys
from psfgen import PsfGen


'''
usage: python psfgen_RNA.py input_pdb_file output_psf_file
note: 1. change the directory of top_RNA.inp
      2. change the segid, which is the name of RNA
'''

pdb = sys.argv[1]
psf = sys.argv[2]

gen = PsfGen()
gen.read_topology('../force_fields/top_RNA.inp')
gen.add_segment(segid='2KOC', pdbfile=pdb, auto_angles=False, auto_dihedrals=False)
gen.write_psf(filename=psf)

