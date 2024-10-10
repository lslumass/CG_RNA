import sys
import numpy as np
import MDAnalysis as mda


'''
Usage: python at2cg_RNA.py input_AA_pdb_file output_CG_pdb_file
'''

inp = sys.argv[1]   #input atomistic pdb file
out = sys.argv[2]   #output hyres pdb file
f = open(out, 'w')

# output in pdb-format
def printcg(atom, f):
    f.write("%4s  %5d %2s   %3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f      %4s\n" % (atom[0],int(atom[1]),atom[2],atom[3][:3],atom[4], int(atom[5]),float(atom[6]),float(atom[7]),float(atom[8]),float(atom[9]),float(atom[10]), atom[11]))

u = mda.Universe(inp)
idx = 1
with open(out, 'w') as f:
    for res in u.residues:
        resname = res.resname
        seg = res.segid
        resid = res.resid
        sel = "(name P O1P O2P O5' and resid "+str(resid)+") or (name O3' and resid "+str(resid-1)+")"
        sel = u.select_atoms(sel)
        com = sel.center_of_mass()
        atom = ['ATOM', idx, 'P', res.resname, 'X', resid, com[0], com[1], com[2], 1.00, 0.00, seg]
        printcg(atom, f)
        idx += 1
        sel = u.select_atoms("name C4' and resid "+str(resid))
        com = sel.center_of_mass()
        atom = ['ATOM', idx, 'C1', res.resname, 'X', resid, com[0], com[1], com[2], 1.00, 0.00, seg]
        printcg(atom, f)
        idx += 1
        sel = u.select_atoms("name C1' and resid "+str(resid))
        com = sel.center_of_mass()
        atom = ['ATOM', idx, 'C2', res.resname, 'X', resid, com[0], com[1], com[2], 1.00, 0.00, seg]
        printcg(atom, f)
        idx += 1

        if resname == 'ADE':
            sel = u.select_atoms("name N9 C4 and resid "+str(resid))
            com = sel.center_of_mass()
            atom = ['ATOM', idx, 'NA', resname, 'X', resid, com[0], com[1], com[2], 1.00, 0.00, seg]
            printcg(atom, f)
            idx += 1
            sel = u.select_atoms("name C8 H8 N7 C5 and resid "+str(resid))
            com = sel.center_of_mass()
            atom = ['ATOM', idx, 'NB', resname, 'X', resid, com[0], com[1], com[2], 1.00, 0.00, seg]
            printcg(atom, f)
            idx += 1
            sel = u.select_atoms("name C6 N1 N6 H61 H62 and resid "+str(resid))
            com = sel.center_of_mass()
            atom = ['ATOM', idx, 'NC', resname, 'X', resid, com[0], com[1], com[2], 1.00, 0.00, seg]
            printcg(atom, f)
            idx += 1
            sel = u.select_atoms("name C2 H2 N3 and resid "+str(resid))
            com = sel.center_of_mass()
            atom = ['ATOM', idx, 'ND', resname, 'X', resid, com[0], com[1], com[2], 1.00, 0.00, seg]
            printcg(atom, f)
            idx += 1
        elif resname == 'GUA':
            sel = u.select_atoms("name N9 C4 and resid "+str(resid))
            com = sel.center_of_mass()
            atom = ['ATOM', idx, 'NA', resname, 'X', resid, com[0], com[1], com[2], 1.00, 0.00, seg]
            printcg(atom, f)
            idx += 1
            sel = u.select_atoms("name C8 H8 N7 C5 and resid "+str(resid))
            com = sel.center_of_mass()
            atom = ['ATOM', idx, 'NB', resname, 'X', resid, com[0], com[1], com[2], 1.00, 0.00, seg]
            printcg(atom, f)
            idx += 1
            sel = u.select_atoms("name C6 N1 N6 H1 O6 and resid "+str(resid))
            com = sel.center_of_mass()
            atom = ['ATOM', idx, 'NC', resname, 'X', resid, com[0], com[1], com[2], 1.00, 0.00, seg]
            printcg(atom, f)
            idx += 1
            sel = u.select_atoms("name C2 N2 H21 H22 N3 and resid "+str(resid))
            com = sel.center_of_mass()
            atom = ['ATOM', idx, 'ND', resname, 'X', resid, com[0], com[1], com[2], 1.00, 0.00, seg]
            printcg(atom, f)
            idx += 1
        elif resname == 'CYT':
            sel = u.select_atoms("name N1 C5 H5 C6 H6 and resid "+str(resid))
            com = sel.center_of_mass()
            atom = ['ATOM', idx, 'NA', resname, 'X', resid, com[0], com[1], com[2], 1.00, 0.00, seg]
            printcg(atom, f)
            idx += 1
            sel = u.select_atoms("name C4 N4 H41 H42 N3 and resid "+str(resid))
            com = sel.center_of_mass()
            atom = ['ATOM', idx, 'NB', resname, 'X', resid, com[0], com[1], com[2], 1.00, 0.00, seg]
            printcg(atom, f)
            idx += 1
            sel = u.select_atoms("name C2 O2 and resid "+str(resid))
            com = sel.center_of_mass()
            atom = ['ATOM', idx, 'NC', resname, 'X', resid, com[0], com[1], com[2], 1.00, 0.00, seg]
            printcg(atom, f)
            idx += 1
        elif resname == 'URA':
            sel = u.select_atoms("name N1 C5 H5 C6 H6 and resid "+str(resid))
            com = sel.center_of_mass()
            atom = ['ATOM', idx, 'NA', resname, 'X', resid, com[0], com[1], com[2], 1.00, 0.00, seg]
            printcg(atom, f)
            idx += 1
            sel = u.select_atoms("name C4 O4 H3 H3 and resid "+str(resid))
            com = sel.center_of_mass()
            atom = ['ATOM', idx, 'NB', resname, 'X', resid, com[0], com[1], com[2], 1.00, 0.00, seg]
            printcg(atom, f)
            idx += 1
            sel = u.select_atoms("name C2 O2 and resid "+str(resid))
            com = sel.center_of_mass()
            atom = ['ATOM', idx, 'NC', resname, 'X', resid, com[0], com[1], com[2], 1.00, 0.00, seg]
            printcg(atom, f)
            idx += 1
    print('END', file=f)
quit()
