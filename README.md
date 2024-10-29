# CG RNA model tutorial 

This is a tutorial for running CG RNA simulations. Below is a short description of all folders in this tutorial:   

**force_fields**: CG RNA force field files (written in the CHARMM format)   

**scripts:** python scripts for building CG model and running simulations   

**examples**: examples for 2KOC RNA   



## Build environment  
well-tested in: OpenMM 8.1+, python 3.8/3.9
1. build a new conda environment   
conda create -n cgrna python=3.9   
2. install openmm 8.1+   
reference: [OpenMM Installation](http://docs.openmm.org/latest/userguide/application/01_getting_started.html#installing-openmm). Use cudatoolkit 11.8 with openmm 8.1.2.      
3. dependencies:  
a. [Mdanalysis](https://www.mdanalysis.org/)   
b. [psfgen-python](https://psfgen.robinbetz.com/) (install through "conda install -c conda-forge psfgen", python < 3.10)   
c. [HyresBuilder](https://github.com/lslumass/HyresBuilder)   


## 1. Prepare PDB file for CG RNA

### A. disordered RNA
For disordered RNA or starting simulation from sequence, please using [build_CG_RNA.py](./scripts/build_CG_RNA.py) to generate a fully linear chain.   
**Usage:** ```python RNA_builder.py output_pdb_file sequence```   
**Example:** ```python RNA_builder.py 2koc_cg.pdb GGCACUUCGGUGCC```   
[2koc_cg.pdb](./examples/disorder/2koc_cg.pdb) will be generated.   

### B. structured RNA
For simulation of structured RNA, one should first download the atomistic PDB file from [RCSB](https://www.rcsb.org/), and then convert it to CG pdb file using [at2cg_RNA.py](./scripts/at2cg_RNA.py).  
**Usage:** ```python at2cg.py input_AA_pdb_file output_CG_pdb_file```   
**Example:** ```python at2cg.py 2KOC.pdb 2koc_cg.pdb```   
[2KOC.pdb](./examples/from-structured/2KOC.pdb) was CHARMM-style PDB file obtained from [CHARMM-GUI](https://www.charmm-gui.org/), and then was converted to [2koc_cg.pdb](./examples/from-structured/2koc_cg.pdb)   
**Note:** up to now, only atomistic pdb file from CHARMM-GUI can be converted correctly. No matter where you get the pdb file from, just convert it to CHARMM-style using CHARMM-GUI (**important:** choose "5PHO" and "3TER" for terminal groups).    

after getting the pdb file, one can create the simulation box using [Packmol](https://m3g.github.io/packmol/) or other tools. Here, we use Packmol to pack the single chain RNA in a box sized 10 nm with [packmol.inp](./examples/disorder/packmol.inp) as an example.   


## 2. Prepare CG PSF file

[psfgen_RNA.py](./scripts/psfgen_RNA.py) was used to create psf file for CG RNA model.   
**Usage:** ```python psfgen_RNA.py input_pdb_file output_psf.file```   
**Example:** ```python psfgen_RNA.py conf.pdb conf.psf```   
[conf.psf](./examples/disorder/conf.psf) was created based on [conf.pdb](./examples/disorder/conf.pdb).   
**Note:** [top_RNA.inp](./force_fields/top_RNA.inp) is necessary for psfgen_RNA.py   


## 3. Run simulation    
After getting pdb and psf files, the simulation was performed using OpenMM through [run.py](./scripts/run.py).   
The correct paths of [param_RNA.inp](./force_fields/param_RNA.inp) and [top_RNA.inp](./force_fields/top_RNA.inp) should be given in run.py.   
**Usage:** ```python run.py -c pdb_file -p psf_file -t temperature -b box -s slat_concentration -m Mg2+ concentraton```   
where salt concentration (monvalent ions only) is in mM, default value is 150 mM. The box is a three float number, e.g., 50 50 50, temperature is in K, default value is 303 K. Mg2+ concentration is in mM, default value is 0.0 mM.    
**Example:** ```python run.py -c conf.pdb -p conf.psf -t 303 -b 50 50 50 -s 150.0 -m 10.0```   
**Note:** default value of GPU_id is 0. Change it on demand.   
