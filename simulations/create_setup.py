#!/usr/bin/env python3

hole_radius = 4 # nm 
rates = [0.5,5,10,15,30,50] #bar/ns
rates = [0.5] #bar/ns
polarizations = ["00", "04", "06", "07", "08", "10"]

from os import mkdir, chdir, path
from shutil import copy
from subprocess import run

import numpy as np
from MDAnalysis import Universe

def freplace(fold, fnew, old, new):
    """All occurrences of substring `old` in file `fold` are replaced by `new` and
    written to `fnew`. Old and new can also be lists of the same lengths."""

    if type(old) not in [list, tuple, np.ndarray]:
        old = [old]

    if type(new) not in [list, tuple, np.ndarray]:
        new = [new]

    if len(old) != len(new):
        raise ValueError("Old and new are not of the same length.")

    with open(fold,"r") as f:
        string = f.read()

    for i in range(len(old)):
        string = string.replace(str(old[i]), str(new[i]))

    with open(fnew,"w") as f:
        f.write(string)

def mkdir_force(path):
    """Create a directory named path without raising an error if the directory
    already exists."""
    try:
        mkdir(path)
    except OSError:
        pass
    

script_path = path.dirname(path.realpath(__file__))
dirname = f"r_{hole_radius}"
mkdir_force(dirname)

steps = 1/np.array(rates) # ns/bar
steps *= 1000 # ns to ps
steps /= 0.002 #stepsize in ps

# Select all SAMs not within given hole_readius of residue 220
for gro in ["restraint.gro", "conf.gro"]:
    u = Universe(f"ref_coordinates/{gro}")
    XX1 = u.select_atoms("resname XX1")
    XX2 = u.select_atoms("resname XX2")
    SOL = u.select_atoms("resname SOL")
    
    ref_atoms = XX1.select_atoms("name OAA")
    XX1 = u.select_atoms("resname XX1")
    # Always take residues from restrained gro file
    if gro == "restraint.gro":
        res_resids = ref_atoms.select_atoms(f"sphlayer {10*hole_radius} {np.inf} resid 101").resindices
    XX1_hole = XX1.residues[res_resids].atoms
    new_u = XX1_hole + XX2 + SOL
    new_u.write(f"{dirname}/{gro}")

for pol in polarizations:
    pol_dirname = f"{dirname}/{pol}_10"
    mkdir_force(pol_dirname)

    for itp in [f"C10OH-XX1_{pol}", "C10OH-XX2_10", "posreCH3"]:
        copy(f"topology/{itp}.itp", f"{pol_dirname}/{itp}.itp")

    freplace("topology/topol.top",
            f"{pol_dirname}/topol.top",
            ["aa", "n_XX1_residues"],
            [pol, XX1_hole.n_residues])
    
    for rate, n_steps in zip(rates, steps):
        rate_dirname = f"{pol_dirname}/{rate}_bar_ns"
        mkdir_force(rate_dirname)
        freplace("grompp.mdp",
                f"{rate_dirname}/grompp.mdp",
                "//STEPS//",
                int(n_steps))
        chdir(rate_dirname)
        run(f"{script_path}/create_realizations.sh")
        chdir(script_path)
