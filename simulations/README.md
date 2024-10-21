# Simulation input

These subolders contain structure, topology and parameter files to run the constant
presure rate simulations. Since GROMACS does not support constant pressure rates the
file `srun_pull_template.sh` is the actual script to perform the pressure ramp
simulations by running sequential simulations while reducing the pressure and reducing
the pressure by 1 bar in each cycle. Different rates are realized by changing the number
os simulations steps for each simulation run inside `grompp.mdp`.

The subfolder `r_x` contain the structure files for different surface polarzbolities.
The `x` denotes the defect size in nm. Values of the files in `r_0`/`r_1`/`r_2` are the
polarizations of the SAM heads which can be mapped to the surface tension according to

    0.0 -> 134.69
    0.4 -> 119.50
    0.6 -> 96.93
    0.7 -> 76.17
    0.76 -> 60.00
    0.8 -> 45.04
    0.87 -> 0.00