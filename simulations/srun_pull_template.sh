#!/bin/bash

#SBATCH --job-name=
#SBATCH --mail-user=ploche@physik.fu-berlin.de
#SBATCH --mail-type=FAIL

#SBATCH --output=job_%x.out
#SBATCH --ntasks-per-node=16
#SBATCH --mem-per-cpu=1G
#SBATCH --time=4-00:00:00
#SBATCH --nodes=1
#SBATCH --partition=gpu
#SBATCH --gres=gpu:rtx2080ti:4

set -e
module load gromacs/single/cuda/2019.6

unset OMP_NUM_THREADS

Lz0=`tail -1 ../../../conf.gro |awk '{print $3}'`
echo $Lz0

GRO=../../../../conf.gro
# Simulation run
for i in {310..5000}; do
    p=$((0-$i))

    echo $p
    mkdir -p p_$p
    cd p_$p

    # Set pressure
    sed -e "s|//PRESSURE//|$p|g" ../../grompp.mdp > grompp.mdp

    # Run simulations
    gmx --quiet grompp \
	    -c $GRO \
	    -p ../../../topol.top \
	    -r ../../../../restraint.gro \
	    -n ../../../../index.ndx \
	    -maxwarn 3

    gmx --quiet mdrun -cpi state.cpt

    # Analyze
    echo 19 | gmx --quiet energy -f ener.edr 
    Lz=`tail -1 energy.xvg |awk '{print $2}'`  
 
    yes=`echo 1 | awk '{if ('$Lz'>1.5*'$Lz0') print 1;else print 0;}'`
    echo current size = $Lz
    echo yes = $yes

    # if the file doesn't exsist --> crashed due to explosion!
    [ -f confout.gro ] || yes='1'

    # if box became very large --> Lx, Ly, Lz fused in GRO file
    fused=`tail -1 confout.gro |wc -w`
    if [ $fused != '3' ]
    then
	yes='1'
	echo "Sizes fused in GRO"
    fi

    # which serial simulations = which entry in Pmax.dat
    prewords=0
    [ -f ../Pmax.dat ] && prewords=`cat ../Pmax.dat |wc -l`
    words=$(($prewords+1))

    if [ $yes == '1' ]
    then
	echo $p >> ../Pmax.dat
	words=`cat ../Pmax.dat |wc -l`
	echo now quitting!
        cd ..
	break
    fi
 
    GRO=../p\_$p/confout.gro
    cd ..
done

# Delete trajectories
echo "Cleanup"
for j in $(seq 100 $(($i - 10))); do
  p=$((0-$j))
  if [[ -f p_$p/traj_comp.xtc ]]; then
    rm p_$p/traj_comp.xtc
  fi
done
