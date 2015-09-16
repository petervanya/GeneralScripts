#!/bin/bash

Usage(){
  echo "Usage: 
    $(basename $0) <prog> <file> <nc>
  
Script to submit Gaussian adsorption runs to the Cottrell

Arguments:
    <prog>        Programme to use, 'gaussian' or 'lammps'
    <file>        Absolute path to the file with NO extension"
exit 0
}

if [ "$1" = "-h" ]; then
  Usage
fi

if [ -z "$2" ]; then
  Usage
else
  filepath=$2
fi

if [ -z "$3" ]; then
  Ncores=16
else
  Ncores=$3
fi

# Generic Cottrell header
# =======================
#$ -cwd                                                                         
#$ -j y                                                                         
#$ -S /bin/bash                                                                 

# Gaussian commands                                                             
g09root="/home/Gaussian"
GAUSS_SCRDIR="/state/partition1/Gaussian_scratch"
GAUSS_EXEDIR="/home/Gaussian/g09/bsd:/home/Gaussian/g09/private:/home/Gaussian/g09"
export g09root GAUSS_SCRDIR GAUSS_EXEDIR
. $g09root/g09/bsd/g09.profile
# =======================
#/home/gaussian/g09/g09 < $filepath.gjf > $filepath.out

if [ "$1" == "gaussian" ]; then
  progpath="/home/Gaussian/g09/g09"
  $progpath < $filepath.gjf > $filepath.out
elif [ "$1" == "lammps" ]; then
  progpath="/home/pv278/sw/bin/lmp_mpi"
  mpi_run="/opt/openmpi/bin/mpirun"
  $mpi_run -np $Ncores $progpath < $filepath
else
  Usage
fi
