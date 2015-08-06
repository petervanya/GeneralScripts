#!/bin/bash

usage(){
  echo "Usage: 
    $(basename $0) <dir> <file>
  
Script to submit Gaussian adsorption runs to the Cottrell

Arguments:
    1  directory name, e.g. "Water" or "Plain" + "Pt9_10_9/Eta_*/S_*"
    2  file name without extension" 
exit 0
}

if [ "$1" = "-h" ]; then
  usage
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

dir=$1
file=$2

cd /home/pv278/Platinum/$dir
/home/Gaussian/g09/g09 < $file.gjf > $file.out
