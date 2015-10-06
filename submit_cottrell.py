#!/usr/bin/env python
"""Usage:
    submit_cottrell.py gaussian <node> <dir> <fname> (--bash | --direct) [--cores <nc>] [--dry]
    submit_cottrell.py lammps <node> <dir> <fname> (--bash | --direct) [--cores <nc>] [--dry]

A script to easily submit Gaussian and LAMMPS jobs to the Cottrell cluster

Arguments:
    <node>        Cluster node, from 0 to 10
    <dir>         Directory of the gjf file
    <fname>       Name of the gjf file, NO EXTENSION
    --direct      Submit directly using Gaussian form
    --bash        Submit via bash script cottrell.sh

Options:
    --cores <nc>  Number of cores [default: 16]
    --dry         Dry run, print file cmd on screen

pv278@cam.ac.uk, 16/05/15
"""
from docopt import docopt
from schema import Schema, And, Use, Optional
import os
import subprocess
import sys

def sup(node):
    """Get supervisor name"""
    if node < 10:
        return "jae"
    else:
        return "pdb"


if __name__ == "__main__":
    args = docopt(__doc__, version=0.1)
    schema = Schema({"<node>" : And(Use(int), lambda n: 0 <= n <= 19),
                     "<dir>"  : And(str, len),
                     "<fname>": And(str, len),
                     "--cores": And(Use(int), lambda n: 1 <= n <= 16),
                     "--dry"  : bool,
                     "--bash" : bool,
                     "--direct": bool,
                     "lammps"  : bool,
                     "gaussian": bool
                     })
    args = schema.validate(args)
#    print args
       
    node = args["<node>"]
    server = sup(node) + ".q@compute-0-" + str(node)
    cores = args["--cores"]

    ddir = args["<dir>"]
    filename  = args["<fname>"]
    bashscript = "~/GeneralScripts/cottrell.sh"

    if args["gaussian"]:
        prog = "gaussian"
        filepath = os.path.join(os.getcwd(), ddir, filename)
        filepath = filepath.rstrip(".xyz")
        if args["--bash"]:
            submit_string = "qsub -q " + server + " -pe orte " + str(cores) + \
                            " " + bashscript + " " + prog + " " + filepath
            if args["--dry"]:
                print submit_string
                sys.exit()
            subprocess.call(submit_string, shell=True)
        elif args["--direct"]:
            ext = ".gjf"
            base_dir = os.path.abspath(os.path.join(os.getcwd(), ".."))
            filedir = os.path.join(base_dir, dir)
            infilepath = filedir + "/" + filename + ext
            outfilepath = filedir + "/" + filename + ".out"
            
            exports = """
            g09root=\"/home/Gaussian/\"
            GAUSS_SCRDIR=\"/state/partition1/Gaussian_scratch\"
            GAUSS_EXEDIR=\"/home/Gaussian/g09/bsd:/home/Gaussian/g09/private:/home/Gaussian/g09\"
            export g09root GAUSS_SCRDIR GAUSS_EXEDIR
            . $g09root/g09/bsd/g09.profile"""
            subprocess.call(exports, shell=True)
            
            gaussianbin = "/home/Gaussian/g09/g09"
            cmd = gaussianbin + " < " + infilepath + " > " + outfilepath
            submit_string = "qsub -b y -q " + server + " -pe orte " + str(cores) + " " + cmd
 
            if args["--dry"]:
                print submit_string
                sys.exit()
            subprocess.call(submit_string, shell=True)

    elif args["lammps"]:
        prog = "lammps"
        filepath = os.path.join(os.getcwd(), ddir, filename)
        submit_string = "qsub -q " + server + " -pe orte " + str(cores) + \
                         " " + bashscript + " " + \
                         prog + " " + filepath + " " + str(cores)
        if args["--dry"]:
            print submit_string
            sys.exit()
        subprocess.call(submit_string, shell=True)
        


