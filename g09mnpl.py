#!/usr/bin/env python
"""Usage:
    g09mnpl.py gen_header [--nproc <np>] [--method <method>] [--basis <bs>] 
                          [--opt] [--freq] [--scfsteps <scf>]
    g09mnpl.py parse_coords <file> (--input | --final)
    g09mnpl.py gen_gjf from_xyz <files> [--header <header>] [--save]
    g09mnpl.py gen_gjf from_out <file>
    g09mnpl.py extract energies <files>  [--save <datafile>]

Gaussian 09 manipulate.
Collection of method to manipulate files entering or leaving Gaussian.
* construct gjf files
* parse xyz coordinates from outfiles
* generate headers
* extract data from outfiles

Options:
    --nproc <np>         Number of processors to use [default: 16]
    --method <method>    Qchem method [default: B3LYP]
    --basis <bs>         Basis set [default: 6-31g*]
    --opt                Produce "Opt" flag
    --freq               Produce "Freq" flag
    --scfsteps <scf>     Number of self-consistent steps [default: 1000]

pv278@cam.ac.uk, 06/10/15
"""
import numpy as np
import pandas as pd
import glob, os, sys, re
from docopt import docopt
from xyzlib import Atoms


def parse_input_coords(infile):
    """
    Parse the input atom names with xyz coords
    TODO: parse charge and multiplicity
    """
    regex = "Charge =(.*?)\\n \\n"
    with open(infile) as f:
        match = re.findall(regex, f.read(), re.S)   # WHAT IS re.S?
        match = match[-1].split("\n")[1:]
        match = [line.split() for line in match]
        names = [line[0] for line in match]
        coords = np.array([line[1:4] for line in match]).astype(float)
        return Atoms(names, coords)


def parse_last_coords(infile):
    """Parse last produced coords in a g09 simulation"""
    names = parse_input_coords(infile).names
    regex = "Standard orientation(.*?)Rotational constants"
    with open(infile) as f:
        match = re.findall(regex, f.read(), re.S)
        match = match[-1].split("\n")[5:-2]
        match = [line.split() for line in match]
        coords = np.array(match).astype(float)[:, -3:]
        return Atoms(names, coords)


def gen_header(p):
    """Generate Gaussian header, arguments: method, basis, opt"""
    header = "%nproc=" + str(p.get("np")) + "\n"
    header += "#T " + p.get("method") + "/" + p.get("basis") \
           + " Test " + p.get("opt") + " " + p.get("freq") \
           + " " + p.get("scf")
    header += "\n\n" + "Blabla" + "\n\n"
    return header


if __name__ == "__main__":
    args = docopt(__doc__)
#    print args
    default_header = "%nproc=16\n#T B3LYP/6-31G* Test\n\nBlabla\n\n"
    params = {}
    params["np"] = args["--nproc"]
    params["method"] = args["--method"]
    params["basis"] = args["--basis"]
    params["opt"] = "Opt" if args["--opt"] else ""
    params["freq"] = "Freq" if args["--freq"] else ""
    params["scf"] = "scf=(direct, maxcycle=%s)" % args["--scfsteps"] if args["--scfsteps"] else ""

    if args["gen_header"]:
        header = gen_header(params)
        print header
    
    if args["parse_coords"]:    #TO TEST
        infile = args["<file>"]
        if args["--input"]:
            A = parse_input_coords(infile)
            if args["--save"]:
                A.save(args["--save"], args["--vmd"])
            else:
                print A
        if args["--final"]:
            B = parse_last_coords(infile)
            if args["--save"]:
                B.save(args["--save"], args["--vmd"])
            else:
                print B

    elif args["gen_gjf"]:
        if args["from_xyz"]:
            infiles = glob.glob(args["<files>"])
            A = Atoms().read(infiles[0])
            for infile in infiles[1:]:
                A = A + Atoms().read(infile)
            string = default_header + str(A) + "\n\n"   # THINK ABOUT GENERATING HEADERS
            if args["--save"]:
                fname = infiles[0].rstrip("xyz") + "gjf"
                open(fname, "w").write(string)
                print "gjf file written into", fname
            else:
                print string


        if args["from_out"]:  #TO TEST
            infile = args["<files>"]
            assert(infile[-3:] == "out")
            A = parse_last_coords(infile)
            gen_g09_script(header, str(A), outfile)

    elif args["extract"] and args["energies"]:
        outfiles = glob.glob(args["<files>"])
        print outfiles
        positions, energies = [], []
        
        for outfile in outfiles:
            outfile = os.getcwd() + "/" + outfile
            E_line = [line for line in open(outfile, "r").readlines() if "SCF Done" in line]
            if E_line:       # if energy line exists
                E_line = E_line[0]
                energies.append(float(E_line.split()[4]))
                pos = re.search(r"(?<=_d)[^}]*(?=.out)", outfile)   # WHAT IS THIS DOING?
                pos = pos.group()
                positions.append(float(pos))
        
        if not energies:     # if the array is empty
            sys.exit()
        tbl = pd.DataFrame(energies, index=positions).sort_index()
        tbl.columns = ["E"]
        print tbl
        if args["<datafile>"]:
            datapath = args["<datafile>"]
            tbl.to_csv(datapath, sep="\t", header=False)
            print "Table saved in ", datapath


