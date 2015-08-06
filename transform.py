#!/usr/bin/env python
"""Usage:
    transform.py <file> [-c <c>] [-a <a>] [-r <r>] [-s <s>] [-f <f>]
                      [--save <out>]

Centre, shift, align or rotate a molecule.

Arguments:
    <file>                Input file

Options:
    -c <c>,--centre <c>   Atom number (1-N) to centre the molecule around.
    -a <a>,--align <a>    Align s.t. given two atoms "n1 n2" lie on x-axis
    -r <r>,--rotate <r>   Rotate by angles theta and phi
    -s <s>,--shift <s>    Final shift of all atoms by a "x y z"
    -f <f>,--flip <f>     Flip atoms "n1 n2" by rotating whole molecule
    --save <out>          Save final coords into <out>  

pv278@cam.ac.uk 16/04/15
"""
from docopt import docopt
import numpy as np
from numpy.linalg import norm
from numpy.matlib import repmat
from scipy.linalg import expm
from math import *
from iolib import save_xyz, read_xyz, print_xyz, rotate_theta, rotate_phi, shift

def rotate(A, omega, alpha):
    """rotate coordinates A around vector omega by angle alpha"""
    R = np.matrix([[0,        -omega[2], omega[1]],
                   [ omega[2], 0,       -omega[0]],
                   [-omega[1], omega[0], 0       ]])
    R *= alpha
  
    for i in range(N):
        A[i,:] = np.dot(A[i,:], expm(R))
    return A

def flip(A, n1, n2):
    """flip atoms at positions n1, n2 by rotating"""
    centre = A[n1-1,:]
    dist = A[n2-1,:] - A[n1-1,:]
    A = shift(A, -(centre + dist/2))
    omega = np.array([0, dist[2], -dist[1]])
    A = rotate(A, omega, pi)
    A = shift(A, centre + dist/2)
    return A

def align(A, n1, n2):
    """align the coords so that atoms n1, n2 lie on x-axis"""
    if (n1 not in range(1,N+1)) or (n2 not in range(1,N+1)):
        raise ValueError
    vec12 = A[n1-1,:] - A[n2-1,:]
    axis = np.array([1,0,0])
    alpha = acos( np.dot(vec12,axis)/norm(vec12) )
    omega = np.cross(vec12, axis)
    if norm(omega) < 1e-12:
        print "Atoms",n1,"and",n2,"already aligned"
        return A

    omega /= norm(omega)
    R = np.matrix([[0,        -omega[2], omega[1]],
                   [ omega[2], 0,       -omega[0]],
                   [-omega[1], omega[0], 0       ]])
    R *= -alpha
  
    for i in range(N):
        A[i,:] = np.dot(A[i,:],expm(R))
    print "Atoms",n1,"and",n2,"aligned on x-axis"
    return A
    

if __name__ == "__main__":
    args = docopt(__doc__,version=1.0)
#    print args
    filename = args["<file>"]
    
    names, A = read_xyz(filename)
    N = len(A)

    if args["--centre"]:
         n = int(args["--centre"])
         A = shift(A, -A[n-1,:])
   
    if args["--align"]:
        n1, n2 = [int(i) for i in args["--align"].split()]
        if n1 != n2:
            A = align(A, n1, n2)

    if args["--flip"]:
        n1, n2 = [int(i) for i in args["--flip"].split()]
        if n1 != n2:
            A = flip(A, n1, n2)
    
    if args["--rotate"]:
        theta, phi = [radians(float(i)) for i in args["--rotate"].split()]
        A = rotate_theta(A, theta)
        print "Rotated by theta =",degrees(theta)
        A = rotate_phi(A, phi)
        print "Rotated by phi =",degrees(phi)
    
    if args["--shift"]:
      s = [float(i) for i in args["--shift"].split()]
      A = shift(A, s)
    
    if args["--save"]:
        fname = args["--save"]
        save_xyz(A, names, fname)
    else:
        print_xyz(A, names)
   

