#!/usr/bin/env python
"""
Collection of often used functions for
* input/output
* rotating, translating and printing data in xyz format

2nd version
pv278@cam.ac.uk, 17/06/15
"""
import numpy as np
from numpy.matlib import repmat
from math import *
import os


def read_Pt_spins(file="../Files/lowest_Pt_spins.txt"):
    """
    Read lowest spin states into a dictionary,
    either WA Goddard's or our results
    NOT USED NOW
    """
    states = [line.split(":") for line in open(file, "r")]
    states = dict(states)
    for i in states.keys():
        states[i] = map(int, states[i].split())
    return states    


class Atoms:
    """
    Class to produce and manipulate xyz coords, contains
    * atom names
    * atom coordinates
    * charge
    * spin state
    """
    def __init__(self, names=[], coords=np.array([]), charge=0, spin=1):
        self.names = list(names)
        self.coords = np.array(coords)
        self.charge = charge
        self.spin = spin

    def __len__(self):
        """Number of atoms"""
        return len(self.coords)

    def __repr__(self):
        """
        >>> from numpy import array
        >>> m = Atoms(["Ar"], [(0, 0, 0)])
        >>> repr(m)
        "Atoms(names=['Ar'], coords=array([[0, 0, 0]]), charge=0, spin=1)"
        >>> str(m) == str(eval(repr(m)))
        True
        """
        kwargs = [(attr, getattr(self, attr))
                  for attr in ['names', 'coords', 'charge', 'spin']
                  if hasattr(self, attr)]
        return "Atoms(%s)" % ', '.join('%s=%r' % kw for kw in kwargs)

    def __str__(self):
        """Print xyz onto screen"""
        if self.coords.any():
            M, N = self.coords.shape
            line = ""
            line += str(self.charge) + " " + str(self.spin) + "\n"
            for i in range(M):
                line += self.names[i] + "\t"
                for j in range(N):
                    line += "%.6f" % self.coords[i, j] + "\t"
                line += "\n"
            return line.rstrip("\n")
        else:
            return ""

    def __add__(self, other):
        return Atoms(self.names + other.names, \
                     np.vstack((self.coords, other.coords)), \
                     self.charge, \
                     self.spin)

    def save(self, filename, vmd=False):
        """
        Save xyz coords, charge and spin into file
        Also can save in VMD-compatible format w/ 1st line 
        with number of atoms and 2nd line a comment
        """
        with open(filename, "w") as f:
            M, N = self.coords.shape
            if vmd:
                f.write(str(M) + "\nBlabla\n")
            else:
                f.write(str(self.charge) + " " + str(self.spin) + "\n")
            for i in range(M):
                line = str(self.names[i]) + "\t"
                for j in range(N):
                    line += "%.6f" % self.coords[i, j] + "\t"
                line += "\n"
                f.write(line)
            f.close()
            print "Coords saved to", filename

    def read(self, fname):
        """Read the structure from an xyz file"""
        assert fname[-3:] == "xyz"
        with open(fname) as f:
            firstline = f.readline().split()
            if len(firstline) == 2:
                charge, spin = [int(s) for s in firstline]
            elif len(firstline) == 1:   # VMD format
                charge, spin = (0, 1)
                comment = f.readline()  # 2nd comment line
            M = np.array([line.split() for line in f.readlines()])
            names = M[:, 0]
            coords = M[:, 1:4].astype(np.float)
        return Atoms(names, coords, charge, spin)

    def shift(self, s):
        """Shift all atoms by a vector s"""
        assert len(s) == 3
        self.coords = self.coords + s    
        # self.coords += s RAISES ERROR! DeprecationWarning: Implicitly casting between incompatible kinds

    def rotate(self, theta=0, phi=0):
        """Rotate all atoms by angles theta and phi respectively"""
        N = len(self)
        Rtheta = np.array([[cos(theta),0,-sin(theta)],
                           [0,         1, 0         ],
                           [sin(theta),0, cos(theta)]])
        Rphi = np.array([[cos(phi),-sin(phi),0],
                         [sin(phi), cos(phi),0],
                         [0,        0,       1]])
        for i in range(N):
            self.coords[i, :] = np.dot(Rtheta, self.coords[i, :])
        for i in range(N):
            self.coords[i, :] = np.dot(Rphi, self.coords[i, :])

    def com(self):
        """Return centre of mass of atoms"""
        f = open(os.path.join(os.path.dirname(os.path.realpath(__file__)), "atomic_numbers.txt"), "r")
        weights = eval(f.read())
        N = len(self)
        com = np.zeros(3)
        for i in range(N):
            com += self.coords[i, :]*weights[self.names[i]]
        denom = sum([weights[i] for i in self.names])
        com /= denom
        return com

    def shift_com(self):
        """Shift all atoms to their centre of mass"""
        com = self.com()
        self.shift(-com)


if __name__ == "__main__":
    print "===== Testing the Atoms class ====="
    coords = np.array([1,0,0] + [0,1,0] + [0,0,1]).reshape((3, 3))
    atoms = Atoms(["Pt"]*3, coords)
    print atoms
    
    s = np.arange(3)
    print "* Shifting atoms by", s
    atoms.shift(s)
    print atoms

    atoms.coords = np.copy(coords)
    theta, phi = 90, 90
    print "* Rotating atoms by theta=%i, phi=%i:" % (theta, phi)
    atoms.rotate(radians(theta), radians(phi))
    print atoms

    atoms.coords = np.copy(coords)
    print "* Shifting atoms to its centre of mass"
    atoms.shift_com()
    print atoms


