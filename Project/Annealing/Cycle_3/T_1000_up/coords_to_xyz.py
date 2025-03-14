#!/usr/bin/env python2
import sys

BOHR_TO_ANGSTROM = 0.52917721067

def get_lines(filename):
    """ Reads lines from the given filename and returns them as a list. """
    with open(filename, "r") as f:
        lines = f.readlines()
    return lines

def coordlines_to_data(line):
    """ Converts coordinate lines from bohr to angstrom and extracts atom type. """
    tokens = line.split()

    x = float(tokens[0]) * BOHR_TO_ANGSTROM
    y = float(tokens[1]) * BOHR_TO_ANGSTROM
    z = float(tokens[2]) * BOHR_TO_ANGSTROM
    atom_type = tokens[3].upper()

    return (atom_type, x, y, z)

if __name__ == "__main__":

    coord_filename = "scoord.1"
    output_filename = "output.xyz"

    lines = get_lines(coord_filename)

    num_atoms = len(lines) - 2  # Exclude the first and last line

    with open(output_filename, "w") as xyz_file:
        xyz_file.write("{}\n".format(num_atoms))
        xyz_file.write(".xyz output generated from file: {}\n".format(coord_filename))

        for line in lines[1:-1]:  # Exclude first and last lines
            atom_type, x, y, z = coordlines_to_data(line)
            xyz_file.write("{:2s}    {:20.12f} {:20.12f} {:20.12f}\n".format(atom_type, x, y, z))

    print("XYZ file successfully written to {}".format(output_filename))
