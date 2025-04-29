#!/usr/bin/env python3
import sys

BOHR_TO_ANGSTROM = 0.52917721067

def get_lines(filename):
    with open(filename, "r") as f:
        lines = f.readlines()
    return lines

def coordlines_to_data(line):
    tokens = line.split()
    x = float(tokens[0]) * BOHR_TO_ANGSTROM
    y = float(tokens[1]) * BOHR_TO_ANGSTROM
    z = float(tokens[2]) * BOHR_TO_ANGSTROM
    atom_type = tokens[3].upper()
    return (atom_type, x, y, z)

if __name__ == "__main__":
    coord_filename = sys.argv[1]
    lines = get_lines(coord_filename)

    atom_count = len(lines) - 2
    output_filename = coord_filename + ".xyz"

    with open(output_filename, "w") as out:
        out.write(f"{atom_count}\n")
        out.write(f".xyz output generated from file: {coord_filename}\n")

        for line in lines[1:-1]:
            (atom_type, x, y, z) = coordlines_to_data(line)
            out.write(" %-2s    %20.12f %20.12f %20.12f\n" % (atom_type, x, y, z))

    print(f"Wrote XYZ file to: {output_filename}")
