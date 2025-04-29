import subprocess
import os
import numpy as np
import pandas as pd
from openbabel import OBMol, OBConversion, OBMolAtomIter


def stack_dimer(input_mol_file, stack_distance):
    obConversion = OBConversion()
    obConversion.SetInAndOutFormats("mol", "xyz")

    mol = OBMol()
    obConversion.ReadFile(mol, input_mol_file)

    mol_copy = OBMol(mol)

    # Translate second molecule along Z
    for atom in OBMolAtomIter(mol_copy):
        x, y, z = atom.GetX(), atom.GetY(), atom.GetZ()
        atom.SetVector(x, y, z + stack_distance)

    # Combine molecules
    mol += mol_copy
    mol.PerceiveBondOrders()

    output_xyz = f"dimer_initial_{stack_distance:.2f}.xyz"
    obConversion.WriteFile(mol, output_xyz)

    return output_xyz


def optimize_geometry(xyz_file, charge, uhf, log_file):
    optimized_xyz = xyz_file.replace("initial", "optimized")
    result = subprocess.run(
        ['xtb', xyz_file, '--opt', '--chrg', str(charge), '--uhf', str(uhf), '--gfnff'],
        capture_output=True, text=True
    )

    with open(log_file, "a") as log:
        log.write(result.stdout)
        log.write(result.stderr)

    # Rename optimized file from xtb output
    os.rename("xtbopt.xyz", optimized_xyz)

    return optimized_xyz


def translate_and_save(xyz_file, translation, new_xyz_file):
    obConversion = OBConversion()
    obConversion.SetInAndOutFormats("xyz", "xyz")

    mol = OBMol()
    obConversion.ReadFile(mol, xyz_file)

    atoms = list(OBMolAtomIter(mol))
    n_atoms = len(atoms) // 2  # assuming dimer

    # Translate only second molecule
    for atom in atoms[n_atoms:]:
        x, y, z = atom.GetX(), atom.GetY(), atom.GetZ()
        atom.SetVector(x, y, z + translation)

    obConversion.WriteFile(mol, new_xyz_file)

def calculate_transfer_integral(xyz_file, charge, uhf, log_file):
    result = subprocess.run(
        ['xtb', xyz_file, '--dipro', '0.3', '--chrg', str(charge), '--uhf', str(uhf), '--spinpol', '--tblite', '--input', 'frag.inp'],
        capture_output=True, text=True
    )

    with open(log_file, "a") as log:
        log.write(result.stdout)
        log.write(result.stderr)

    hole_transport = electron_transport = None

    for line in result.stdout.splitlines():
        if "total |J(AB,eff)| for hole transport" in line:
            try:
                hole_transport = float(line.split(":")[-2].strip().split()[0])
            except (IndexError, ValueError) as e:
                print(f"Failed to parse hole J from line: {line.strip()} — {e}")
                hole_transport = None
        elif "total |J(AB,eff)| for charge transport" in line:
            try:
                electron_transport = float(line.split(":")[-2].strip().split()[0])
            except (IndexError, ValueError) as e:
                print(f"Failed to parse electron J from line: {line.strip()} — {e}")
                electron_transport = None

    return hole_transport, electron_transport


if __name__ == "__main__":
    charge, uhf = 0, 2
    input_mol_file = "radical_monomer.mol"
    initial_distance = 6.0

    with open("calculation_log.txt", "w") as log_file:
        # Initial stacking and optimization
        initial_xyz = stack_dimer(input_mol_file, initial_distance)
        optimized_xyz = optimize_geometry(initial_xyz, charge, uhf, "calculation_log.txt")

        # Scan distances around optimized structure
        scan_range = np.arange(0, 5.0, 0.05)
        results = []

        for delta in scan_range:
            scan_distance = initial_distance + delta
            scan_xyz_file = f"dimer_scan_{scan_distance:.2f}.xyz"
            translate_and_save(optimized_xyz, delta, scan_xyz_file)

            hole_J, electron_J = calculate_transfer_integral(
                scan_xyz_file, charge, uhf, "calculation_log.txt"
            )

            results.append({
                'Distance (Å)': scan_distance,
                'Hole Transfer Integral (eV)': hole_J,
                'Electron Transfer Integral (eV)': electron_J
            })

            log_file.write(f"Distance: {scan_distance} Å, Hole J: {hole_J} eV, Electron J: {electron_J} eV\n")

        df = pd.DataFrame(results)
        df_output = "optimized_transfer_integrals_vs_distance.csv"
        df.to_csv(df_output, index=False)
        log_file.write(f"Results saved in {df_output}\n")

    print(f"Calculation completed. Results saved in {df_output}")
