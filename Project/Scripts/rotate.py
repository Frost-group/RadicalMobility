import os
import subprocess
import math
import re
import numpy as np
from openbabel import openbabel

# ========== CONFIGURABLE PARAMETERS ==========
rotation_angle = 30  # Rotation angle in degrees (Modify this if needed)
monomer_file = "radical_monomer.mol"
log_file_path = "J_values_log.txt"

# ========== FUNCTION DEFINITIONS ==========
def stack_dimer(input_mol_file, stack_distance, rotation_angle, output_mol_file):
    """Stacks the radical monomer at the specified distance and rotation angle."""
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("mol", "mol")

    mol = openbabel.OBMol()
    obConversion.ReadFile(mol, input_mol_file)

    mol_copy = openbabel.OBMol(mol)
    angle_rad = math.radians(rotation_angle)
    cos_theta = math.cos(angle_rad)
    sin_theta = math.sin(angle_rad)

    for atom in openbabel.OBMolAtomIter(mol_copy):
        x, y, z = atom.GetX(), atom.GetY(), atom.GetZ()
        new_x = cos_theta * x - sin_theta * y
        new_y = sin_theta * x + cos_theta * y
        new_z = z + stack_distance
        atom.SetVector(new_x, new_y, new_z)

    mol += mol_copy
    mol.PerceiveBondOrders()

    obConversion.WriteFile(mol, output_mol_file)
    return output_mol_file

def calculate_transfer_integral(xyz_file, charge, uhf, log_file_path):
    """Runs xTB DIPRO to calculate transfer integrals."""
    result = subprocess.run(
        ['xtb', xyz_file, '--dipro', '0.3', '--chrg', str(charge), '--uhf', str(uhf)],
        check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )

    with open(log_file_path, 'a') as log:
        log.write(result.stdout)  

    hole_transport = None
    charge_transport = None

    for line in result.stdout.splitlines():
        if 'hole transport (occ. MOs)' in line:
            parts = line.split()
            try:
                hole_transport = float(parts[-3])
            except ValueError:
                log.write(f"Error: Could not parse hole transport from line: {line}\n")
                continue

        if 'charge transport (unocc. MOs)' in line:
            parts = line.split()
            try:
                charge_transport = float(parts[-3])
            except ValueError:
                log.write(f"Error: Could not parse charge transport from line: {line}\n")
                continue

    return hole_transport, charge_transport

# ========== MAIN SCRIPT ==========
if __name__ == "__main__":
    radical_dimer_charge = -2
    radical_dimer_uhf = 4

    with open(log_file_path, 'w') as log:
        log.write("Transfer Integral Calculation Log\n")
        log.write(f"Using rotation angle: {rotation_angle}°\n\n")

    current_directory = os.getcwd()
    directories = [d for d in os.listdir(current_directory) if os.path.isdir(d)]

    for directory in directories:
        # Extract stacking distance from directory name
        match = re.search(r'(\d+\.?\d*)', directory)
        if match:
            stack_distance = float(match.group(1))
            print(f"\nProcessing directory: {directory} with stack distance {stack_distance} Å")

            os.chdir(directory)  # Move into the directory

            # Stack the dimer
            stacked_dimer_file = f"stacked_dimer_{stack_distance}.mol"
            stack_dimer(monomer_file, stack_distance, rotation_angle, stacked_dimer_file)

            # Run xTB DIPRO to calculate transfer integrals
            J_hole, J_electron = calculate_transfer_integral(stacked_dimer_file, radical_dimer_charge, radical_dimer_uhf, log_file_path)

            # Save results to the log file
            with open(log_file_path, 'a') as log:
                log.write(f"Stack distance: {stack_distance} Å\n")
                log.write(f"Electron Transfer integral J: {J_electron} eV\n")
                log.write(f"Hole Transfer integral J: {J_hole} eV\n\n")

            os.chdir(current_directory)  # Return to the main directory

    print("\nAll calculations completed. Check J_values_log.txt for results.")
