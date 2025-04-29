import subprocess
import os
import numpy as np
import pandas as pd
import math
from openbabel import OBMol, OBConversion, OBMolAtomIter


def stack_dimer(input_mol_file, stack_distance, rotation_angle=0):
    obConversion = OBConversion()
    obConversion.SetInAndOutFormats("mol", "xyz")

    mol = OBMol()
    obConversion.ReadFile(mol, input_mol_file)

    mol_copy = OBMol(mol)

    rad = math.radians(rotation_angle)
    cos_theta = math.cos(rad)
    sin_theta = math.sin(rad)

    for atom in OBMolAtomIter(mol_copy):
        x, y, z = atom.GetX(), atom.GetY(), atom.GetZ()
        new_x = x * cos_theta - y * sin_theta
        new_y = x * sin_theta + y * cos_theta
        new_z = z + stack_distance
        atom.SetVector(new_x, new_y, new_z)

    mol += mol_copy
    mol.PerceiveBondOrders()

    output_xyz = f"dimer_{stack_distance:.2f}.xyz"
    obConversion.WriteFile(mol, output_xyz)

    return output_xyz


def calculate_transfer_integral(xyz_file, charge, uhf, log_file):
    result = subprocess.run(
        ['xtb', xyz_file, '--dipro', '0.3', '--chrg', str(charge), '--uhf', str(uhf), '--input', 'frag.inp', '--spinpol', '--tblite'],
        capture_output=True, text=True
    )

    log_file.write(result.stdout)
    log_file.write(result.stderr)

    hole_transport = None
    electron_transport = None

    # Robust parsing
    for line in result.stdout.splitlines():
        if "total |J(AB,eff)| for hole transport (occ. MOs)" in line:
            try:
                hole_transport = float(
                    line.split("hole transport (occ. MOs) :")[1].split()[0]
                )
            except Exception as e:
                log_file.write(f"Error parsing hole_transport: {e}\n")

        if "total |J(AB,eff)| for charge transport (unocc. MOs)" in line:
            try:
                electron_transport = float(line.split("charge transport (unocc. MOs) :")[1].split()[0])
            except Exception as e:
                log_file.write(f"Error parsing electron_transport: {e}\n")

    if hole_transport is None and electron_transport is None:
        raise ValueError(f"Could not parse J values from {xyz_file}")

    return hole_transport, electron_transport



if __name__ == "__main__":
    radical_dimer_charge = 0
    radical_dimer_uhf = 2
    input_mol_file = "radical_monomer.mol"

    distances = [round(x, 2) for x in np.arange(4, 10, 0.05)]

    results = []

    with open("calculation_log.txt", "w") as log_file:
        for dist in distances:
            log_file.write(f"Stacking distance: {dist} Å\n")
            stacked_xyz_file = stack_dimer(input_mol_file, dist)

            try:
                hole_J, electron_J = calculate_transfer_integral(
                    stacked_xyz_file, radical_dimer_charge, radical_dimer_uhf, log_file
                )
                results.append({
                    'Distance (Å)': dist,
                    'Hole Transfer Integral (eV)': hole_J if hole_J is not None else 'NaN',
                    'Electron Transfer Integral (eV)': electron_J if electron_J is not None else 'NaN'
                })

                log_file.write(
                    f"Calculated J at {dist} Å: Hole={hole_J if hole_J else 'NaN'} eV, "
                    f"Electron={electron_J if electron_J else 'NaN'} eV\n"
                )

            except Exception as e:
                log_file.write(f"Calculation failed for distance {dist} Å: {e}\n")

        df = pd.DataFrame(results)
        df_output = "transfer_integrals_vs_distance.csv"
        df.to_csv(df_output, index=False)
        log_file.write(f"Results saved in {df_output}\n")

    print(f"Results saved in {df_output}")
