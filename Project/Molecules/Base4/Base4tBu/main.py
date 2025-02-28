import subprocess
import os
from openbabel import openbabel
import math


def calculate_single_molecule_energy(optimized_xyz_file, charge, uhf, log_file_path):
    with open(log_file_path, 'a') as log:
        result = subprocess.run(
            ['xtb', optimized_xyz_file, '--sp', '--chrg', str(charge), '--uhf', str(uhf), '--spinpol', '--tblite'], 
            check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
        )
        log.write(result.stdout)
        energy_line = next(line for line in result.stdout.splitlines() if 'TOTAL ENERGY' in line)
        energy = float(energy_line.split()[3])
        log.write(f"Final Energy: {energy} Ha\n")
        return energy

def optimize_geometry(input_mol_file, output_mol_file, charge, uhf, log_file_path):
    with open(log_file_path, 'a') as log:
        subprocess.run(['xtb', input_mol_file, '--opt', 'tight', '--cycles', '500', '--chrg', str(charge), '--uhf', str(uhf), '--spinpol', '--tblite'], 
                       check=True, stdout=log, stderr=log)
        subprocess.run(['cp', 'xtbopt.mol', output_mol_file], check=True)
    return output_mol_file

def optimize_geometry_gfnff(input_xyz_file, output_mol_file, charge, uhf, log_file_path):
    with open(log_file_path, 'a') as log:
        log.write(f"Optimizing {input_xyz_file} using GFNFF\n")
        
        subprocess.run(['xtb', input_xyz_file, '--opt', '--cycles', '500', '--gfnff', '--chrg', str(charge), '--uhf', str(uhf)], 
                       check=True, stdout=log, stderr=log)
        
        subprocess.run(['cp', 'xtbopt.mol', output_mol_file], check=True)

    return output_mol_file


def stack_dimer(input_mol_file, stack_distance, rotation_angle):
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

    output_mol_file = "stacked_dimer.mol"
    obConversion.WriteFile(mol, output_mol_file)
    return output_mol_file

def calculate_reorganization_energy(neutral_file, anion_file, cation_file, 
                                    radical_monomer_charge, radical_monomer_uhf, 
                                    anion_monomer_charge, anion_monomer_uhf, 
                                    cation_monomer_charge, cation_monomer_uhf, 
                                    log_file_path):
    conversion_factor = 27.2114

    with open(log_file_path, 'a') as log:
        log.write(f"Calculating reorganization energy\n")

        # Calculate energies for electron reorganization energy
        E_neutral_optimized = calculate_single_molecule_energy(neutral_file, radical_monomer_charge, radical_monomer_uhf, log_file_path)
        E_anion_optimized = calculate_single_molecule_energy(anion_file, anion_monomer_charge, anion_monomer_uhf, log_file_path)
        E_neutral_in_anion_geom = calculate_single_molecule_energy(anion_file, radical_monomer_charge, radical_monomer_uhf, log_file_path)
        E_anion_in_neutral_geom = calculate_single_molecule_energy(neutral_file, anion_monomer_charge, anion_monomer_uhf, log_file_path)

        reorg_energy_electron = (E_neutral_in_anion_geom - E_neutral_optimized) + (E_anion_in_neutral_geom - E_anion_optimized)
        reorg_energy_electron *= conversion_factor

        # Calculate energies for hole reorganization energy
        E_cation_optimized = calculate_single_molecule_energy(cation_file, cation_monomer_charge, cation_monomer_uhf, log_file_path)
        E_neutral_in_cation_geom = calculate_single_molecule_energy(cation_file, radical_monomer_charge, radical_monomer_uhf, log_file_path)
        E_cation_in_neutral_geom = calculate_single_molecule_energy(neutral_file, cation_monomer_charge, cation_monomer_uhf, log_file_path)

        reorg_energy_hole = (E_neutral_in_cation_geom - E_neutral_optimized) + (E_cation_in_neutral_geom - E_cation_optimized)
        reorg_energy_hole *= conversion_factor

        # Log results
        log.write(f"Reorganization energy (electron): {reorg_energy_electron:.6f} eV\n")
        log.write(f"Reorganization energy (hole): {reorg_energy_hole:.6f} eV\n")

    return reorg_energy_electron, reorg_energy_hole


def calculate_transfer_integral(xyz_file, charge, uhf, log_file_path):
    result = subprocess.run(
        ['xtb', xyz_file, '--dipro', '0.3', '--chrg', str(charge), '--uhf', str(uhf), '--spinpol', '--tblite'],
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
                hole_transport = float(parts[-3])  # Extracting value
            except ValueError:
                log.write(f"Error: Could not parse hole transport from line: {line}\n")
                continue
        
        if 'charge transport (unocc. MOs)' in line:
            parts = line.split()
            try:
                charge_transport = float(parts[-3])  # Extracting value
            except ValueError:
                log.write(f"Error: Could not parse charge transport from line: {line}\n")
                continue

    if hole_transport is None or charge_transport is None:
        log.write("Error: Hole or charge transport values not found in the output.\n")
        return None, None
    else:
        return hole_transport, charge_transport

if __name__ == "__main__":
    input_mol_file = "Base4tBu.mol"
    radical_monomer_charge = 0
    radical_monomer_uhf = 1
    anion_monomer_charge = -1
    anion_monomer_uhf = 0
    cation_monomer_charge = 1
    cation_monomer_uhf = 0
    radical_dimer_charge = 0
    radical_dimer_uhf = 2
    stack_distance = 6
    log_file_path = "output_log.txt"
    
    with open(log_file_path, 'w') as log:
        log.write("Optimization of monomers started\n")
    
    opt_radical_monomer = optimize_geometry(input_mol_file, 'radical_monomer.mol', radical_monomer_charge, radical_monomer_uhf, log_file_path)
    opt_anion_monomer = optimize_geometry(input_mol_file, 'anion_monomer.mol', anion_monomer_charge, anion_monomer_uhf, log_file_path)
    opt_cation_monomer = optimize_geometry(input_mol_file, 'cation_monomer.mol', cation_monomer_charge, cation_monomer_uhf, log_file_path)
    
    reorg_energy_electron, reorg_energy_hole = calculate_reorganization_energy(opt_radical_monomer, opt_anion_monomer, opt_cation_monomer,radical_monomer_charge, radical_monomer_uhf, anion_monomer_charge, anion_monomer_uhf,cation_monomer_charge,cation_monomer_uhf, log_file_path)
    
    dimer = stack_dimer(opt_radical_monomer, stack_distance, rotation_angle=0)
    loose_opt_dimer = optimize_geometry_gfnff(dimer, 'loose_opt_dimer.mol', radical_dimer_charge, radical_dimer_uhf, log_file_path)
    opt_dimer = optimize_geometry(loose_opt_dimer, 'opt_dimer.mol', radical_dimer_charge, radical_dimer_uhf, log_file_path)
    
    J_hole, J_electron = calculate_transfer_integral(opt_dimer, radical_dimer_charge,radical_dimer_uhf, log_file_path)
    
    with open(log_file_path, 'a') as log:
        log.write(f"Reorganization energy (electron): {reorg_energy_electron} eV\n")
        log.write(f"Reorganization energy (hole): {reorg_energy_hole} eV\n")
        log.write(f"Electron Transfer integral J: {J_electron} eV\n")
        log.write(f"Hole Transfer integral J: {J_hole} eV\n")
