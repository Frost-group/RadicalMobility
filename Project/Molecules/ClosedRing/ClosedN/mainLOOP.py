import subprocess
import os
from openbabel import openbabel


import subprocess

def optimize_geometry(input_mol_file, output_mol_file, charge, multiplicity, log_file_path):
    # Optimize the geometry using xtb and log the output
    with open(log_file_path, 'a') as log:  # Append to log file
        log.write(f"Running optimization for {input_mol_file}\n")
        
        # Run the xtb optimization with specified charge and multiplicity
        subprocess.run(['xtb', input_mol_file, '--opt', '--gfn', '2', '--chrg', str(charge), '--uhf', str(multiplicity-1), '--tight'], 
                       check=True, stdout=log, stderr=log)
        
        # Ensure the optimized geometry is saved as a unique .mol file
        subprocess.run(['cp', 'xtbtopo.mol', output_mol_file], check=True)
    
    return output_mol_file

def optimize_geometry_gfnff(input_xyz_file, output_mol_file, charge, multiplicity, log_file_path):
    # Optimize the geometry using xtb and log the output with GFNFF
    with open(log_file_path, 'a') as log:  # Append to log file
        log.write(f"Running optimization for {input_xyz_file} using GFNFF\n")
        
        # Run the xtb optimization with specified charge, multiplicity and using the GFNFF method
        subprocess.run(['xtb', input_xyz_file, '--opt', '--gfnff', '--chrg', str(charge), '--uhf', str(multiplicity-1)], 
                       check=True, stdout=log, stderr=log)
        
        # Ensure the optimized geometry is saved as a unique .mol file
        subprocess.run(['cp', 'xtbopt.mol', output_mol_file], check=True)
    
    return output_mol_file

def calculate_single_molecule_energy(optimized_xyz_file, charge, multiplicity, log_file_path):
    with open(log_file_path, 'a') as log:
        log.write(f"Running single point energy calculation for {optimized_xyz_file}\n")
        result = subprocess.run(
            ['xtb', optimized_xyz_file, '--gfn2', '--sp', '--chrg', str(charge), '--uhf', str(multiplicity - 1)], 
            check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
        )
        log.write(result.stdout)  # Log the output
        energy_line = next(line for line in result.stdout.splitlines() if 'TOTAL ENERGY' in line)
        energy = float(energy_line.split()[3])
        log.write(f"Energy for {optimized_xyz_file}: {energy}\n")
        return energy


def calculate_reorganization_energy(neutral_file, anion_file, cation_file, log_file_path):
    conversion_factor = 27.2114
    with open(log_file_path, 'a') as log:
        log.write(f"Calculating reorganization energy\n")
        E_neutral_optimized = calculate_single_molecule_energy(neutral_file, radical_monomer_charge, radical_monomer_multiplicity, log_file_path)
        E_anion_optimized = calculate_single_molecule_energy(anion_file, anion_monomer_charge, anion_monomer_multiplicity, log_file_path)
        E_neutral_in_anion_geom = calculate_single_molecule_energy(anion_file, radical_monomer_charge, anion_monomer_multiplicity, log_file_path)
        E_anion_in_neutral_geom = calculate_single_molecule_energy(neutral_file, anion_monomer_charge, radical_monomer_multiplicity, log_file_path)
        reorg_energy_electron_hartree = (E_neutral_in_anion_geom - E_neutral_optimized) + (E_anion_in_neutral_geom - E_anion_optimized)
        reorg_energy_electron_eV = reorg_energy_electron_hartree * conversion_factor
        E_cation_optimized = calculate_single_molecule_energy(cation_file,cation_monomer_charge, cation_monomer_multiplicity, log_file_path)
        E_neutral_in_cation_geom = calculate_single_molecule_energy(cation_file, radical_monomer_charge, cation_monomer_multiplicity, log_file_path)
        E_cation_in_neutral_geom = calculate_single_molecule_energy(neutral_file, cation_monomer_charge, radical_monomer_multiplicity, log_file_path)
        reorg_energy_hole_hartree = (E_neutral_in_cation_geom - E_neutral_optimized) + (E_cation_in_neutral_geom - E_cation_optimized)
        reorg_energy_hole_eV = reorg_energy_hole_hartree * conversion_factor
    return reorg_energy_electron_eV, reorg_energy_hole_eV

def rotate_molecule(mol, angles):
    """
    Rotates the molecule by given angles around all three axes (x, y, z).
    """
    rad_angles = np.radians(angles)
    cos_a, sin_a = np.cos(rad_angles), np.sin(rad_angles)
    
    rot_x = np.array([[1, 0, 0], [0, cos_a[0], -sin_a[0]], [0, sin_a[0], cos_a[0]]])
    rot_y = np.array([[cos_a[1], 0, sin_a[1]], [0, 1, 0], [-sin_a[1], 0, cos_a[1]]])
    rot_z = np.array([[cos_a[2], -sin_a[2], 0], [sin_a[2], cos_a[2], 0], [0, 0, 1]])
    
    rot_matrix = rot_x @ rot_y @ rot_z  # Combined rotation matrix
    
    for atom in openbabel.OBMolAtomIter(mol):
        x, y, z = atom.GetX(), atom.GetY(), atom.GetZ()
        new_coords = np.dot(rot_matrix, np.array([x, y, z]))
        atom.SetVector(*new_coords)

def stack_dimer(input_mol_file, translation, rotation_angles):
    """
    Stacks a molecule by translating along (x, y, z) and rotating by (x, y, z) angles.
    """
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("mol", "mol")
    
    mol = openbabel.OBMol()
    obConversion.ReadFile(mol, input_mol_file)
    
    mol_copy = openbabel.OBMol(mol)
    
    # Apply rotation
    if any(angle != 0 for angle in rotation_angles):
        rotate_molecule(mol_copy, rotation_angles)
    
    # Apply translation
    for atom in openbabel.OBMolAtomIter(mol_copy):
        x, y, z = atom.GetX(), atom.GetY(), atom.GetZ()
        atom.SetVector(x + translation[0], y + translation[1], z + translation[2])
    
    mol += mol_copy
    mol.PerceiveBondOrders()
    
    output_mol_file = "stacked_dimer.mol"
    obConversion.WriteFile(mol, output_mol_file)
    
    return output_mol_file

def calculate_transfer_integral(xyz_file, charge, multiplicity, log_file_path):
    result = subprocess.run(
        ['xtb', xyz_file, '--dipro', '0.3', '--gfn', '2', '--chrg', str(charge), '--uhf', str(multiplicity-1)],
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
                hole_transport = float(parts[-3])  # Modified the index to match the value
            except ValueError:
                log.write(f"Error: Could not parse hole transport from line: {line}\n")
                continue
        
        if 'charge transport (unocc. MOs)' in line:
            parts = line.split()
            try:
                charge_transport = float(parts[-3])  # Modified the index to match the value
            except ValueError:
                log.write(f"Error: Could not parse charge transport from line: {line}\n")
                continue

    if hole_transport is None or charge_transport is None:
        log.write("Error: Hole or charge transport values not found in the output.\n")
        return None, None
    else:
        return hole_transport, charge_transport

def calculate_transfer_rate(J, reorg_energy, temperature=300):
    """
    Calculates the rate of electron/hole transfer using Marcus theory.
    """
    h_bar = 4.135667696E-15  # Planck's constant (eV·s)
    k_B = 0.000086173  # Boltzmann constant (eV/K)
    lambda_ = reorg_energy + 0.3  # Adding empirical correction term
    
    prefactor = (J**2) / h_bar
    denominator = np.sqrt(np.pi * lambda_ * k_B * temperature)
    exponent = np.exp(-lambda_ / (4 * k_B * temperature))
    
    rate = prefactor * (1 / denominator) * exponent
    return rate

def calculate_mobility(transfer_rate, displacement=1e-8, temperature=300):
    """
    Calculates charge carrier mobility using the diffusion relation: D = k * A^2 and the Einstein relation.
    """
    k_B = 0.000086173  # Boltzmann constant (eV/K)
    diffusion_coefficient = transfer_rate * displacement**2  # D = k * A^2
    mobility = diffusion_coefficient / (k_B * temperature)  # μ = D / (k_B * T)
    return mobility



import csv
import os
import numpy as np

def log_to_text(log_file, message):
    """Append a log message to the text log file."""
    with open(log_file, 'a') as log:
        log.write(message + '\n')

if __name__ == "__main__":
    input_mol_file = "ClosedN.mol"
    
    radical_monomer_charge = 0
    radical_monomer_multiplicity = 1
    anion_monomer_charge = -1
    anion_monomer_multiplicity = 1
    cation_monomer_charge = 1
    cation_monomer_multiplicity = 1
    radical_dimer_charge = 0
    radical_dimer_multiplicity = 1
    displacement = 4e-10  # Distance for diffusion calculation
    
    log_txt_path = "output_log.txt"
    log_csv_path = "molecular_data_log.csv"
    header = ["Molecule", "dx", "dy", "dz", "rx", "ry", "rz", "Electron Rate (s^-1)", "Hole Rate (s^-1)"]

    # Create the CSV file with headers if it does not exist
    if not os.path.exists(log_csv_path):
        with open(log_csv_path, 'w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(header)

    # Initialize the text log file
    with open(log_txt_path, 'w') as log:
        log.write("Optimization and calculation log:\n\n")

    with open(log_csv_path, 'a', newline='') as file:  # Append data to CSV
        writer = csv.writer(file)

        opt_radical_monomer = optimize_geometry(input_mol_file, 'radical_monomer.mol', radical_monomer_charge, radical_monomer_multiplicity, log_txt_path)
        opt_anion_monomer = optimize_geometry(input_mol_file, 'anion_monomer.mol', anion_monomer_charge, anion_monomer_multiplicity, log_txt_path)
        opt_cation_monomer = optimize_geometry(input_mol_file, 'cation_monomer.mol', cation_monomer_charge, cation_monomer_multiplicity, log_txt_path)
        reorg_energy_electron, reorg_energy_hole = calculate_reorganization_energy(opt_radical_monomer, opt_anion_monomer, opt_cation_monomer, log_txt_path)
        
        for dx in range(4, 5, 1):
            for dy in range(4, 5, 1):
                for dz in range(4, 5, 1):
                    for rx in range(0, 1, 1):
                        for ry in range(0, 1, 1):
                            for rz in range(0, 1, 1):
                                translation = [dx, dy, dz]
                                rotation_angles = [rx, ry, rz]
                                
                                log_to_text(log_txt_path, f"Processing: Translation {translation}, Rotation {rotation_angles}")
                                dimer = stack_dimer(opt_radical_monomer, translation, rotation_angles)
                                
                                J_hole, J_electron = calculate_transfer_integral(dimer, radical_dimer_charge, radical_dimer_multiplicity, log_txt_path)
                                
                                rate_electron = calculate_transfer_rate(J_electron, reorg_energy_electron)
                                rate_hole = calculate_transfer_rate(J_hole, reorg_energy_hole)
                                
                                # Write results to text log
                                log_to_text(log_txt_path, f"Electron Rate: {rate_electron} s^-1")
                                log_to_text(log_txt_path, f"Hole Rate: {rate_hole} s^-1")
                                
                                # Write final numerical results to CSV file
                                writer.writerow([input_mol_file, dx, dy, dz, rx, ry, rz, rate_electron, rate_hole])

                                log_to_text(log_txt_path, f"Logged to database: {input_mol_file}, {dx}, {dy}, {dz}, {rx}, {ry}, {rz}, {rate_electron}, {rate_hole}\n")
