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
        subprocess.run(['xtb', input_xyz_file, '--opt', '--gfnff', '--chrg', str(charge), '--uhf', str(multiplicity-1), '--loose'], 
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


def stack_dimer(input_mol_file, stack_distance):

    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("mol", "mol")

    mol = openbabel.OBMol()
    obConversion.ReadFile(mol, input_mol_file)

    mol_copy = openbabel.OBMol(mol)
    for atom in openbabel.OBMolAtomIter(mol_copy):
        x, y, z = atom.GetX(), atom.GetY(), atom.GetZ()
        atom.SetVector(x, y, z + stack_distance)

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


if __name__ == "__main__":
    input_mol_file = "Katherine3.mol"  # Changed to .xyz

    radical_monomer_charge = 0
    radical_monomer_multiplicity = 2
    anion_monomer_charge = -1
    anion_monomer_multiplicity = 1
    cation_monomer_charge = 1
    cation_monomer_multiplicity = 1
    radical_dimer_charge = 0
    radical_dimer_multiplicity = 3
    stack_distance = 6
    
    log_file_path = "output_log.txt"
    
    with open(log_file_path, 'w') as log:  # Write the log header
        log.write("Optimization of monomers started\n")
        opt_radical_monomer = optimize_geometry(input_mol_file, 'radical_monomer.mol', radical_monomer_charge, radical_monomer_multiplicity, log_file_path)
        opt_anion_monomer = optimize_geometry(input_mol_file, 'anion_monomer.mol', anion_monomer_charge, anion_monomer_multiplicity, log_file_path)
        opt_cation_monomer = optimize_geometry(input_mol_file, 'cation_monomer.mol', cation_monomer_charge, cation_monomer_multiplicity, log_file_path)
        
        reorg_energy_electron, reorg_energy_hole = calculate_reorganization_energy(opt_radical_monomer, opt_anion_monomer, opt_cation_monomer, log_file_path)
        
        log.write("Stacking dimer\n")
        dimer = stack_dimer(opt_radical_monomer, stack_distance)
        loose_opt_dimer = optimize_geometry_gfnff(dimer, 'loose_opt_dimer.mol', radical_dimer_charge, radical_dimer_multiplicity, log_file_path)
        opt_dimer = optimize_geometry(loose_opt_dimer, 'opt_dimer.mol', radical_dimer_charge, radical_dimer_multiplicity, log_file_path)
        
        log.write("Calculating transfer integral\n")
        J_hole, J_electron = calculate_transfer_integral(opt_dimer, radical_dimer_charge, radical_dimer_multiplicity, log_file_path)
        
        log.write(f"Reorganization energy (electron): {reorg_energy_electron} eV\n")
        log.write(f"Reorganization energy (hole): {reorg_energy_hole} eV\n")
        log.write(f"Electron Transfer integral J: {J_electron} eV\n")
        log.write(f"Hole Transfer integral J: {J_hole} eV\n")
    
    print(f"Reorganization energy (electron): {reorg_energy_electron} eV")
    print(f"Reorganization energy (hole): {reorg_energy_hole} eV")
    print(f"Electron Transfer integral J: {J_electron} eV")
    print(f"Hole Transfer integral J: {J_hole} eV")
