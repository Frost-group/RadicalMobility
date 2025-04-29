import subprocess
import os
from openbabel import openbabel
import math
import numpy as np

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

    output_mol_file = f"stacked_dimer_{stack_distance}.mol"
    obConversion.WriteFile(mol, output_mol_file)
    return output_mol_file

def calculate_transfer_integral(xyz_file, charge, uhf, log_file_path):
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

    if hole_transport is None or charge_transport is None:
        log.write("Error: Hole or charge transport values not found in the output.\n")
        return None, None
    else:
        return hole_transport, charge_transport
        
def calculate_rate_and_mobility(J, reorg_energy, log_file_path):
    h = 4.135667696E-15  # Planck constant in eV·s
    kT = 0.000086173 * 300  # Thermal energy in eV at 300K
    distance = 4  # Distance between dimers in nm

    rate = (((J**2) / h) * (1 / (math.sqrt(math.pi * (reorg_energy + 0.3) * kT))) *
            math.exp(-(reorg_energy + 0.3) / (4 * kT)))
    
    mobility = (((1.6E-19) * (distance * 1E-7) ** 2) * rate) / (2 * 300 * kT * 1.6E-19)
    
    with open(log_file_path, 'a') as log:
        log.write(f"Rate: {rate:.6e} s^-1\n")
        log.write(f"Mobility: {mobility:.6e} cm^2/Vs\n")
    
    return rate, mobility

if __name__ == "__main__":
    input_mol_file = "Katherine1.mol"
    radical_monomer_charge = -1
    radical_monomer_uhf = 2
    radical_dimer_charge = -2
    radical_dimer_uhf = 4
    log_file_path = "output_log.txt"
    stack_distances = np.arange(6.0, 10.5, 0.5)  # Example range from 6.0 to 10.0
    
    with open(log_file_path, 'w') as log:
        log.write("Optimization of monomers started\n")
    
    opt_radical_monomer = optimize_geometry(input_mol_file, 'radical_monomer.mol', radical_monomer_charge, radical_monomer_uhf, log_file_path)
    
    for stack_distance in stack_distances:
        dimer = stack_dimer(opt_radical_monomer, stack_distance, rotation_angle=0)
        loose_opt_dimer = optimize_geometry_gfnff(dimer, f'loose_opt_dimer_{stack_distance}.mol', radical_dimer_charge, radical_dimer_uhf, log_file_path)
        opt_dimer = optimize_geometry(loose_opt_dimer, f'opt_dimer_{stack_distance}.mol', radical_dimer_charge, radical_dimer_uhf, log_file_path)
        J_hole, J_electron = calculate_transfer_integral(opt_dimer, radical_dimer_charge, radical_dimer_uhf, log_file_path)
        
        # Compute rate and mobility
        rate_hole, mobility_hole = calculate_rate_and_mobility(J_hole, reorg_energy_hole, log_file_path)
        rate_electron, mobility_electron = calculate_rate_and_mobility(J_electron, reorg_energy_electron, log_file_path)

        
        with open(log_file_path, 'a') as log:
            log.write(f"Stack distance: {stack_distance} Å\n")
            log.write(f"Electron Transfer integral J: {J_electron} eV\n")
            log.write(f"Hole Transfer integral J: {J_hole} eV\n")
            log.write("Electron Mobility: {mobility_electron}")
            log.write("Hole Mobility: {mobility_hole}")