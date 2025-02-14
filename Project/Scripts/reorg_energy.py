import subprocess

def calculate_single_molecule_energy(optimized_xyz_file, charge, multiplicity, log_file):
    with open(log_file, 'w') as log:
        result = subprocess.run(
            ['xtb', optimized_xyz_file, '--gfn2', '--sp', '--chrg', str(charge), '--uhf', str(multiplicity - 1)], 
            check=True, stdout=log, stderr=log, text=True
        )
    with open(log_file, 'r') as log:
        energy_line = next(line for line in log if 'TOTAL ENERGY' in line)
    return float(energy_line.split()[3])

def calculate_reorganization_energy(neutral_xyz_file, anion_xyz_file, cation_xyz_file):
    conversion_factor = 27.2114
    E_neutral_optimized = calculate_single_molecule_energy(neutral_xyz_file, 0, 2, "neutral_energy.log")
    E_anion_optimized = calculate_single_molecule_energy(anion_xyz_file, -1, 1, "anion_energy.log")
    E_neutral_in_anion_geom = calculate_single_molecule_energy(anion_xyz_file, 0, 1, "neutral_in_anion_energy.log")
    E_anion_in_neutral_geom = calculate_single_molecule_energy(neutral_xyz_file, -1, 2, "anion_in_neutral_energy.log")
    reorg_energy_electron_hartree = (E_neutral_in_anion_geom - E_neutral_optimized) + (E_anion_in_neutral_geom - E_anion_optimized)
    reorg_energy_electron_eV = reorg_energy_electron_hartree * conversion_factor
    E_cation_optimized = calculate_single_molecule_energy(cation_xyz_file, 1, 1, "cation_energy.log")
    E_neutral_in_cation_geom = calculate_single_molecule_energy(cation_xyz_file, 0, 1, "neutral_in_cation_energy.log")
    E_cation_in_neutral_geom = calculate_single_molecule_energy(neutral_xyz_file, 1, 2, "cation_in_neutral_energy.log")
    reorg_energy_hole_hartree = (E_neutral_in_cation_geom - E_neutral_optimized) + (E_cation_in_neutral_geom - E_cation_optimized)
    reorg_energy_hole_eV = reorg_energy_hole_hartree * conversion_factor
    return reorg_energy_electron_eV, reorg_energy_hole_eV

if __name__ == "__main__":
    neutral_xyz_file = "opt_radical_Base1CO2Me.mol"
    anion_xyz_file = "opt_anion_Base1CO2Me.mol"
    cation_xyz_file = "opt_cation_Base1CO2Me.mol"
    reorg_energy_electron, reorg_energy_hole = calculate_reorganization_energy(neutral_xyz_file, anion_xyz_file, cation_xyz_file)
    print(f"Electron Reorganization Energy: {reorg_energy_electron:.6f} eV")
    print(f"Hole Reorganization Energy: {reorg_energy_hole:.6f} eV")
