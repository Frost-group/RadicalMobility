import os
import subprocess
import csv
import numpy as np

def read_xyz_box(filename, atoms_per_mol=58):
    with open(filename, 'r') as f:
        lines = f.readlines()

    total_atoms = int(lines[0])
    all_atoms = lines[2:]
    n_mols = total_atoms // atoms_per_mol
    molecules = []

    for i in range(n_mols):
        start = i * atoms_per_mol
        end = (i + 1) * atoms_per_mol
        mol_lines = all_atoms[start:end]
        molecules.append(mol_lines)

    return molecules

def write_dimer_xyz(mol1, mol2, filename):
    with open(filename, 'w') as f:
        f.write(f"{len(mol1) + len(mol2)}\n")
        f.write("Dimer\n")
        f.writelines(mol1)
        f.writelines(mol2)

def parse_xyz_coords(lines):
    coords = []
    for line in lines:
        parts = line.strip().split()
        if len(parts) >= 4:
            coords.append([float(x) for x in parts[1:4]])
    return np.array(coords)

def calculate_com_distance(mol1, mol2):
    coords1 = parse_xyz_coords(mol1)
    coords2 = parse_xyz_coords(mol2)
    com1 = np.mean(coords1, axis=0)
    com2 = np.mean(coords2, axis=0)
    return np.linalg.norm(com2 - com1)

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
                    hole_transport = float(parts[-3])
                except ValueError:
                    log.write(f"Error: Could not parse hole transport from line: {line}\n")
            elif 'charge transport (unocc. MOs)' in line:
                parts = line.split()
                try:
                    charge_transport = float(parts[-3])
                except ValueError:
                    log.write(f"Error: Could not parse charge transport from line: {line}\n")

        if hole_transport is None or charge_transport is None:
            log.write(f"Error: Hole or charge transport values not found for {xyz_file}.\n")

        return hole_transport, charge_transport

def run_dipro_on_all_dimers(box_xyz, atoms_per_mol=58, charge=0, uhf=1):
    molecules = read_xyz_box(box_xyz, atoms_per_mol)
    output_dir = "dipro_dimers"
    os.makedirs(output_dir, exist_ok=True)
    log_file = os.path.join(output_dir, "dipro_log.txt")
    result_file = os.path.join(output_dir, "transfer_integrals.csv")

    raw_results = []

    for i in range(len(molecules)):
        for j in range(i + 1, len(molecules)):
            dimer_name = f"dimer_{i}_{j}"
            dimer_xyz = os.path.join(output_dir, f"{dimer_name}.xyz")

            write_dimer_xyz(molecules[i], molecules[j], dimer_xyz)
            distance = calculate_com_distance(molecules[i], molecules[j])

            try:
                hole_J, charge_J = calculate_transfer_integral(
                    dimer_xyz, charge=charge, uhf=uhf, log_file_path=log_file
                )

                # Skip if either transfer integral is 0
                if hole_J == 0.0 or charge_J == 0.0:
                    continue

                raw_results.append([i, j, hole_J, charge_J, distance])

            except subprocess.CalledProcessError as e:
                with open(log_file, 'a') as log:
                    log.write(f"Failed to run DIPRO on {dimer_xyz}\n{e.stderr}\n")

    # Write cleaned results
    with open(result_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['MolA', 'MolB', 'Hole_Transfer', 'Charge_Transfer', 'Distance'])
        writer.writerows(raw_results)

    print(f"Finished all dimers. Final results saved to: {result_file}")

# === Run it ===
run_dipro_on_all_dimers("scoord.19.xyz", atoms_per_mol=58, charge=0, uhf=2)
