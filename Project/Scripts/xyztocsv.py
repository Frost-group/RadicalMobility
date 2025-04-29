import pandas as pd
import numpy as np

BOHR_TO_ANGSTROM = 0.52917721067

# Define atomic masses (simplified common elements, fallback to 12.0)
mass_table = {
    'H': 1.008, 'C': 12.01, 'N': 14.01, 'O': 16.00,
    'F': 18.998, 'S': 32.06, 'Cl': 35.45, 'Br': 79.90, 'I': 126.90,
    'B': 10.81
}

def load_xyz_and_compute_COM(xyz_file, n_molecules=100, atoms_per_molecule=58):
    with open(xyz_file, 'r') as f:
        lines = f.readlines()[2:]  # Skip first two header lines

    # Extract atom types and positions, converted from Bohr to Å
    data = []
    for line in lines:
        tokens = line.strip().split()
        if len(tokens) < 4:
            continue
        atom = tokens[0].capitalize()
        x, y, z = map(float, tokens[1:4])
        x *= BOHR_TO_ANGSTROM
        y *= BOHR_TO_ANGSTROM
        z *= BOHR_TO_ANGSTROM
        data.append((atom, x, y, z))

    assert len(data) == n_molecules * atoms_per_molecule, "Mismatch between molecule count and total atom count."

    coms = []

    for i in range(n_molecules):
        start = i * atoms_per_molecule
        end = start + atoms_per_molecule
        mol_atoms = data[start:end]

        coords = np.array([[x, y, z] for _, x, y, z in mol_atoms])
        masses = np.array([mass_table.get(atom, 12.0) for atom, _, _, _ in mol_atoms])
        total_mass = np.sum(masses)
        com = np.sum(coords * masses[:, None], axis=0) / total_mass

        coms.append({'Mol': i, 'x': com[0], 'y': com[1], 'z': com[2]})

    return pd.DataFrame(coms)

# Run the conversion
xyz_file = 'scoord.19.xyz'
output_csv = 'scoord.19_com.csv'
com_df = load_xyz_and_compute_COM(xyz_file, n_molecules=100, atoms_per_molecule=58)
com_df.to_csv(output_csv, index=False)

print(f"✅ Saved molecular COM positions to: {output_csv}")
