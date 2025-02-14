import numpy as np
from scipy.spatial.distance import pdist, squareform

def read_mol(file_path):
    atoms = []
    coordinates = []

    with open(file_path, 'r') as file:
        lines = file.readlines()
        
        atom_count = int(lines[3][0:3].strip())
        bond_count = int(lines[3][3:6].strip())
        
        for i in range(4, 4 + atom_count):
            parts = lines[i].split()
            x, y, z = map(float, parts[:3])
            atom = parts[3]
            atoms.append(atom)
            coordinates.append([x, y, z])
    
    return atoms, np.array(coordinates)

def identify_fragments(atoms, coordinates, bond_threshold=1.3):
    distance_matrix = squareform(pdist(coordinates))
    n_atoms = len(atoms)
    visited = [False] * n_atoms
    fragments = []

    def dfs(atom_index, fragment):
        visited[atom_index] = True
        fragment.append(atom_index)
        for j in range(n_atoms):
            if not visited[j] and distance_matrix[atom_index][j] < bond_threshold:
                dfs(j, fragment)

    for i in range(n_atoms):
        if not visited[i]:
            fragment = []
            dfs(i, fragment)
            fragments.append(fragment)

    return fragments

def format_fragment_indices(indices):
    ranges = []
    start = indices[0]
    end = indices[0]

    for i in range(1, len(indices)):
        if indices[i] == end + 1:
            end = indices[i]
        else:
            if start == end:
                ranges.append(f"{start}")
            else:
                ranges.append(f"{start}-{end}")
            start = indices[i]
            end = indices[i]
    
    if start == end:
        ranges.append(f"{start}")
    else:
        ranges.append(f"{start}-{end}")
    
    return ", ".join(ranges)

def print_fragments(fragments, atoms):
    for i, fragment in enumerate(fragments):
        fragment_indices = sorted([idx + 1 for idx in fragment]) 
        formatted_indices = format_fragment_indices(fragment_indices)
        fragment_atoms = [atoms[idx] for idx in fragment]
        print(f"Fragment {i + 1}: {formatted_indices} ({', '.join(fragment_atoms)})")

mol_file = 'opt_radical_stacked_dimer_blatter1.mol'  
atoms, coordinates = read_mol(mol_file)
fragments = identify_fragments(atoms, coordinates)
print_fragments(fragments, atoms)
