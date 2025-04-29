import csv

def calculate_diffusion_coefficient(rates, displacements):
    """
    Calculates the overall diffusion coefficient using the formula:
    D_overall = (1/6) * (kx * ax^2 + ky * ay^2 + kz * az^2)
    where kx, ky, kz are transfer rates and ax, ay, az are displacements.
    For fewer than 3 variables, the formula adapts dynamically.
    """
    num_variables = len(rates)
    coefficient = 1 / (2 * num_variables)  # Adapts for 1, 2, or 3 variables
    D_overall = coefficient * sum(rate * disp**2 for rate, disp in zip(rates, displacements))
    return D_overall

def calculate_mobility_from_diffusion(diffusion_coefficient, temperature=300):
    """
    Calculates charge carrier mobility using the Einstein relation:
    μ = D / (k_B * T).
    """
    k_B = 0.000086173  # Boltzmann constant (eV/K)
    mobility = diffusion_coefficient / (k_B * temperature)  # μ = D / (k_B * T)
    return mobility

def find_entry_and_calculate_mobility(database_file, molecule_name, packing_type, displacements, rotations):
    """
    Searches the database for the matching molecule, displacement, and rotation.
    Calculates and returns mobility based on the packing type.
    """
    electron_rates = []
    hole_rates = []
    
    with open(database_file, 'r') as file:
        reader = csv.DictReader(file)
        for i in range(len(displacements)):
            dx, dy, dz = displacements[i]
            rx, ry, rz = rotations[i]
            
            for row in reader:
                if (row['Molecule'] == molecule_name and
                    float(row['dx']) == dx and float(row['dy']) == dy and float(row['dz']) == dz and
                    float(row['rx']) == rx and float(row['ry']) == ry and float(row['rz']) == rz):
                    
                    electron_rate = float(row['Electron Rate (s^-1)'])
                    hole_rate = float(row['Hole Rate (s^-1)'])
                    electron_rates.append(electron_rate)
                    hole_rates.append(hole_rate)
                    break  # Exit loop once matching entry is found
    
    if electron_rates and hole_rates:
        displacements_magnitudes = [sum(d**2 for d in disp)**0.5 for disp in displacements]
        D_electron = calculate_diffusion_coefficient(electron_rates, displacements_magnitudes)
        D_hole = calculate_diffusion_coefficient(hole_rates, displacements_magnitudes)
        
        mobility_electron = calculate_mobility_from_diffusion(D_electron)
        mobility_hole = calculate_mobility_from_diffusion(D_hole)
        
        return mobility_electron, mobility_hole
    
    return None, None  # Return None if no match found

# Example Usage
database_path = "molecular_data_log.csv"
molecule_name = "Base3tBu.mol"
packing_type = "herringbone"  # Options: "slipped stack", "slipped pi stack", "brick wall", "herringbone"

displacements = [(2, 2, 0), (3, 3, 0), (4, 4, 0)] if packing_type == "herringbone" else \
               [(2, 2, 0), (3, 3, 0)] if packing_type == "brick wall" else \
               [(2, 2, 0)]  # Default for slipped stack

rotations = [(30, 30, 0), (45, 45, 0), (60, 60, 0)] if packing_type == "herringbone" else \
            [(30, 30, 0), (45, 45, 0)] if packing_type == "brick wall" else \
            [(30, 30, 0)]  # Default for slipped stack

mobility_electron, mobility_hole = find_entry_and_calculate_mobility(database_path, molecule_name, packing_type, displacements, rotations)

if mobility_electron is not None and mobility_hole is not None:
    print(f"Electron Mobility: {mobility_electron} cm^2/Vs")
    print(f"Hole Mobility: {mobility_hole} cm^2/Vs")
else:
    print("No matching entry found in the database.")