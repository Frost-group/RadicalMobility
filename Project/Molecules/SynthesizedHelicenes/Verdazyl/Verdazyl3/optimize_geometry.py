import subprocess
import os

def determine_charge_multiplicity(state):
    if state == "anion":
        charge = -2
        multiplicity = 3
    elif state == "cation":
        charge = 2
        multiplicity = 3
    else:
        charge = 0
        multiplicity = 5
    return charge, multiplicity

def optimize_geometry(input_xyz_file, output_xyz_file, charge, multiplicity, log_file):
    with open(log_file, 'w') as log:
        subprocess.run(['xtb', input_xyz_file, '--opt', '--gfn2', '--chrg', str(charge), '--uhf', str(multiplicity-1), '--tight'], check=True, stdout=log, stderr=log)
    subprocess.run(['cp', 'xtbopt.mol', output_xyz_file], check=True)
    return output_xyz_file

def optimize_for_states(input_xyz_file):
    states = ["radical", "anion", "cation"]
    output_files = {}
    
    for state in states:
        charge, multiplicity = determine_charge_multiplicity(state)
        output_xyz_file = f"opt_{state}_{os.path.splitext(input_xyz_file)[0]}.mol"
        log_file = f"{state}_dimer_optimization.log"
        optimize_geometry(input_xyz_file, output_xyz_file, charge, multiplicity, log_file)
        output_files[state] = output_xyz_file
        print(f"Optimized geometry for {state} saved in: {output_xyz_file}")
    
    return output_files

if __name__ == "__main__":
    input_xyz_file = "stacked_dimer_Verdazyl3.mol"
    optimized_files = optimize_for_states(input_xyz_file)

    print("\n--- Optimized Files ---")
    for state, file in optimized_files.items():
        print(f"{state.capitalize()} optimized file: {file}")
