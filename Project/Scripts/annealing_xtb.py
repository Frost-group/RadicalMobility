import os
import shutil
import subprocess
import re
import sys

# ===================== CONFIGURATION =====================
n_cycles = 10   # Number of full 300 → 800 → 300K annealing cycles
initial_temp = 300   # Starting temperature
max_temp = 1000       # Peak temperature
heat_step = 10       # Increment for heating phase (in K)
cool_step = 10       # Decrement for cooling phase (in K)
base_xyz = "packed_box.xyz"   # Initial structure file
md_input_template = "md_settings.inp"  # MD settings template file
log_filename = "md_log.txt"  # Name of the log file for each temperature step

# ===================== CHECK FILES =====================
if not os.path.isfile(base_xyz):
    print(f"Error: Initial structure file '{base_xyz}' not found.")
    sys.exit(1)
if not os.path.isfile(md_input_template):
    print(f"Error: MD settings file '{md_input_template}' not found.")
    sys.exit(1)

# Read the MD settings template into memory for editing
with open(md_input_template, "r") as f:
    md_settings_lines = f.readlines()

# Generate temperature sequence for heating (up) and cooling (down)
heating_temps = list(range(initial_temp, max_temp + 1, heat_step))
cooling_temps = list(range(max_temp - cool_step, initial_temp - 1, -cool_step))

# Full temperature sequence for one annealing cycle
temperature_sequence = heating_temps + cooling_temps

# ===================== RUN MULTIPLE CYCLES =====================
prev_structure_file = None  # Track last structure file

for cycle in range(1, n_cycles + 1):
    cycle_dir = f"Cycle_{cycle}"  # Create a new directory per cycle
    os.makedirs(cycle_dir, exist_ok=True)
    
    print(f"\n==================== Starting Annealing Cycle {cycle}/{n_cycles} ====================")
    
    for i, temp in enumerate(temperature_sequence):
        # Distinguish between heating and cooling to prevent overwriting
        if i < len(heating_temps):
            step_label = "up"
        else:
            step_label = "down"

        dir_name = os.path.join(cycle_dir, f"T_{temp}_{step_label}")  # Unique directory names
        os.makedirs(dir_name, exist_ok=True)

        # Determine input structure file
        if prev_structure_file is None:
            structure_file = base_xyz  # First step: use the base XYZ file
        else:
            structure_file = prev_structure_file  # Subsequent steps use the last scoord file

        # Copy structure file into the new directory
        shutil.copy(structure_file, os.path.join(dir_name, os.path.basename(structure_file)))

        # Modify MD settings with the correct temperature
        new_md_settings = []
        for line in md_settings_lines:
            if re.match(r'\s*temp\s*=', line):  # Match only the temp line
                new_md_settings.append(f"temp = {temp}\n")  # Replace correctly
            else:
                new_md_settings.append(line)  # Keep all other lines unchanged

        # Save the updated md_settings.inp file in the new directory
        md_settings_path = os.path.join(dir_name, os.path.basename(md_input_template))
        with open(md_settings_path, "w") as f:
            f.writelines(new_md_settings)

        # Run xTB MD simulation and log output
        print(f"Running MD at {temp} K (Cycle {cycle}, Step {step_label})...")
        cmd = ["xtb", os.path.basename(structure_file), "--md", "--gfnff", "--input", os.path.basename(md_input_template)]
        
        log_path = os.path.join(dir_name, log_filename)
        with open(log_path, "w") as log_file:
            result = subprocess.run(cmd, cwd=dir_name, stdout=log_file, stderr=subprocess.STDOUT)

        # Check if xTB run was successful
        if result.returncode != 0:
            print(f"Error: xTB MD at {temp}K in Cycle {cycle} failed (exit code {result.returncode}).")
            sys.exit(1)

        # Find the highest-numbered scoord file in this directory
        scoord_files = [f for f in os.listdir(dir_name) if f.startswith("scoord.")]
        if not scoord_files:
            print(f"Error: No scoord.* files found after {temp}K MD run in Cycle {cycle}.")
            sys.exit(1)

        # Get the file with the highest numbered extension
        max_index = -1
        latest_file = None
        for fname in scoord_files:
            try:
                idx = int(fname.split('.')[-1])
            except ValueError:
                continue
            if idx > max_index:
                max_index = idx
                latest_file = fname

        if latest_file is None:
            latest_file = sorted(scoord_files)[-1]

        # Update prev_structure_file for the next step
        prev_structure_file = os.path.join(dir_name, latest_file)

print("All annealing cycles completed successfully.")
