import os
import shutil
import subprocess
import re
import sys
import numpy as np

# ===================== CONFIGURATION =====================
n_cycles = 1   # Number of full annealing cycles
initial_temp = 300   # Minimum temperature (final state)
max_temp = 800       # Peak temperature
heat_step = 50       # Increment for heating phase (in K)
cooling_lambda = 0.1  # Exponential decay rate (adjust for desired speed)
base_xyz = "opt_dimer.mol"   # Initial structure file
md_input_template = "md_settings.inp"  # MD settings template file
log_filename = "md_log.txt"  # Name of the log file for each temperature step
energy_tolerance = 1e-5  # Convergence threshold for energy (Eh)
k_B = 3.166811563e-6  # Boltzmann constant in Eh/K

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

# **Generate Temperature Sequences**
heating_temps = list(range(initial_temp, max_temp + 1, heat_step))  # Heating up in fixed steps

# **Exponential Cooling Sequence**
cooling_temps = []
n = 0
while True:
    temp = initial_temp + (max_temp - initial_temp) * np.exp(-cooling_lambda * n)  # Exponential decay
    temp = int(round(temp))  # Convert to integer
    if temp <= initial_temp:
        break  # Stop when reaching 300K
    cooling_temps.append(temp)
    n += 1

cooling_temps = sorted(set(cooling_temps), reverse=True)  # Ensure unique & descending temperatures

# **Full Annealing Schedule**
temperature_sequence = heating_temps + cooling_temps

# ===================== RUN MULTIPLE CYCLES =====================
prev_structure_file = None  # Track last structure file

for cycle in range(1, n_cycles + 1):
    cycle_dir = f"Cycle_{cycle}"  # Create a new directory per cycle
    os.makedirs(cycle_dir, exist_ok=True)
    
    print(f"\n==================== Starting Annealing Cycle {cycle}/{n_cycles} ====================")

    prev_energy = None  # Track last accepted energy

    i = 0  # Temperature index
    while i < len(temperature_sequence):
        temp = temperature_sequence[i]
        step_label = "up" if temp in heating_temps else "down"
        dir_name = os.path.join(cycle_dir, f"T_{temp}_{step_label}")  # Unique directory names
        os.makedirs(dir_name, exist_ok=True)

        # Determine input structure file
        structure_file = base_xyz if prev_structure_file is None else prev_structure_file

        # Copy structure file into the new directory
        destination_path = os.path.join(dir_name, os.path.basename(structure_file))

        # Only copy if the source and destination paths are different
        if os.path.abspath(structure_file) != os.path.abspath(destination_path):
            shutil.copy(structure_file, destination_path)

        # Modify MD settings with the correct temperature
        new_md_settings = []
        for line in md_settings_lines:
            if re.match(r'\s*temp\s*=', line):
                new_md_settings.append(f"temp = {temp}\n")
            else:
                new_md_settings.append(line)

        md_settings_path = os.path.join(dir_name, os.path.basename(md_input_template))
        with open(md_settings_path, "w") as f:
            f.writelines(new_md_settings)

        # ===================== RUN MD SIMULATION =====================
        print(f"Running MD at {temp} K (Cycle {cycle}, Step {step_label})...")

        cmd = ["xtb", os.path.basename(structure_file), "--md", "--gfnff", "--input", os.path.basename(md_input_template)]
        log_path = os.path.join(dir_name, log_filename)
        with open(log_path, "w") as log_file:
            result = subprocess.run(cmd, cwd=dir_name, stdout=log_file, stderr=subprocess.STDOUT)

        if result.returncode != 0:
            print(f"Error: xTB MD at {temp}K in Cycle {cycle} failed.")
            sys.exit(1)

        # Extract total energy (Etot) from md_log.txt
        total_energy = None
        with open(log_path, "r") as f:
            for line in f:
                if "TOTAL ENERGY" in line:
                    match = re.search(r"TOTAL ENERGY\s+([-\d.]+)\s+Eh", line)
                    if match:
                        total_energy = float(match.group(1))
                        break

        if total_energy is None:
            print(f"Warning: Could not extract Etot from {log_filename} in {dir_name}")

        # Find the highest-numbered scoord file
        scoord_files = [f for f in os.listdir(dir_name) if f.startswith("scoord.")]
        latest_file = sorted(scoord_files, key=lambda x: int(x.split('.')[-1]))[-1] if scoord_files else None

        if latest_file is None:
            print(f"Error: No scoord.* files found after {temp}K MD run in Cycle {cycle}.")
            sys.exit(1)

        # Ensure we use the latest scoord file, even for rejected steps
        latest_structure_file = os.path.join(dir_name, latest_file)

        # ===================== METROPOLIS CRITERION (ONLY FOR COOLING) =====================
        accept_move = True  # Default accept for heating steps

        if step_label == "down":  # Apply Metropolis only for cooling steps
            if prev_energy is not None and total_energy > prev_energy:
                delta_E = total_energy - prev_energy
                metropolis_prob = np.exp(-delta_E / (k_B * temp))
                accept_move = np.random.rand() < metropolis_prob  # Accept with probability e^(-ΔE/kT)
                
                if not accept_move:
                    print(f"Cooling step rejected at {temp}K (ΔE = {delta_E:.6f} Eh, P = {metropolis_prob:.4f}). Retrying with latest structure.")
                    prev_structure_file = latest_structure_file  # Keep rejected step's own geometry
                    continue  # Retry at the same temperature with new geometry

        if accept_move:
            prev_energy = total_energy  # Update last accepted energy
            prev_structure_file = latest_structure_file  # Update structure for next step
            with open(os.path.join(cycle_dir, "energy_log.txt"), "a") as energy_log:
                energy_log.write(f"{temp}K - Etot: {total_energy} Eh\n")

        i += 1  # Move to the next temperature **only if the step is accepted**
    
print("All annealing cycles completed successfully.")
