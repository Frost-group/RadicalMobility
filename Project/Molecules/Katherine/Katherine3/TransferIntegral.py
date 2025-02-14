import subprocess

def calculate_transfer_integral(xyz_file):
    result = subprocess.run(['xtb', xyz_file, '--dipro', '--gfn', '1', '--chrg', '0', '--blatter1.UHFfrag'], 
                                check=True, capture_output=True, text=True)
    for line in result.stdout.splitlines():
        if 'total |J(AB,eff)| for charge transfer (CT)' in line:
            parts = line.split()
            for part in parts:
                try:
                    transfer_integral = float(part)
                    return transfer_integral
                except ValueError:
                    continue
    return None

if __name__ == "__main__":
    xyz_file = "OptimizedStack.xyz"
    transfer_integral = calculate_transfer_integral(xyz_file)
    print(f"Transfer Integral: {transfer_integral} eV")
