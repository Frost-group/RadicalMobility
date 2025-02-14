import numpy as np

def calculate_charge_carrier_mobility(reorganization_energy, transfer_integral, temperature):
    e, k_B, hbar = 1.602e-19, 8.6173e-5, 6.582e-16
    prefactor = (transfer_integral**2 / hbar**2) * (np.pi * hbar / (reorganization_energy + 0.3))
    exponential_term = np.exp(-reorganization_energy / (4 * k_B * temperature))
    D = prefactor * np.sqrt(np.pi / (4 * reorganization_energy * k_B * temperature)) * exponential_term
    return (e / (k_B * temperature)) * D * 1e4

if __name__ == "__main__":
    reorganization_energy = 0.2 
    transfer_integral = 0.05     
    mobility = calculate_charge_carrier_mobility(reorganization_energy, transfer_integral, 300)
    print(f"Charge Carrier Mobility: {mobility:.6f} cm²/V·s")
