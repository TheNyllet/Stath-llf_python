# Hugo Nylander & Lucas Molander

import numpy as np
import matplotlib.pyplot as plt

# Data
M_max = 518.9e6 # Maximalt böjmoment [Nmm], från GC4c
I_y = 758986666.7  # Yttröghetsmoment [mm^4], från GC5a
h_total = 560  # Total höjd av I-tvärsnittet [mm] (520 mm liv + 2 * 20 mm flänsar)
d = np.linspace(-h_total/2, h_total/2, 100)  # Avstånd från mittpunkt [mm]
z = np.linspace(0, h_total, 100)  # z-koordinat [mm]

# Normalspänning
sigma_z = (M_max * d) / I_y  # Normalspänning [MPa]

# Maximal drag- och tryckspänning
sigma_max_tension = np.max(sigma_z)  # Maximal dragspänning
sigma_max_compression = np.min(sigma_z)  # Maximal tryckspänning

print(f"Maximal dragspänning: {sigma_max_tension:.2f} MPa vid z = {z[np.argmax(sigma_z)]:.2f} mm")
print(f"Maximal tryckspänning: {sigma_max_compression:.2f} MPa vid z = {z[np.argmin(sigma_z)]:.2f} mm")

# Plotta normalspänning
plt.plot(z, sigma_z, label="Normalspänning")
plt.xlabel("z-koordinat (mm)")
plt.ylabel("Normalspänning (MPa)")
plt.title("Normalspänning över tvärsnittshöjd")
plt.grid()
plt.legend()
plt.show()