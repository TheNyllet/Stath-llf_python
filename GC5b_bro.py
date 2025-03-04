import numpy as np
import matplotlib.pyplot as plt

# Data
M_max = 510.5e6  # Maximalt böjmoment [Nmm], från GC4c
I_y = 712506666.7  # Yttröghetsmoment [mm^4], från GC5a
h_total = 560  # Total höjd av I-tvärsnittet [mm] (520 mm liv + 2 * 20 mm flänsar)
z = np.linspace(-h_total/2, h_total/2, 100)  # z-koordinater [mm]

# Normalspänning
sigma_z = (M_max * z) / I_y  # Normalspänning [MPa]

# Plotta normalspänning
plt.plot(sigma_z, z, label="Normalspänning")
plt.xlabel("Normalspänning (MPa)")
plt.ylabel("z-koordinat (mm)")
plt.title("Normalspänning över tvärsnittshöjd")
plt.grid()
plt.legend()
plt.show()

# Maximal drag- och tryckspänning
sigma_max_tension = np.max(sigma_z)  # Maximal dragspänning
sigma_max_compression = np.min(sigma_z)  # Maximal tryckspänning

print(f"Maximal dragspänning: {sigma_max_tension:.2f} MPa")
print(f"Maximal tryckspänning: {sigma_max_compression:.2f} MPa")