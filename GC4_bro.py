import numpy as np
import matplotlib.pyplot as plt

# Given data (Givna data)
L = 14.45  # Längd på balken i meter
W_s = 6.1     #kN/m
W_d = 4.47    #kN/m
W_t = 12.3    #kN/m
omega = 0.52
eta = 0.31
P = 30      #kN
total_load = L*(W_d+W_s+omega*W_t) + P

R_a = (L/2) * (W_d + W_s + omega*omega*W_t)+ P* (1-eta)
R_b = total_load - R_a
print(f'Rekationskraft vid stöd A: {R_a:.2f} kN')
print(f'Rekationskraft vid stöd B: {R_b:.2f} kN')
# 1. Kraftjämvikt: R_A + R_B == total_load
kraft_jamvikt = R_a + R_b
print(f"Kraftjämvikt: R_A + R_B = {kraft_jamvikt:.2f} kN, förväntat: {total_load:.2f} kN")

Moment_A =  (L*L) * ((W_d + W_s)/2 + W_t*omega*(1-(omega/2))) + P * eta * L - R_b*L
Moment_B = L**2/2 * (W_d + W_s + omega**2*W_t) + P*L*(1-eta) - R_a*L
print(f'Resulterande moment vid punkt A: {Moment_A} kNm')
print(f'Resulterande moment vid punkt B: {Moment_B} kNm')

x_P = L * 0.31  # Punktlasten verkar vid x = L*0,31
x_Wt = L* (1-omega) #Trängsellasten börjar verka vid x = L*(1-omega)
# Definiera x-axeln (position längs balken)
x = np.linspace(0, L, 1000)  # 1000 punkter från 0 till L

# Initiera arrayer för tvärkraft (T) och böjmoment (M)
T = np.zeros_like(x)  # Tvärkraft
M = np.zeros_like(x)  # Böjmoment

# Beräkna tvärkraft (T) och böjmoment (M) längs balken
for i, xi in enumerate(x):
    if xi < eta * L:
        # Region 1: 0 <= x < 0.31L (endast W_d och W_s)
        T[i] = R_a - (W_d + W_s) * xi
        M[i] = R_a * xi - 0.5 * (W_d + W_s) * xi**2
    elif xi < 0.48 * L:
        # Region 2: 0.31L <= x < 0.48L (W_d, W_s och P)
        T[i] = R_a - (W_d + W_s) * xi - P
        M[i] = R_a * xi - 0.5 * (W_d + W_s) * xi**2 - P * (xi - eta * L)
    else:
        # Region 3: 0.48L <= x <= L (W_d, W_s, W_t och P)
        T[i] = R_a - (W_d + W_s) * xi - W_t * (xi - 0.48 * L) - P
        M[i] = R_a * xi - 0.5 * (W_d + W_s) * xi**2 - 0.5 * W_t * (xi - 0.48 * L)**2 - P * (xi - eta * L)

# Kontrollera global jämvikt (summan av krafter och moment bör vara noll)
s = R_a + R_b - (W_t+W_s+W_d) * L - P  # Summa av krafter i vertikal led
sum_moments = R_b * L - (W_t+W_s+W_d) * L * (L / 2) - P * (L / 2)  # Summa av moment kring stöd A


# Hitta maximalt böjmoment och dess position
max_M = np.max(M)  # Maximalt böjmoment
x_max_M = x[np.argmax(M)]  # Position för maximalt böjmoment

print(f"Maximalt böjmoment: {max_M:.2f} kNm vid x = {x_max_M:.2f} m")

# Rita tvärkraftdiagram (T-diagram)

plt.plot(x, T, label="Tvärkraft (T) [kN]", color="blue")
plt.axhline(0, color="black", linewidth=0.5, linestyle="--")  # Rita en horisontell linje vid 0
plt.title("Tvärkraftdiagram (T-diagram)")
plt.xlabel("Position längs balken (x) [m]")
plt.ylabel("Tvärkraft (T) [kN]")
plt.grid()  # Lägg till ett rutnät
plt.legend()  # Visa legend
plt.show()

# Rita böjmomentdiagram (M-diagram)

plt.plot(x, M, label="Böjmoment (M) [kNm]", color="red")
plt.axhline(0, color="black", linewidth=0.5, linestyle="--")  # Rita en horisontell linje vid 0
plt.title("Böjmomentdiagram (M-diagram)")
plt.xlabel("Position längs balken (x) [m]")
plt.ylabel("Böjmoment (M) [kNm]")
plt.grid()  # Lägg till ett rutnät
plt.legend()  # Visa legend
plt.show()
