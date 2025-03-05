# Program som beräknar maximala tvärkraften och böjmomentet i en balk
# Grupp nummer 87
# Författare: Lucas Molander & Hugo Nylander

import numpy as np
import matplotlib.pyplot as plt

# Given data 
L = 14.45  # balkens längd
W_s = 6.1     #kN/m (från uppgift 1)
W_d = 4.47    #kN/m(från uppgift 1)
W_t = 12.3    #kN/m(från uppgift 1)
omega = 0.52
eta = 0.31
P = 30      #kN
total_load = L*(W_d+W_s+omega*W_t) + P

R_a = (L/2) * (W_d + W_s + omega*omega*W_t)+ P* (1-eta) # Beräkning av R_a enligt ekvation från uppgift 3a
R_b = total_load - R_a # ENligt ekvation från uppgift 3a
print(f'Reakationskraft vid stöd A: {R_a:.2f} kN')
print(f'Reaktionskraft vid stöd B: {R_b:.2f} kN')
# 1. Kraftjämvikt: R_A + R_B == total_load
kraft_jamvikt = R_a + R_b - total_load
print(f"Kraftjämvikt: R_A + R_B - totala lasten = {kraft_jamvikt:.2f} kN (bör vara 0)")

Moment_A =  (L*L) * ((W_d + W_s)/2 + W_t*omega*(1-(omega/2))) + P * eta * L - R_b*L #beräkning av momenten vid a och b för att säkerställa jämvikt
Moment_B = L**2/2 * (W_d + W_s + omega**2*W_t) + P*L*(1-eta) - R_a*L
print(f'Resulterande moment vid punkt A: {Moment_A} kNm')
print(f'Resulterande moment vid punkt B: {Moment_B} kNm')

x_P = L * eta  # Punktlasten verkar vid x = eta*L
x_Wt = L* (1-omega) #Trängsellasten börjar verka vid x = L*(1-omega)
x = np.linspace(0, L, 1000)  # 1000 punkter jämnt fördelade mellan 0 och L

T = np.zeros_like(x)  # Array för tvärkrafterna 
M = np.zeros_like(x)  # Array för böjmommentem

for i, xi in enumerate(x): #Beräknar tvärkraft och böjmoment längs balken
    if xi < eta * L:
        # Region 1: 0 <= x < eta*L (innan punktlasten)
        T[i] = (W_d + W_s) * xi - R_a
        M[i] = 0.5 * (W_d + W_s) * xi**2 - R_a * xi
    elif xi < 0.48 * L:
        # Region 2: eta*L <= x < (1-omega)*L (innan trängsel lasten)
        T[i] = (W_d + W_s) * xi - R_a + P
        M[i] = 0.5 * (W_d + W_s) * xi**2 - R_a * xi + P * (xi - eta * L)
    else:
        # Region 3: Efter trängsellasten till slutet av balken
        T[i] = (W_d + W_s) * xi + W_t * (xi - 0.48 * L) + P - R_a
        M[i] = 0.5 * (W_d + W_s) * xi**2 + 0.5 * W_t * (xi - 0.48 * L)**2 - R_a * xi + P * (xi - eta * L)

#Hittar maximala tvärkraften och böjmomenten
tmax = np.max(abs(T))
print(f'Maximala tvärkraften är {tmax:.2f} kN ')
max_M = np.max(abs(M))  # Maximalt böjmoment
x_max_M = x[np.argmax(abs(M))]  # Position för maximalt böjmoment
print(f"Maximalt böjmoment: {max_M:.2f} kNm vid x = {x_max_M:.2f} m")

#plottar figurerna enligt beskrivning i uppgiften
plt.plot(x, T)
plt.axhline(0, color="black", linewidth=0.5, linestyle="--") 
plt.title("Tvärkraftdiagram (T-diagram)")
plt.xlabel("Position längs balken x [m]")
plt.ylabel("Tvärkraft T(x) [kN]")
plt.grid()
plt.show()

plt.plot(x, M)
plt.axhline(0, color="black", linewidth=0.5, linestyle="--")
plt.title("Böjmomentdiagram (M-diagram)")
plt.xlabel("Position längs balken [m]")
plt.ylabel("Böjmoment [kNm]")
plt.grid()
plt.show()