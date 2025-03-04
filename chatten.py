import numpy as np
import matplotlib.pyplot as plt

# Givna värden
W_b = 4.47  # kN/m, egenvikt
W_s = 6.1   # kN/m, snölast
P = 30      # kN, punktlast
W_t = 12.3  # kN/m, tillkommande linjelast
L = 14.45   # m, balkens längd
x_p = 4.4795  # m, punktlastens position
x_t = 6.936   # m, start för tillkommande linjelast

# Totala krafter
F_total = (W_b * L) + (W_s * L) + P + (W_t * (L - x_t))

# Momentberäkning om vänstra stödet för att räkna ut R_1
M_1 = W_b * L * (L / 2)  # Moment från egenvikt
M_2 = W_s * L * (L / 2)  # Moment från snölast
M_3 = P * x_p           # Moment från punktlast
M_4 = W_t * (L - x_t) * (L + x_t) / 2  # Moment från tillkommande linjelast

# Beräkna R_1 (stödkraft vid vänstra stödet)
R_1 = (M_2 + M_3 + M_4) / L

# Beräkna R_2 (stödkraft vid högra stödet)
R_2 = F_total - R_1

# Skriv ut R_1 och R_2
print(f"R_1 (stödkraft vid vänstra stödet): {R_1:.2f} kN")
print(f"R_2 (stödkraft vid högra stödet): {R_2:.2f} kN")

# Funktion för tvärkraften
def tvarkraft(x):
    if x < x_p:
        V = R_1 - (W_b + W_s) * x
    elif x < x_t:
        V = R_1 - (W_b + W_s) * x - P
    else:
        V = R_1 - (W_b + W_s) * x - P - W_t * (x - x_t)
    
    return V

# Funktion för momentet
def moment(x):
    # Numerisk integration av tvärkraften
    integral = 0
    step = 0.01  # Stegstorlek för numerisk integration
    for i in np.arange(0, x, step):
        integral += tvarkraft(i) * step
    return integral

# Skapa ett intervall av x-värden
x_values = np.linspace(0, 14.45, 500)

# Beräkna tvärkraften och momentet för varje x-värde
V_values = [tvarkraft(x) for x in x_values]
M_values = [moment(x) for x in x_values]

# Skapa grafen
fig, ax1 = plt.subplots()

# Plot tvärkraften på vänstra y-axeln
ax1.plot(x_values, V_values, label="Tvärkraft V(x)", color='b')
ax1.set_xlabel("Position längs balken (m)")
ax1.set_ylabel("Tvärkraft (kN)", color='b')
ax1.tick_params(axis='y', labelcolor='b')

# Skapa en andra y-axel för momentet
ax2 = ax1.twinx()
ax2.plot(x_values, M_values, label="Moment M(x)", color='r')
ax2.set_ylabel("Moment (kNm)", color='r')
ax2.tick_params(axis='y', labelcolor='r')

# Lägg till grid och legend
ax1.grid(True)
fig.tight_layout()
plt.title("Tvärkraft och Moment längs balken")
plt.show()