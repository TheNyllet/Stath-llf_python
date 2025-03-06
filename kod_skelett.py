import numpy as np
import matplotlib.pyplot as plt
from utils import *

# Indata
E = 2.1e11  # Pa
A0 = 7.85e-3  # m^2
P = 1.5e5  # N
L = 2  # m
sigma_s = 2.3e8 # Pa

# Topology
Edof = np.array([
    [3, 4, 7, 8],
    [3, 4, 5, 6],
    [1, 2, 7, 8],
    [1, 2, 5, 6],
    [5, 6, 7, 8],
    [5, 6, 9, 10],
    [7, 8, 9, 10],
    [7, 8, 11, 12],
    [9, 10, 11, 12],
    [9, 10, 13, 14],
    [11, 12, 13, 14]
])

# Koordinater för varje nod i numerisk ordning:
Coord = L * np.array([
    [0.0, 0.0],
    [1.0, 0.0],
    [0.0, 2.0],
    [1.0, 2.0],
    [1.0, 3.0],
    [3.0, 2.0],
    [3.0, 3.0]
])

# x-koordinater för varje element i numerisk ordning
Ex = L*np.array([
    [1.0,1.0], 
    [1.0,0.0], 
    [0.0,1.0], 
    [0.0,0.0], 
    [0.0,1.0], 
    [0.0,1.0], 
    [1.0,1.0], 
    [1.0,3.0], 
    [1.0,3.0], 
    [1.0,3.0], 
    [3.0,3.0]
])

# y-koordinater för varje element i numerisk ordning
Ey = L*np.array([
    [0.0,2.0], 
    [0.0,2.0], 
    [0.0,2.0], 
    [0.0,2.0], 
    [2.0,2.0], 
    [2.0,3.0], 
    [2.0,3.0], 
    [2.0,2.0], 
    [3.0,2.0], 
    [3.0,3.0], 
    [2.0,3.0]
])

nel = len(Ex)  # Antal element
ndofs = 2 * len(Coord)  # Totalt antal frihetsgrader

# Plot mesh
eldraw2(Ex, Ey)

# Fördefinera styvhetsmatrisen och kraftvektorn
K = np.zeros((ndofs, ndofs))
f = np.zeros(ndofs)

# Assemblera element
for el in range(nel):
    # Hämta elementets koordinater och längd
    x1, x2 = Ex[el]
    y1, y2 = Ey[el]
    L_el = np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)

    # Beräkna vinkeln theta och cos/sin
    theta = np.arctan2(y2 - y1, x2 - x1)
    c = np.cos(theta)
    s = np.sin(theta)

    # Bestäm tvärsnittsarean A beroende på vilket element det är
    if el + 1 in [4, 6, 9]:
        A = A0 * 2
    else:
        A = A0

    # Räkna ut styvhetsmatrisen (Se föreläsningar eller Ekvation 11.18 i kursboken [Hållfasthetslära, Allmänna tillstånd])
    Ke = E * A / L_el * np.array([
        [c ** 2, c * s, -c ** 2, -c * s],
        [c * s, s ** 2, -c * s, -s ** 2],
        [-c ** 2, -c * s, c ** 2, c * s],
        [-c * s, -s ** 2, c * s, s ** 2]
    ])

    # Assemblera in element styvhetsmatrisen och globala matrisen
    K = assem(Edof[el], K, Ke)

# Lägg till kraften P i lastvektorn:
f[11] = -P

bcdofs = np.array([i for i in range(1,5)])  # Definera de låsta frihetsgraderna
bcvals = np.zeros(len(bcdofs))  # Se till att frihetsgraderna är låsta

# Lös ekvationssystemet (: använd solveq i utils.py)
a = solveq(K, f, bcdofs, bcvals)[0]

# Extrahera elementförskjutningar
Ed = extract_eldisp(Edof, a)

# Definera variabler för uppg. 3-5
max_drag = 0
max_tryck = 0
min_sigma = 0
max_sigma = 0

# Räkna ut krafter och spänningar i varje element
for el in range(nel):
    # Ändra tvärsnittsarean på vissa element
    if el + 1 in [4, 6, 9]:
        A = A0 * 2
        width = 2
    else:
        A = A0
        width = 1

    # Materialegenskaper
    ep = [E, A]
    
    # Beräkna snittkrafter
    N = bar2s(Ex[el], Ey[el], ep, Ed[el])[0]

    # Beräkna spänningen manuellt (sigma = N / A)
    sigma = N / A

    # Spara beloppet av det högsta och lägsta sigma samt vilket element som har det
    if sigma > max_sigma: 
        max_sigma = sigma
        max_tryck = el+1
    elif sigma < min_sigma: 
        min_sigma = sigma
        max_drag = el+1

    # Runda och omvandla enheter
    N = round(N/1e3,2) # kN
    sigma = round(sigma/1e6,2) # MPa 

    print(f"Element {el+1}: Kraft = {N} KN, Spänning = {sigma} MPa")

    # Färga elementet beroende på om det är drag eller tryck
    if sigma > 0: color = 'r'
    elif sigma < 0: color = 'b'
    else: color = 'k'

    # Plotta elementen
    ex =  Ex[el,:] + Ed[el,[0,2]]
    ey =  Ey[el,:] + Ed[el,[1,3]]
    plt.plot(ex, ey, color=color, linewidth=width)

# Skriv ut resultat
print(f'\nStången med högst dragspänning är stång nr. {max_drag}')
print(f'Stången med högst tryckspänning är stång nr. {max_tryck}')

# Beräkna maximala spänningen
max_sigma = max(abs(max_sigma),abs(min_sigma))

# Beräkna maximala kraften
max_P = abs(P*sigma_s/max_sigma)
print(f'Fackverket deformerar plastiskt vid P = {round(max_P/1e3)} kN')

# Uppgift 5
min_A0 = A0*P/max_P*1e4
print(f'Med P = 150 kN börjar fackverket deformera vid A0 = {round(min_A0,1)} cm^2')

plt.show()