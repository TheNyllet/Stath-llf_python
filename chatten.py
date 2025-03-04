import numpy as np
import matplotlib.pyplot as plt
from utils import *

## Indata (Geometry, material)
E = 2.1e11  # N/m^2
A0 = 7.85e-3  # m^2
P = 1.5e5  # N
L = 2  # m
sigma_s = 2.3e8 #Pa

## Topology
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
], dtype=int)

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

# x-koordinater för varje element
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

# y-koordinater för varje element
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

# Hjälpvariabler:
nel = len(Ex)  # Antal element
ndofs = 2 * len(Coord)  # Totalt antal frihetsgrader

# Plot mesh (tips: använd eldraw2 i utils.py)
eldraw2(Ex, Ey)

# Fördefinera styvhetsmatrisen och kraftvektorn
K = np.zeros((ndofs, ndofs))
f = np.zeros(ndofs)

# Assemblera element
for el in range(nel):
    x1, x2 = Ex[el]
    y1, y2 = Ey[el]
    L_el = np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)

    theta = np.arctan2(y2 - y1, x2 - x1)
    c = np.cos(theta)
    s = np.sin(theta)

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

# Bestäm bcdofs och bcvals
bcdofs = np.array([1, 2, 3, 4])  # Exempel på frihetsgrader med randvillkor
bcvals = np.array([0.0, 0.0, 0.0, 0.0])  # Värden för randvillkoren

# Lös ekvationssystemet (: använd solveq i utils.py)
a, r = solveq(K, f, bcdofs, bcvals)

# Plotta deformerad mesh (: använd eldisp2 i utils.py)
Ed = extract_eldisp(Edof, a)  # Extrahera elementförskjutningar

# Definera variabler för uppg. 3-5
max_drag = 0
max_tryck = 0
min_sigma = 0
max_sigma = 0

# Räkna ut krafter och spänningar i varje element
for el in range(nel):
    ep = [E, A0 if el + 1 not in [4, 6, 9] else A0 * 2]  # Materialegenskaper
    ed = extract_eldisp(Edof[el], a)  # Extrahera elementförskjutningar

    # Se till att ed är en 1D-array med längd 4
    if ed.ndim == 2:  # Om ed är 2D, ta första raden
        ed = ed.flatten()  # Förenkla till 1D-array
    if len(ed) != 4:  # Om längden inte är 4, fyll med nollor eller justera
        ed = np.zeros(4)  # Standardvärde om något går fel

    # Anropa bar2s
    es = bar2s(Ex[el], Ey[el], ep, ed)  # Beräkna snittkrafter

    # Om bar2s returnerar en array, ta det första värdet
    if isinstance(es, np.ndarray):
        N = es[0]  # Ta första värdet som snittkraft
    else:
        N = es  # Om det redan är ett enskilt värde

    if el + 1 in [4, 6, 9]:
        A = A0 * 2
    else:
        A = A0

    # Beräkna spänningen manuellt (sigma = N / A)
    sigma = N / A  # ep[1] är tvärsnittsarean A

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

    if sigma > 0: color = 'r'
    elif sigma < 0: color = 'b'
    else: color = 'k'
    ex =  Ex[el,:] + Ed[el,[0,2]]
    ey =  Ey[el,:] + Ed[el,[1,3]]
    plt.plot(ex, ey, color=color)

print(f'\nStången med högst dragspänning är stång nr. {max_drag}')
print(f'Stången med högst tryckspänning är stång nr. {max_tryck}')

max_sigma = max(abs(max_sigma),abs(min_sigma))

max_P = round(abs(P*sigma_s/max_sigma)) #kN
print(f'Fackverket deformerar plastiskt vid P = {max_P/1e3} kN')

min_A0 = round(A0*P/max_P*1e4,1)
print(f'Med P = 150 kN börjar fackverket deformera vid A0 = {min_A0} cm^2')

plt.show()