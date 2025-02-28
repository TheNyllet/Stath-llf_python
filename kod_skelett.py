import numpy as np
import matplotlib.pyplot as plt
from utils import *

## Indata (Geometry, meterial)
E = 2.1e11 # N/m^2
A0 = 7.85e-3 # m^2
P = 1.5e5 # N
L = 2 # m

## Topology
Edof = np.array([
    [[3, 4], [7, 8]], 
    [[3, 4], [5, 6]], 
    [[1, 2], [7, 8]], 
    [[1, 2], [5, 6]], 
    [[5, 6], [7, 8]], 
    [[5, 6], [9, 10]], 
    [[7, 8], [9, 10]], 
    [[7, 8], [11, 12]], 
    [[9, 10], [11, 12]], 
    [[9, 10], [13, 14]], 
    [[11, 12], [13, 14]]
], dtype=int )

# Koordinater för varje nod:
Coord = L*np.array([
    [0.0,0.0], 
    [1.0,0.0], 
    [0.0,2.0], 
    [1.0,2.0], 
    [1.0,3.0], 
    [3.0,2.0], 
    [3.0,3.0]
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

#Hjälpvariabler:
nel = len(Ex)  # Antal element
ndofs = 2*len(Coord) #Totalt antal frihetsgrader

# Plot mesh (tips: använd eldraw2 i utils.py)


# Fördefinera styvhetsmatrisen och kraftvektorn
K = np.zeros((ndofs,ndofs))
f = np.zeros(ndofs)

# Assemblera elemented
for el in range(nel):
    x1, x2 = Ex[el]
    y1, y2 = Ey[el]
    L = np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)

    theta = np.arctan2(y2 - y1, x2 - x1)
    c = np.cos(theta)
    s = np.sin(theta)

    if el + 1 in [4,6,9]: A = A0 * 2
    else: A = A0
    
    #Räkna ut styvhetsmatrisen (Se föreläsningar eller Ekvation 11.18 i kursboken [Hållfasthetslära, Allmänna tillstånd])
    Ke = E*A/L*np.array([
        [1,-1],
        [-1,1]
    ])

    """
        [c**2, c*s, -c**2, -c*s],
        [c*s, s**2, -c*s, -s**2],
        [-c**2, -c*s, c**2, c*s],
        [-c*s, -s**2, c*s, s**2]
    """

    #Assemblera in element styvhetsmatrisen och globala matrisen
    K = assem(Edof[el], K, Ke)

    print(K)

# Lägg till kraften P i lastvektorn:
f[12] = -P

# Bestäm bcdofs och bcvals
bcdofs = np.array([i for i in range(1,15)])  # Exempel på frihetsgrader med randvillkor
bcvals = np.array([0.0, 0.0])  # Värden för randvillkoren

# Lös ekvationssystemet (: använd solveq i utils.py)
a, r = solveq(K, f, bcdofs, f)

# Plotta deformerad mesh (: använd eldisp2 i utils.py)
eldisp2(Ex, Ey, a)


# Räkna ut krafter och spänningar i varje element
for el in range(nel):
    N, sigma = bar2s(Ex[el], Ey[el], a[Edof[el, :]])
    print(f"Element {el}: Kraft = {N}, Spänning = {sigma}")
    #... tips: använd bar2s i utils.py
