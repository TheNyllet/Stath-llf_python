import numpy as np
import matplotlib.pyplot as plt
from utils import *

## Indata (Geometry, meterial)
E = 2.1e11 # N/m^2
A0 = 7.85e-3 # m^2
P = 1,5e5 # N
L = 2 # m

## Topology
Edof = np.array([
    [0, 0, 2, 2], 
    [0, 0, 2, 2], 
    [0, 0, 2, 2], 
    [0, 0, 2, 2], 
    [2, 2, 2, 2], 
    [2, 2, 2, 2], 
    [2, 2, 2, 2], 
    [2, 2, 2, 2], 
    [2, 2, 2, 2], 
    [2, 2, 2, 2], 
    [2, 2, 2, 2]
], dtype=int )

# Koordinater för varje nod:
Coord = np.array([
    [0.0,0.0], [1.0,0.0], [0.0,2.0], [1.0,2.0], [1.0,3.0], [3.0,2.0], [3.0,3.0]
])

# x-koordinater för varje element
Ex = np.array([
    [1.0,1.0], [1.0,0.0], [0.0,1.0], [0.0,0.0], [0.0,1.0], [0.0,1.0], [1.0,1.0], [1.0,3.0], [1.0,3.0], [1.0,3.0], [3.0,3.0]
])

# y-koordinater för varje element
Ey = np.array([
    [0.0,2.0], [0.0,2.0], [0.0,2.0], [0.0,2.0], [2.0,2.0], [2.0,3.0], [2.0,3.0], [2.0,2.0], [3.0,2.0], [3.0,3.0], [2.0,3.0]
])

#Hjälpvariabler:
nel = len(Ex)  # Antal element
ndofs = len(Coord) - 4 #Totalt antal frihetsgrader

# Plot mesh (tips: använd eldraw2 i utils.py)


# Fördefinera styvhetsmatrisen och kraftvektorn
K = np.zeros((4,4))
f = np.zeros(4)

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
        [c**2, c*s, -c**2, -c*s],
        [c*s, s**2, -c*s, -s**2],
        [-c**2, -c*s, c**2, c*s],
        [-c*s, -s**2, c*s, s**2]
    ])

    #Assemblera in element styvhetsmatrisen och globala matrisen
    K = assem(Edof[el, :], K, Ke)

    print(K)

# Lägg till kraften P i lastvektorn:
...


# Lös ekvations systemet (: använd solveq i utils.py)
...


# Plotta deformerad mesh (: använd eldisp2 i utils.py)
...

# Räkna ut krafter och spänningar i varje element
for el in range(nel):
    print('hej')
    #... tips: använd bar2s i utils.py

#osv...