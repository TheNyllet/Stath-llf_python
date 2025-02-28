import numpy as np
import matplotlib.pyplot as plt
from utils import *

# Kontrollera att coordxtr finns i utils
try:
    from utils import coordxtr
except ImportError:
    raise ImportError("Funktionen coordxtr saknas i utils. Se till att den finns definierad i utils.py.")

# Indata
E = 210e9  # Elasticitetsmodul (Pa)
A0 = 78.5e-4  # Tvärsnittsarea (m^2)
P = 150e3  # Last (N)
L = 2.0  # Baslängd (m)

# Topologi (stänger och knutpunkter)
Edof = np.array([
    [1, 2, 3, 4],
    [3, 4, 5, 6],
    [5, 6, 7, 8],
    [1, 2, 5, 6],
    [3, 4, 7, 8],
    [5, 6, 9, 10],
    [7, 8, 11, 12],
    [9, 10, 11, 12]
], dtype=int)

# Nodkoordinater
Coord = np.array([
    [0, 0], [L, 0], [2*L, 0], [0.5*L, L], [1.5*L, L], [0, 2*L], [2*L, 2*L]
])

# Importera coordxtr från utils
Ex, Ey = coordxtr(Edof, Coord, np.arange(1, len(Coord) * 2 + 1).reshape(-1, 2))

# Antal element och frihetsgrader
nel = len(Edof)
ndofs = 14

# Initialisera styvhetsmatris och lastvektor
K = np.zeros((ndofs, ndofs))
f = np.zeros(ndofs)

# Assemblera elementens styvhetsmatriser
for el in range(nel):
    x1, x2 = Ex[el]
    y1, y2 = Ey[el]
    L_el = np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)
    theta = np.arctan2(y2 - y1, x2 - x1)
    
    c = np.cos(theta)
    s = np.sin(theta)
    
    T = np.array([
        [c, s, 0, 0],
        [-s, c, 0, 0],
        [0, 0, c, s],
        [0, 0, -s, c]
    ])
    
    k_local = (E * A0 / L_el) * np.array([
        [1, -1],
        [-1, 1]
    ])
    
    Ke = T.T @ np.block([
        [k_local, np.zeros((2, 2))],
        [np.zeros((2, 2)), k_local]
    ]) @ T
    
    K = assem(Edof[el, :], K, Ke)

# Applicera last
f[4] = -P  # Last vid nod 3 i y-led

# Randvillkor (fixerade noder)
bcdofs = np.array([1, 2, 3, 4])  # Fixera båda noderna i botten (x och y-led)
bcvals = np.zeros_like(bcdofs, dtype=float)

# Lös ekvationssystemet
U, R = solveq(K, f, bcdofs, bcvals)

# Plotta deformerad struktur
Ed = extract_eldisp(Edof, U)
eldisp2(Ex, Ey, Ed, sfac=1)
plt.show()

# Beräkna stångkrafter och spänningar
for el in range(nel):
    ep = [E, A0]
    es = bar2s(Ex[el], Ey[el], ep, Ed[el])
    print(f"Stång {el+1}: Kraft = {es:.2f} N")
