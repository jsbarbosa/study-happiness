import numpy as np
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET
from mpl_toolkits.mplot3d import axes3d
from matplotlib.animation import FuncAnimation


def sphere(x, y, z, radious, N=50):
    r = radious
    theta = np.linspace(0, np.pi, N)
    phi = np.linspace(0, 2*np.pi, N)
    theta, phi = np.meshgrid(theta, phi)
    x = r*np.sin(theta)*np.cos(phi) + x
    y = r*np.sin(theta)*np.sin(phi) + y
    z = r*np.cos(theta) + z
    return x, y, z

tree = ET.parse('BINAPcomplex_(1).mrv')
root = tree.getroot()[1][0][0]
atoms_ = root[0]
bonds_ = root[1]
N_atoms = 0
for item in atoms_:
    N_atoms += 1
N_bonds = 0
for item in bonds_:
    N_bonds += 1

atoms = {}
bonds = {}
colors = {"C": "black", "P": "sandybrown", "Cl": "green", "O": "red", "Ru": "blue", "R": "black", "H": "grey"}
sizes = {"C": 67, "P": 98, "Cl": 79, "O": 48, "Ru": 178, "R": 100, "H": 53}
for i in range(N_atoms):
    attrib = atoms_[i].attrib
    iden = attrib["id"]
    atoms[iden] = attrib
for i in range(N_bonds):
    attrib = bonds_[i].attrib
    iden = attrib["id"]
    bonds[iden] = attrib

bond_values = []
for (i, bond) in enumerate(bonds):
    bond = bonds[bond]
    ref = bond['atomRefs2']
    a1, a2 = ref.split(' ')
    a1 = atoms[a1]
    a2 = atoms[a2]
    x = [float(a1['x3']), float(a2['x3'])]
    y = [float(a1['y3']), float(a2['y3'])]
    z = [float(a1['z3']), float(a2['z3'])]
    bond_values.append([x, y, z])

fig, ax = plt.subplots(subplot_kw=dict(projection='3d'))
ax.set_axis_off()

for bond in bond_values:
    ax.plot(bond[0], bond[1], bond[2], color="k")
    
fig.tight_layout()
fig.subplots_adjust(left=-0.3, right=1.2, bottom=-0.25, top=1.2)

for atom in atoms:
    atom = atoms[atom]
    element = atom["elementType"]
    color = colors[element]
    size = sizes[element]/5
    x, y, z = float(atom['x3']), float(atom['y3']), float(atom['z3'])
    x, y, z = [x], [y], [z]
    ax.plot(x, y, z, "o", ms = size, color = color)

def update(i):
    angle = i*360/frames
    ax.view_init(30, angle)
    
frames = 60
ani = FuncAnimation(fig, update, frames=frames)
ani.save("noyori.gif", writer="imagemagick")
plt.show()
