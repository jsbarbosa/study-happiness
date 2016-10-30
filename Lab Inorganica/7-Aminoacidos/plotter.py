import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

data = np.genfromtxt('pKa.txt')

pH = data[:,0]
formation = data[:,1]

inter = interp1d(formation, pH)
x = inter(50)

plt.plot(pH, formation)
plt.plot(x, 50, "o", label="$pK_a=%.2f$"%x)
plt.xlabel('pH')
plt.ylabel('% formacion')
plt.grid()
plt.legend()
plt.savefig('images/pka_sim.pdf')
