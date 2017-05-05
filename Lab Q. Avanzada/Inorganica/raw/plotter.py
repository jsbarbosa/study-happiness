import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt("spectrum.txt")
fig, ax = plt.subplots()
ax.plot(data[:,0], data[:,1])
ax.set_xlim(0, 10)
ax.invert_xaxis()
plt.show()
