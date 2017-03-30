import numpy as np
import matplotlib.pyplot as plt

Z = 2
U, R, P, Q, V = np.genfromtxt("data.dat").T
logU = np.log(U)
logV = np.log(-V)
less = 2/U
greater = 2/(Z*U)
less = np.log(less)
greater = np.log(greater)

fig, axes = plt.subplots(2, 2)
axes = axes.reshape(4)

axes[0].plot(U, R, label = "Radial")
axes[1].plot(U, P, label = "Probability")
axes[2].plot(U, Q, label = "Charge")
axes[3].plot(logU, logV, label = "Potential")
axes[3].plot(logU, less, "--", c="g")
axes[3].plot(logU, greater, "--", c="g")
for i in range(4):
    axes[i].legend()
plt.show()
