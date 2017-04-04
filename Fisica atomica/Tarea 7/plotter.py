import numpy as np
import matplotlib.pyplot as plt

Z = 2

energies = np.genfromtxt("energies.dat")
N = len(energies)+1
n = np.arange(1, N)

fig, (ax, ax1) = plt.subplots(2)
for i in range(1, 20, 4):
    data = np.genfromtxt("%d_data.dat"%i).T
    n_points = data.shape[1]
    data = data[:, ::int(n_points/100)]
    U, R, P, Q, V = data
    ax.plot(U, np.log(-V), label = "%d"%(i-1))

ax.set_xlabel("$r/a_0$")
ax.set_ylabel("$\log(V/E_0)$")

ax1.plot(n, energies, "-o", ms=2)
ax1.set_xlabel("Iteration")
ax1.set_ylabel("Energy/$E_0$")

ax.legend()
fig.tight_layout()
plt.savefig("change.pdf")
    
logU = np.log(U)
logV = np.log(-V)
less = 2/U
greater = 2/(Z*U)
less = np.log(less)
greater = np.log(greater)

fig, axes = plt.subplots(2, 2)
axes = axes.reshape(4)

axes[0].plot(U, R)
axes[1].plot(U, P)
axes[2].plot(U, Q)
axes[3].plot(logU, logV)
axes[3].plot(logU, less, "--", c="g")
axes[3].plot(logU, greater, "--", c="g")

axes[0].set_xlabel("$r/a_0$")
axes[1].set_xlabel("$r/a_0$")
axes[2].set_xlabel("$r/a_0$")
axes[3].set_xlabel("$\log(r/a_0)$")
axes[0].set_ylabel("$R$")
axes[1].set_ylabel("Probability")
axes[2].set_ylabel("Charge density")
axes[3].set_ylabel("$\log(V/E_0)$")
fig.tight_layout()
fig.savefig("complete.pdf")
#plt.show()
