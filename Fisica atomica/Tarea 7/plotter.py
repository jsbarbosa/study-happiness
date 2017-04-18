import numpy as np
import matplotlib.pyplot as plt

Z = 2

fig, ax = plt.subplots()
for i in range(1, 11):
    data = np.genfromtxt("%d_data.dat"%i).T
    N = data.shape[1]
    data = data[:, ::int(N/100)]
    U, R, P, Q, V = data
    ax.plot(U, np.log(-V), label = "%d"%i)

ax.set_xlabel("$r/a_0$")
ax.set_ylabel("$\log(V/E_0)$")

plt.legend()
plt.savefig("change.pdf")
#plt.show()
    
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
