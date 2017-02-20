import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt("functions.dat")
energies = np.genfromtxt("energies.dat")
omegas = np.genfromtxt("omegas.dat")

N = len(omegas)
n = 10000
dx = 0.001
x = np.linspace(0, (n-1)*dx, n)
data_split = np.split(data, N)
energies_split = np.split(energies, N)
answers = int(data_split[0].shape[0]/n)

fig, axes = plt.subplots(2, 2, sharex=True, sharey=True)
axes = axes.reshape(4)
for j in range(N):
    data = np.split(data_split[j], answers)
    for i in range(answers):
        p = axes[j].plot(x, data[i], label="E = %.3f"%energies_split[j][i], lw = 1)[0]
        if i%2 == 0:
            axes[j].plot(-x, -data[i], "--", lw = 1, c = p.get_color()) 
        else:
            axes[j].plot(-x, data[i], "--", lw = 1, c = p.get_color()) 
    axes[j].set_ylim(-2, 2)
    axes[j].set_xlim(-5, 5)
    axes[j].legend(loc = 1, fontsize=5)
    if j%2 == 0:
        axes[j].set_ylabel("$\Psi(x)$")
    if j > 1:
        axes[j].set_xlabel("$x$")
fig.tight_layout()
plt.savefig("plot.pdf")


fig = plt.figure()

for i in range(N):
    energies = energies_split[i]
    n = np.arange(1, answers + 1)
    A, B = np.polyfit(n, energies, 1)
    x = np.linspace(1, 7, 10)
    y = A*x + B
    p = plt.plot(n, energies, "o", label = "$\omega = %.1f$"%omegas[i])[0]
    plt.plot(x, y, "--", c=p.get_color())
    plt.text(i+1.2, i+1.5, "$E = %.2fn+%.2f$"%(A, B), color=p.get_color())
    
plt.xlabel("$n$")
plt.ylabel("$E$")
plt.legend()
plt.savefig("energies.pdf")

