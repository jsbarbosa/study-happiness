import numpy as np
import matplotlib.pyplot as plt

N = 1000; n = 5
x = np.zeros(N); v = np.zeros(N)
ks = np.linspace(0.1, 1, n)
dt = 0.1
t = np.linspace(0, (N-1)*dt, N)

x[0] = 1
v[0] = 1

fig, axes = plt.subplots(1, 2, sharey = True, figsize=(8, 5))
ax1, ax2 = axes

def equation(x):
    return -k*x

def solver():
    for i in range(1,N):
        v[i] = v[i-1] + equation(x[i-1])*dt
        x[i] = x[i-1] + v[i]*dt
    ax1.plot(t, x)
    ax2.plot(t*np.sqrt(k), x, label="$%.2f$ "%k)
    
for i in range(n):
    k = ks[i]
    solver()

for ax in axes:
    ax.set_xlim(0, 20)
    ax.grid()

plt.legend(title = "$k$ (s$^{-2}$) values")
ax1.set_ylabel("$x$ (m)")
ax1.set_xlabel("$t$ (s)")
ax2.set_xlabel("$t\sqrt{k}$ (rad)")

n = int(np.ceil(max(t)/np.pi))
xticks = [i*np.pi for i in range(n)]
xticks_labels = ['$0$', '$\pi$'] + ['$%d\pi$'%i for i in range(2, n)]
plt.xticks(xticks,xticks_labels)
ax2.set_xlim(0, 10*np.pi)
plt.savefig("plot.pdf")
