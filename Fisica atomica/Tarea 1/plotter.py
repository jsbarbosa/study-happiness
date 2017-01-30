import numpy as np
import matplotlib.pyplot as plt

N = 1000; n = 5
x = np.zeros(N); v = np.zeros(N)
ks = np.linspace(0.1, 1, n)
dt = 0.1
t = np.linspace(0, (N-1)*dt, N)

x[0] = 1
v[0] = 1

def equation(x):
    return -k*x

def solver():
    for i in range(1,N):
        v[i] = v[i-1] + equation(x[i-1])*dt
        x[i] = x[i-1] + v[i]*dt
    plt.plot(t*k, x, label="$%.2f$ "%k)
    
for i in range(n):
    k = ks[i]
    solver()
plt.legend(title = "$k$ (s$^{-2}$) values")
plt.ylabel("$x$ (m)")
plt.xlabel("$kt$ (s$^{-1}$)")
plt.xlim(0, 10)
plt.grid()
plt.savefig("plot.pdf")
