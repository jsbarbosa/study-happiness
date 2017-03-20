import numpy as np
from scipy.constants import e, c, h
from scipy.interpolate import interp1d 
import matplotlib.pyplot as plt

def solver(y, u, l, epsilon, du):
    L = l*(l+1)
    prime = np.zeros_like(u)
    x = np.ones_like(u)
    prime[0], x[0] = y
    for i in range(N-1):
        derivative = - 2*prime[i]/u[i] - (epsilon + 2/u[i] - L/u[i]**2)*x[i]
        prime[i+1] = prime[i] + derivative*du
        x[i+1] = x[i] + prime[i+1]*du 
    return x

N = 1000
du = 100/N
U = np.linspace(du, N*du, N) 
s = 0.5
epsilons = [[-0.999], [-0.243794, -0.252106], [-0.109864, -0.112379, -0.111361]]
thresh = [[30], [300, 300], [300, 300, 300]]

R_VALUES = np.zeros((3, 3))
plt.figure()
for n in range(1, len(epsilons)+1):
    epsilon = epsilons[n-1]
    for l in range(len(epsilon)):
        if l == 0:
            y = [-1, 1]
        else:
            y = [1, 0]
        e_ = epsilon[l]
        R = solver(y, U, l, e_, du)
        pos = thresh[n-1][l]
        r = R[:pos]
        u = U[:pos]
        P = 4*np.pi*(r*u)**2
        temp = (np.trapz(P, dx = du))
        P *= 1/(np.trapz(P, dx = du))
        r_E = np.trapz(P*u, dx = du)
        R_VALUES[n-1, l] = r_E
        const = 3.622608e-4
        j1 = l+s
        j2 = l-s
        if j2 < 0: j2 = 0
        beta1 = j1*(j1+1)-l*(l+1)-s*(s+1)
        beta2 = j2*(j2+1)-l*(l+1)-s*(s+1)
        e_E = const*np.trapz(P/u**3, dx=du)
        e_E1 = beta1*e_E
        e_E2 = beta2*e_E
        print(r"%d & %d & %.3f & %.3e & %.3e \\"%(n, l, r_E, e_E1, e_E2))
        plt.plot(u, P, label="$n=%d$, $l=%d$"%(n, l))
plt.xlabel("$r/a_0$")
plt.ylabel("Probabilidad Radial ($4\pi R^2r^2$)")
plt.legend()
plt.savefig("probability.pdf")

plt.figure()
ns = np.arange(1, 4)
for l in range(len(epsilon)):
    plt.plot(ns**2, R_VALUES[:, l], "-o", label="$l=%d$"%l)
plt.xlabel("$n^2$")
plt.ylabel("$<r>/a_0$")
plt.legend()
plt.savefig("radius.pdf")
