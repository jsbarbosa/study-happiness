import numpy as np
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
n = 6
nl = 4

energies = np.zeros((nl, n))
functions = np.zeros((nl, n, N))

for l in range(nl):
    if l == 0:
        e = -1.0
        y = [-1, 1]
    else:
        e = E[0]*1.1 
        y = [10**(-l+1), 0]
    E = np.zeros(n)
    R = np.zeros((n, N))
    for i in range(n):
        de = abs(e)*0.01
        last = 0
        current = 10
        thresh = 0.09
        while abs(current) > thresh and e < 0 and de > 1e-9:
            e += de
            r = solver(y, U, l, e, du)
            current = r[-1]
            if current*last < 0:
                e += -de
                de *= 0.1
            last = current
        E[i] = e
        R[i] = r
        e += 20*de
    energies[l] = E
    functions[l] = R

fix, axes = plt.subplots(2, 2, sharex=True, figsize=(12, 6))
axes = axes.reshape(4)
for l in range(nl):
    y = functions[l]
    for i in range(n):
        axes[l].plot(U, y[i], label="$\epsilon=%f$"%energies[l, i])
    if l == 0 or l == 2:
        axes[l].set_ylabel("$R(U)$")
    if l > 1:
        axes[l].set_xlabel("$U$")
    axes[l].set_title("$l=%d$"%l)
    axes[l].set_ylim(min(y[0]), max(y[1]))
    axes[l].legend(fontsize=8)
plt.savefig("R.pdf")

n_n = np.arange(2, n+2)
E_T = -1/n**2

L = np.zeros((nl, n+2))
for l in range(1, 3):
    L[l+1, l:n+l] = energies[1+l]

L[0, :n] = energies[0]
L[1, :n] = energies[1]
L[2, 1:n+1] = energies[2]
L[3, 2:n+2] = energies[3]
for i in range(n+2):
    values = L[:, i]
    values = [str(item) if item != 0 else "" for item in values]
    text = " & ".join(values)        
    print(r"%d & %f & %s \\"%(i+2, -1/(i+2)**2,text))

temps = np.arange(0, 400)
thresh = [0.01, 0.01, 0.25, 5]
fix, axes = plt.subplots(2, 2, sharex=True, figsize=(12, 6))
axes = axes.reshape(4)
for l in range(nl):
    rl = functions[l]
    th = thresh[l]
    for i in range(n):
        r = rl[i]
        pos = np.where(abs(r[400:]) < th)[0] + 400
        pos = np.concatenate((temps, pos))
        u = U[pos]
        r = r[pos]
        P = 4*np.pi*(u*r)**2
        P *= 1/np.trapz(P)
        axes[l].plot(u, P, label="$\epsilon=%f$"%energies[l, i])
    if l == 0 or l == 2:
        axes[l].set_ylabel("$P(U)/\int P(u)dU$")
    if l > 1:
        axes[l].set_xlabel("$U$")
    axes[l].set_title("$l=%d$"%l)
    axes[l].legend(fontsize=8)
plt.savefig("P.pdf")



