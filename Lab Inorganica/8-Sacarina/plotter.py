import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt('UV-vis.csv', delimiter=',')
data = data[3:]

labels = ['Cobalto', 'Cobre']

for i in range(len(labels)):
    x = data[:, 2*i]
    y = data[:, 2*i+1]
    
    try:
        nans = np.argwhere(np.isnan(x))[0][0]
    except:
        nans = -1
    
    x = x[:nans]
    y = y[:nans]
    
    pos = y.argsort()[::-1][0]
    plot = plt.plot(x, y, label = labels[i])[0]
    
    plt.text(x[pos] - 30, y[pos] + 0.01, '%.3f'%x[pos])
    plt.fill_between(x, 0, y, color = plot.get_color(), alpha = 0.3)
    
plt.grid()
plt.legend(loc='upper left')
plt.ylabel('Absorbancia (u.a.)')
plt.xlabel('Longitud de onda (nm)')

plt.savefig('images/absorbances.pdf')
