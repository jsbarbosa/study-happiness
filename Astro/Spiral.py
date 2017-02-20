import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

def galaxy(arms = 5):
    N = 10000
    a = 0.6
    b = 0.5
    t = np.random.randn(N)
    x = np.array([0])
    y = np.array([0])

    for i in range(arms):
        xt = a*np.exp(b*t)*np.cos(t+2*i*np.pi/arms) + np.random.normal(0, a*0.2, N)
        yt = a*np.exp(b*t)*np.sin(t+2*i*np.pi/arms) + np.random.normal(0, a*0.2, N)
        x = np.hstack((x, xt))
        y = np.hstack((y, yt))
        
    return x, y

x, y = galaxy()   
frames = 120
phi = 2*np.pi/frames

fig, ax = plt.subplots(facecolor='black')
ax.set_axis_off()
ax.set_facecolor('black')
points, = ax.plot(x, y, "o", ms=0.1, c="black")
ax.set_xlim(-2, 2)
ax.set_ylim(-2, 2)
def update(i):
    xt = x*np.cos(i*phi) - y*np.sin(i*phi)
    yt = x*np.sin(i*phi) + y*np.cos(i*phi)
    points.set_data(xt, yt)
    return points,
    
ani = FuncAnimation(fig, update, interval=frames, frames = frames)
ani.save('ani.gif', writer='imagemagick', dpi=72)


