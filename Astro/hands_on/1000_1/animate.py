import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.set_aspect("equal")
#ax.set_axis_off()
x, y, z = np.genfromtxt("0_data.dat")
N = int(len(x)/2)
plot1 = ax.plot([], [], [], "o", ms=0.5, c="g", alpha = 0.5)[0]
plot2 = ax.plot([],[],[], "o", ms=0.5, c="b", alpha = 0.5)[0]
text = ax.text(0,0,0, "", color = "r")
ax.set_xlabel("$x$")
ax.set_ylabel("$y$")
ax.set_zlabel("$z$")

frames = 99
ratio = int(frames/frames)

def update(i):
    temp = np.genfromtxt("%d_data.dat"%(i+1))
    plot1.set_data(temp[:N,0], temp[:N,1])
    plot1.set_3d_properties(temp[:N, 2])
    plot2.set_data(temp[N:,0], temp[N:,1])
    plot2.set_3d_properties(temp[N:, 2])
#    text.set_text("%d"%i)
    ax.view_init(0, i*360/frames)

min_value, max_value = min(x), max(x)
min_value *= 2.0
max_value *= 2.0
ax.set_xlim(min_value, max_value)
ax.set_ylim(min_value, max_value)
ax.set_zlim(min_value, max_value)
fig.tight_layout()
ani = FuncAnimation(fig, update, frames=frames, interval=frames)
#ani.save("100.gif", writer='imagemagick')
plt.show()
