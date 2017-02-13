import numpy as np
import matplotlib.pyplot as plt
#from scipy.signal import find_peaks_cwt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset

plt.rcParams.update({'font.size': 8})
plt.rc('grid', linestyle="--", color='black', alpha=0.2)
data = np.genfromtxt("benzoina.txt", delimiter="\t")
x = data[:,0]
y = data[:,1]
dpi = 100#300
fig, ax = plt.subplots(dpi=dpi)

locations = [1, 7, 4]
subplots = []
xlims = [(7.2, 8.0), (3.65, 3.75), (1.18, 1.28)]
ylims = [max(y)*1.01, 500, 800]
     
for i in range(3):
    sub = inset_axes(ax, width=2, height=0.9, bbox_to_anchor=[6*dpi, dpi*(4.5 - 1.10*(i))])
    mark_inset(ax, sub, loc1=4, loc2=3, ec="r", alpha = 0.3, zorder = -1)
    sub.plot(x, y, lw=0.5, c='black')
#    min_index = abs(xlims[i][0] - x).argmin()
#    max_index = abs(xlims[i][1] - x).argmin()
#    pk = find_peaks_cwt(y[min_index:max_index], np.arange(13, 20)) + min_index
#    pos = y[pk].argsort()[::-1]
#    xs = [x[pk[pos[0]]], x[pk[pos[1]]]]
#    ys = [y[pk[pos[0]]]*1.1, y[pk[pos[0]]]*1.1]
#    sub.plot(xs, ys, "-o", ms = 3)
#    sub.text(xs[1], ys[0]*1.2, "%.3f"%(xs[1]-xs[0]))
    sub.set_ylim(0, ylims[i])
    sub.set_xlim(xlims[i][0], xlims[i][1])
    sub.invert_xaxis()
    
ax.plot(x, y, lw=0.5, c="black")
ax.set_xlim(-1, 10.5)
ax.invert_xaxis()
ax.set_ylabel("Intensidad (u.a)")
ax.set_xlabel("$\delta$ (ppm)")
plt.suptitle("$^1$H-RMN Benzoina", fontsize=12)
plt.show()
#fig.savefig("H-Benzoina.png")
