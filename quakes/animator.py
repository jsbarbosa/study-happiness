#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  animator.py
#  
#  Copyright 2016 Juan Barbosa <juan@Lenovo-U410>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.

import numpy as np
from textwrap import wrap
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.basemap import Basemap

def dataloader(name = "signif.txt"):
    data = np.genfromtxt(name, delimiter="\t", usecols=[2, 16, 20, 21], skip_header=True)
    nans = np.where(data[:,1] != data[:,1])[0]
    data[nans, 1] = np.mean(data[:, 1])
    
    years = data[:,0]
    intensities = data[:, 1]
    x = data[:,3]
    y = data[:,2]
    
    n = len(x)

    
    cite = "National Geophysical Data Center / World Data Service (NGDC/WDS): Significant Earthquake Database. National Geophysical Data Center, NOAA. doi:10.7289/V5TD9V7K"
    title = "\n".join(wrap(cite, 70))
    
    figure = plt.figure(figsize=(12,6))
    
    my_map = Basemap(projection='robin', lon_0=0, resolution='l')

    my_map.bluemarble()
    my_map.drawcountries()
    
    xs, ys = my_map(x, y)
    points = my_map.scatter(xs, ys, marker='o', alpha = 0.2, s = intensities*5, c = intensities, lw = 0.1)
    text = plt.annotate('', xy=(0.45, 0.01), xycoords='axes fraction')
    
    points.set_offsets(np.zeros((1,2)))
    fps = 1/20
    duration = 60

    step = int(np.ceil(n*fps/duration))
    frames = int(n/step)
    fps = frames/duration
    
    data = np.hstack((xs[:,np.newaxis], ys[:, np.newaxis]))
    
    def update(i):
        p = i*step
        points.set_offsets(data[:p])
        if years[p] == years[p]:
            text.set_text("Year: %d"%(years[p]))
        return points, text
    
    cbar = plt.colorbar(points)
    cbar.ax.set_xlabel('Intensity')
    
    plt.title(title, fontsize=10)
    plt.tight_layout()

    ani = animation.FuncAnimation(figure, update, frames = frames)
    
    return ani, fps

if __name__ == "__main__":
    ani, fps = dataloader()
    ani.save('Quakes.mp4', dpi=180, codec='h264', fps = fps)
#    ani.save("Quakes.gif", writer = "imagemagick", dpi = 45)
#    plt.show()
