import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

calibration = fits.open('raw/calibracion.FIT')
fig, (ax1, ax2) = plt.subplots(2, sharex=True)

data = calibration[0].data
integral = np.sum(data, axis=0)
ax1.imshow(data, cmap='gray')
ax2.plot(integral)
plt.show()
