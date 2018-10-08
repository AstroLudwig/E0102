# -*- Copyright (c) 2018, Bethany Ann Ludwig, All rights reserved. -*-
"""
NAME:
	Fix Saturation
PURPOSE:
	In the Spitzer 24 micron image a bright point source oversaturations and befuddles convolution efforts.
	To correct for this the region of this point source is set to nan.
"""

from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np

# File Handling
f = fits.open('../Original_Files/Infrared/n76crop_smcsage_24.fits')[0].data
hdr = fits.open('../Original_Files/Infrared/n76crop_smcsage_24.fits')[0].header

# View Image
plt.figure(1)
plt.imshow(np.flipud(f),vmin=0,vmax=10)

# Define Oversaturated Region
x = np.arange(86,105)
y = np.arange(211,228)

x,y = np.meshgrid(x,y)

# Replace with NaN
f[y,x] = np.nan

# View Image Again
plt.figure(2)
plt.imshow(np.flipud(f),vmin=0,vmax=10)

# Save
fits.writeto('../Saturation_Correction/n76crop_smcsage_24_desaturated.fits',f,hdr,overwrite='true')

plt.show()
