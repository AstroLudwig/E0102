# -*- coding: utf-8 -*-
"""
In the new and larger 24 micron image a bright point source oversaturations and befuddles convolution efforts.
To correct for this the region of this point source is set to nan.
"""

from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
f = fits.open('../Original_Files/Infrared/n76crop_smcsage_24.fits')[0].data
hdr = fits.open('../Original_Files/Infrared/n76crop_smcsage_24.fits')[0].header

plt.imshow(np.flipud(f),vmin=0,vmax=10)

x = np.arange(86,105)
y = np.arange(211,228)

x,y = np.meshgrid(x,y)

f[y,x] = np.nan

plt.figure()
plt.imshow(np.flipud(f),vmin=0,vmax=10)
plt.show()

fits.writeto('../Saturation_Correction/n76crop_smcsage_24_desaturated.fits',f,hdr,overwrite='true')
