from astropy.io import fits
import matplotlib.pyplot as plt 
import numpy as np 
import prune

data = prune.Prune('100um_70modeled/100um_70modeled_snr.fits',0,'nan')

plt.figure()
plt.hist(data[np.isfinite(data)])
plt.title("Histogram of pixel intensities in 100 micron img modeled by 70")

im = prune.Prune('100um_70modeled/100um_70modeled_snr.fits',0,'zero')
im[im<0.] = 1000
plt.figure()
plt.imshow(im)
plt.xlim(170,230)
plt.ylim(170,210)
plt.title("Where Negatives Are in 100 Micron SNR")
plt.show()
