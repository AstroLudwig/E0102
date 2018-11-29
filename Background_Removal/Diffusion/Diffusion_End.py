# -*- Copyright (c) 2018, Bethany Ann Ludwig, All rights reserved. -*-
"""
NAME:
    Diffusion End
PURPOSE:
    To decide which step the process should be haulted at. 
NOTE: 
    If starting from scratch you will need to create the appropriate folders for this to run.        
"""
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

f24 = True
f70 = False
get_data = True
noise_file = "../../Sky_Remove/Sigma.txt"

if f24: 
	prefix = "im24/snr24/24um_diff_"
	noise = np.loadtxt(noise_file)[0]
	RowDim = 301; ColDim = 401
if f70:
	prefix = "im70/snr70/70um_diff_"
	noise = np.loadtxt(noise_file)[1]
	RowDim = 462; ColDim = 484
appendix = "_steps_snr.fits"
End = 5000
# Get Data
if get_data:   
    # Create Cube
    SubtractionData = np.zeros(RowDim*ColDim*End).reshape(RowDim,ColDim,End)
    # Fill Cube
    for i in range(5,End,5):
        with fits.open(prefix+str(i)+appendix) as hdu:
        	SubtractionData[:,:,i] = hdu[0].data
# When does homogenization occur?
SumDifference = []
for i in range(5,int(End-5),5):
	# Take the sum of the difference of absolute pixel intensities.
	last = np.abs(SubtractionData[:,:,i+5]); first = np.abs(SubtractionData[:,:,i])
	SumDifference.append(np.nansum(first-last))
FirstSol = np.min(np.where(np.isclose(SumDifference,0,atol=noise))) * 5
# Report
print("Homogenization occurs at step:"+str(FirstSol))
print("The end sum is "+str(np.nansum(np.nansum(SubtractionData[:,:,FirstSol]))))
print(str(np.nansum(SubtractionData[:,:,5])-np.nansum(SubtractionData[:,:,FirstSol]))+" has been subtracted.")
# View
plt.figure()
plt.imshow(np.flipud(SubtractionData[:,:,FirstSol]))
plt.title("Final Subtracted Image")
plt.figure()
X = np.arange(5,End-5,5)
diff = np.abs(np.nansum(SubtractionData[:,:,FirstSol+5])-np.nansum(SubtractionData[:,:,FirstSol]))
plt.plot(X,SumDifference)
plt.scatter(FirstSol,diff,c="red")
plt.xlabel("Steps"); plt.ylabel("Difference from Previous Step")

plt.show()
