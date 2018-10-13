# -*- Copyright (c) 2018, Bethany Ann Ludwig, All rights reserved. -*-
"""
NAME:
    Diffusion Plot Negatives
PURPOSE:
    Quantify Oversubtraction from Diffusion Method  
"""

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from matplotlib.widgets import Slider

#################
# File Handling #
#################

steps_start = 5
steps = 9995
s = np.arange(steps_start,steps,5)

im24 = False
im70 = False

if im24:
    pref_sub = 'im24/snr24/24um_diff_'
    app_sub = '_steps_snr.fits'
    # Define Region for SNR
    xdim = [190,230]; ydim = [130,180]     
if im70:
    pref_sub = 'im70/snr70/70um_diff_'
    app_sub = '_steps_snr.fits'
    # Define Region for SNR
    xdim = [210,275]; ydim = [200,260] 

#############################
# Quantify Oversubtraction  #
#############################

def count_negs(data):
    x = np.arange(xdim[0],xdim[1])
    y = np.arange(ydim[0],ydim[1])
    x,y = np.meshgrid(x,y)
    newdata = data[y,x]

    # Get only finite data
    a,b = np.where(np.isnan(newdata))
    newdata[a,b] = 0.
    
    # Count all the negatives
    e,count = np.shape(np.where(newdata< 0))
    c,d = np.where(newdata < 0)
    # Get the absolute sum
    sum_n = np.sum(np.abs(newdata[c,d]))
    
    # Return Stats
    return count, sum_n, np.mean(newdata)

# Get an array of stats for each time step
neg_sums = []; neg_counts = [];totals = []; mean_region = []
for i in range(steps_start,steps,5):
    with fits.open(pref_sub+str(i)+app_sub) as hdul:
        neg_count,neg_sum, avg = count_negs(hdul[0].data)
        mean_region.append(avg)
        neg_sums.append(neg_sum)
        neg_counts.append(neg_count)
        totals.append(np.nansum(hdul[0].data))

# How many values are there in the region?
pix_count = (xdim[1]-xdim[0])*(ydim[1]-ydim[0])

###############
# PLOT RESULT #
###############

f, (ax,bx,cx) = plt.subplots(1,3)
ax.plot(s,neg_sums,c='b')
ax.plot(s,totals,c='Orange')
ax.grid
ax.set_title("Total Sum(Orange), Negative Absolute Sum (Blue)")
ax.set_xlabel("Time Steps")
bx.plot(s,neg_counts)
bx.set_title("Total Sum(Orange), Negative Count (Blue)")
bx.set_xlabel("Time Steps")
bx.plot(s,totals,c='Orange')
cx.plot(s,neg_counts)
cx.set_title(("Negative Count, {:2f} % max negative pixels in region").format(np.max(neg_counts)/pix_count*100))


plt.show()