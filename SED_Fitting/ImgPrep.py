# -*- Copyright (c) 2018, Bethany Ann Ludwig, All rights reserved. -*-
"""
NAME:
    Image Prep
PURPOSE:
    Before SED Fitting all files must be convolved and regridded to the lower res 160 image. 
    Data is also reduced (pruned) so the only remaining data is specific to the SNR.      
"""
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt 

# My Modules
import sys
sys.path.insert(0, '../Convolve_Regrid')
import Convolve, Regrid
sys.path.insert(0, '../Background_Removal/Bootstrap/')
import prune

##############
## Switches ##
##############
# Testing
plot_ConvolveRegrid = True
plot_Prune = True

###########
## Files ##
###########

# Get Final SNR Images from appropriate folders
# Must run previous background removal steps first.
f24 = '../Background_Removal/Diffusion/im24/snr24/24um_diff_3000_steps_snr.fits'
f70 = '../Background_Removal/Diffusion/im70/snr70/70um_diff_3000_steps_snr.fits'
f100 = '../Background_Removal/Bootstrap/100um/100um_70modeled_snr.fits'
f160 = '../Background_Removal/Bootstrap/160um/160um_70modeled_snr.fits'

Files = [f24,f70,f100,f160]

#######################
# Plotting Parameters #
#######################
xdim = [100,150]
ydim = [100,140]
title = ["24um","70um","100um","160um"]
#######################################
# Image Prep: Convolve/Regrid/Prune   #
#######################################

# Kernels to Convolve with
k24 = '../Convolve_Regrid/Kernels/Kernel_LoRes_MIPS_24_to_PACS_160.fits'
k70 = '../Convolve_Regrid/Kernels/Kernel_LoRes_PACS_70_to_PACS_160.fits'
k100 = '../Convolve_Regrid/Kernels/Kernel_LoRes_PACS_100_to_PACS_160.fits'
Kernels = [k24,k70,k100]

# File Save Names
s24 = '../Final_Files/24/24_SNR_'
s70 = '../Final_Files/70/70_SNR_'
s100 = '../Final_Files/100/100_SNR_'
s160 = '../Final_Files/160/160_SNR_'
Save = [s24,s70,s100,s160]

for i in range(3):
#Convolve
    Convolve.master_convolve(Kernels[i],Files[i], Save[i]+'Convolve.fits')
#Regrid
    Regrid.resample(Save[i]+'Convolve.fits',f160,Save[i]+'Convolve_Regrid.fits')
#Prune
    prune.SavePrune(Save[i]+'Convolve_Regrid.fits',Save[i]+'Convolve_Regrid_Prune.fits')
prune.SavePrune(f160,Save[3]+'Prune.fits')

if plot_ConvolveRegrid:
    f, axes = plt.subplots(1,4,sharey=False)
    
    for i in range(4): 
        if i == 3: 
            file = f160
        else:
            file = Save[i]+'Convolve_Regrid.fits'
        data = fits.open(file)[0].data
        select = data[ydim[0]:ydim[1],xdim[0]:xdim[1]]
        axes[i].imshow(data,vmin=np.nanmin(select),vmax = np.nanmax(select))
        axes[i].set_xlim(xdim); axes[i].set_ylim(ydim)
        axes[i].set_title(title[i])

if plot_Prune:    
    g, axis = plt.subplots(1,4,sharey=False)

    for i in range(4):  
        if i == 3: 
            file = Save[3]+'Prune.fits'
        else:
            file = Save[i]+'Convolve_Regrid_Prune.fits'     
        data = fits.open(file)[0].data
        select = data[ydim[0]:ydim[1],xdim[0]:xdim[1]]
        axis[i].imshow(data,vmin=np.nanmin(select),vmax = np.nanmax(select))
        axis[i].set_xlim(xdim); axis[i].set_ylim(ydim)
        axis[i].set_title(title[i])
plt.show()        