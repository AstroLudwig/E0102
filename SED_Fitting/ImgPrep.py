from astropy import constants as const
from astropy import units as u
from astropy.wcs import WCS, utils
from astropy.io import fits
from astropy.coordinates import SkyCoord

import numpy as np

import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import Axes3D

from scipy.stats import chisquare
from scipy.optimize import curve_fit

import sys
sys.path.insert(0, '../Convolve_Regrid')
import Convolve
import regrid
sys.path.insert(0, '../Prune')
import prune

##############
## Switches ##
##############

convolve_regrid_prune = True

# Testing
test_convolve = True
test_prune = True
###########
## Files ##
###########

f24 = '../Final_Files/24/24um_diff_3000_steps_snr.fits'
f70 = '../Final_Files/70/70um_diff_3000_steps_snr.fits'
f100 = '../Final_Files/100/100um_70modeled_ForceSlope_snr.fits'
f160 = '../Final_Files/160/160um_70modeledFS_prnd_snr.fits'

#kfile = '../Final_Files/KS/snr_24_convolved.fits'
Ofiles = [f24,f70,f100,f160]

#############
# Constants #
#############
# http://docs.astropy.org/en/stable/constants/
# Using CGS units

M_sun = const.M_sun.cgs# g 
pc = const.pc.cgs.value # cm
SMC_d = 61 * np.power(10.,3) * u.pc #Parsecs  
h_ = const.h.cgs
c_ = const.c.cgs
kb_ = const.k_B.cgs
xdim = [100,150]
ydim = [100,140]
title = ["24um","70um","100um","160um"]
#######################################
# Image Prep: Convolve/Regrid/Prune   #
#######################################

if convolve_regrid_prune: 
    # Kernels to Convolve with
    k24_160 = '../Convolve_Regrid/Kernels/Kernel_LoRes_MIPS_24_to_PACS_160.fits'
    k70_160 = '../Convolve_Regrid/Kernels/Kernel_LoRes_PACS_70_to_PACS_160.fits'
    k100_160 = '../Convolve_Regrid/Kernels/Kernel_LoRes_PACS_100_to_PACS_160.fits'

    #Convolve
    Convolve.master_convolve(k24_160, f24, '../Final_Files/24/24_SNR_Convolved.fits')
    Convolve.master_convolve(k70_160, f70, '../Final_Files/70/70_SNR_Convolved.fits')
    Convolve.master_convolve(k100_160, f100, '../Final_Files/100/100_SNR_Convolved.fits')

    #Regrid
    regrid.resample('../Final_Files/24/24_SNR_Convolved.fits',f160,'../Final_Files/24/Final_24_SNR_CR.fits')
    #regrid.resample(kfile,f160,'../Final_Files/24/Final_24_SNR_CR.fits')
    regrid.resample('../Final_Files/70/70_SNR_Convolved.fits',f160,'../Final_Files/70/Final_70_SNR_CR.fits')
    regrid.resample('../Final_Files/100/100_SNR_Convolved.fits',f160,'../Final_Files/100/Final_100_SNR_CR.fits')

    #Prune
    prune.prune('../Final_Files/24/Final_24_SNR_CR.fits','../Final_Files/24/Final_24_SNR_CR_Prnd.fits',np.nan)
    prune.prune('../Final_Files/70/Final_70_SNR_CR.fits','../Final_Files/70/Final_70_SNR_CR_Prnd.fits',np.nan)
    prune.prune('../Final_Files/100/Final_100_SNR_CR.fits','../Final_Files/100/Final_100_SNR_CR_Prnd.fits',np.nan)

RCfiles = ['../Final_Files/24/Final_24_SNR_CR.fits','../Final_Files/70/Final_70_SNR_CR.fits','../Final_Files/100/Final_100_SNR_CR.fits',f160]
files = ['../Final_Files/24/Final_24_SNR_CR_Prnd.fits','../Final_Files/70/Final_70_SNR_CR_Prnd.fits','../Final_Files/100/Final_100_SNR_CR_Prnd.fits',f160]

if test_convolve:
    f, axes = plt.subplots(2,2,sharey=False)

    counter = 0     
    for i in range(2):
        for j in range(2):
            data = fits.open(RCfiles[counter])[0].data
            select = data[ydim[0]:ydim[1],xdim[0]:xdim[1]]
            axes[i,j].imshow(data,vmin=np.nanmin(select),vmax = np.nanmax(select))
            axes[i,j].set_xlim(xdim)
            axes[i,j].set_ylim(ydim)
            axes[i,j].set_title(title[counter])
            counter += 1
if test_prune:    
    g, axes = plt.subplots(2,2,sharey=False)

    counter = 0     
    for i in range(2):
        for j in range(2):
            data = fits.open(files[counter])[0].data
            select = data[ydim[0]:ydim[1],xdim[0]:xdim[1]]
            axes[i,j].imshow(data,vmin=np.nanmin(select),vmax = np.nanmax(select))
            axes[i,j].set_xlim(xdim)
            axes[i,j].set_ylim(ydim)
            axes[i,j].set_title(title[counter])
            counter += 1

###########
# Error   #
###########
# Std in Sky Remove Region, Goes in order 24 70 100 160
Noise = np.loadtxt('../Sky_Remove/Sigma.txt')
Calibr = np.repeat(.1,4)
BkgRemov = np.array([0.004714082,0.078627909,.25,.5])   ## This is Probably Different Now Check With Karin
# Input Measured Value Array in Order: 24,70,100,160
# Outputs Array of Errors in the same order
def error(vals):
    err = []
    for i in range(4):
        X = vals[i]
        quad = (Noise[i])**2 + (Calibr[i]*X)**2 + (BkgRemov[i]*X)**2
        err.append(np.sqrt(quad))
    return np.asarray(err)
################
# Pixel Area   #
################

#Calculate the pixel area using cdelt from the header and x = theta * d
#where d is the distance to the smc (61 kpc)
#Area is then radian^2 (or steradian) * distance_to_the_smc^2

# Units are in sq parsec

def pixarea(filename):                      
    hdr = fits.open(filename)[0].header
    cdelt1 = np.abs(hdr["CDELT1"]) # Assuming Pixel is Square
    cdelt2 = np.abs(hdr["CDELT2"])
    cd1 = (cdelt1 * u.deg).to(u.rad)
    cd2 = (cdelt2 * u.deg).to(u.rad)
    sr = (cd1 * cd2) # Steradian
    area = sr * np.power(SMC_d,2) / u.rad**2
    
    return area.value 
##############
# Get Data   #
##############
data = []; sums = []; areas = []; means = []
for i in range(4):
    data.append(fits.open(files[i])[0].data)
    sums.append(np.nansum(fits.open(files[i])[0].data))
    areas.append(pixarea(files[i]))
    means.append(np.nanmean(fits.open(files[i])[0].data))
##########
# Plot   #
##########

plt.figure(3)
plt.scatter([24,70,100,160],means)
plt.title("Mean Value in each subtraction.")
plt.ylabel("Means")
plt.xlabel("Wavelength (Micron)")

plt.show()