import sys
sys.path.insert(0, '../SED_Fitting/')

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt 
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
import Eqs
import seaborn as sns
from scipy.stats import norm
import matplotlib.mlab as mlab
from scipy.stats import iqr
##############
## Switches ##
##############

show = False
save = True


###########
## Files ##
###########

files = ['../Final_Files/24/24_SNR_Convolve_Regrid_Prune.fits','../Final_Files/70/70_SNR_Convolve_Regrid_Prune.fits',
		'../Final_Files/100/100_SNR_Convolve_Regrid_Prune.fits','../Final_Files/160/160_SNR_Prune.fits']

#############
# Constants #
#############
# http://docs.astropy.org/en/stable/constants/
# Using CGS units

xdim = [120,137]
ydim = [112,128]
title = ["24um","70um","100um","160um"]

##################
# Plotting Range #
##################

lam, k_amc = Eqs.log_extrap('../SED_Fitting/Kappa/kappa_amc.dat')
lam, k_mg2 = Eqs.log_extrap('../SED_Fitting/Kappa/kappa_mg2sio4.dat')
wv = [24,70,100,160]
pix_area = Eqs.pixarea(files[0])

####################
# Parameter Ranges #
####################
# Create a grid to do a fit
ColdTemp = np.arange(2,70,1) # Kelvin
# Fractional Mass of the Sun using Dex
# Cold Mass
ColdMass = 10**np.arange(-4,0.1,.1)
# Warm Mass
WarmMass = 10**np.arange(-8,-3,.1)
# Error Intervals
# 68.3%, 95.4% and 99.73%, or 1, 2 and 3 sigma 
intervals = [1,4,9]
##############
# Get Data   #
##############
# Get SNR data, the physical area of each SNR, and the average pixel intensity
data = []; Areas = []; AverageIntensities = np.zeros(4)
for i in range(4):
	AverageIntensity, Area = Eqs.AverageSED(files[i])
	data.append(fits.open(files[i])[0].data)
	Areas.append(Area)
	AverageIntensities[i] = AverageIntensity

# Get the error associated with each intensity
AverageError = Eqs.error(AverageIntensities)

#############################
# Load pixel by pixel files #
#############################

Pix_Warm_SED =  np.load("../SED_Fitting/Sols/PixbyPix/Warm_SED.npy")
Pix_Cold_SED =  np.load("../SED_Fitting/Sols/PixbyPix/Cold_SED.npy")

# Temperature                                         
Pix_Temp = np.loadtxt("../SED_Fitting/Sols/PixbyPix/Temperature.txt")

# Mass 
Pix_Warm_Mass = np.loadtxt("../SED_Fitting/Sols/PixbyPix/WarmMass.txt")
Pix_Cold_Mass = np.loadtxt("../SED_Fitting/Sols/PixbyPix/ColdMass.txt")   

# Chi Squared Maps                            
Pix_chisqrd = np.loadtxt("../SED_Fitting/Sols/PixbyPix/ChiSquared.txt")

sky_noise = np.loadtxt("../Sky_Remove/Sigma.txt")

Maps = [Pix_Temp,Pix_Warm_Mass,Pix_Cold_Mass]
Titles = ["Temperature","Warm Mass","Cold Mass"]
###################
# Create Template #
###################

# Confidence Interval
# ========================================================
sigma_fit = 0 # 0 is 1 sigma, 1 is 2 sigma, 2 is 3 sigma. 
# ========================================================
# Sets what multiple of the sky's noise to remove. 
# ========================================================
noise_multiplier = 3
# ========================================================
# Initialize array of empty 2d matrices the size of the images.
empty_template = np.zeros(np.shape(data[0])[0]*np.shape(data[0])[1]).reshape(np.shape(data[0])[0],np.shape(data[0])[1])
templates = [np.copy(empty_template),np.copy(empty_template),np.copy(empty_template),np.copy(empty_template)]

for i in range(4):
	for j in range(np.shape(data[i])[0]):
		for k in range(np.shape(data[i])[1]):
			if np.isfinite(data[i][j,k]):
				# If the data is greater than some multiple of the noise, keep it by setting it to 1. 
				if data[i][j,k] > noise_multiplier * sky_noise[i]: 
					templates[i][j,k] = 1

# Multiply all 4 templates together to get a final template				
template = np.ones(np.shape(data[0])[0]*np.shape(data[0])[1]).reshape(np.shape(data[0])[0],np.shape(data[0])[1])
for item in templates:
	template *= item 
# Turn off background by setting 0s's to nan's
template[np.where(template==0)] = np.nan					

##########
# Plot   #
##########

# Template
plt.figure()
plt.imshow(data[3]*template)
plt.xlim(xdim); plt.ylim(ydim)
plt.axis('off')
plt.title("Template * 160 microns",size="x-large")
if save:
	plt.savefig("Plots/Template.png",dpi=200)

# Maps
f, axes = plt.subplots(1,3,figsize=(3.75,1.5))
for i in range(3):
	axes[i].imshow(Maps[i]*template,vmin=np.nanmin(Maps[i]*template),vmax=np.nanmax(Maps[i]*template))
	axes[i].set_xlim(xdim)
	axes[i].set_ylim(ydim)
	axes[i].axis('off')
	axes[i].set_title(Titles[i])
plt.subplots_adjust(top=1, bottom=0, left=0.05, right=0.95, hspace=0, wspace=0.09)
if save:
	plt.savefig("Plots/Maps.png",dpi=200)

# Since we're saving plots we don't always need to open a GUI
if show:
	plt.show()
