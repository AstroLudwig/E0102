# -*- Copyright (c) 2019, Bethany Ann Ludwig, All rights reserved. -*-
"""
NAME:
	Data Visualization - E0102
PURPOSE:
	Plots for background removed, convolved and regridded snrs at all wavelengths.
	First plot shares the same color bar.
	Second plot they each have their own. 
"""
import sys
sys.path.insert(0, '../SED_Fitting/')
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt 

##############
## Switches ##
##############

save = True
show = False
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

xdim = [116,140]
ydim = [110,130]
title = ["24um","70um","100um","160um"]

##############
# Get Data   #
##############

# Get SNR data
data = [];
for i in range(4):

	data.append(fits.open(files[i])[0].data)

titles = ["24 $\mu$m","70 $\mu$m","100 $\mu$m","160 $\mu$m"] 

########################
# Plot Single Colorbar #
########################

# Single Colorbar
f, axes = plt.subplots(2,2,figsize=(5,5))

# Put all images at the same color range.
vmin = 0; vmax = 10
# Set x and y limits
x1 = 116; x2 = 140; y1 = 112; y2 = 128

count = 0
for i in range(2):
	for j in range(2):
		# Plot images, Set X,Y Limits, Set Titles
		cax = axes[i,j].imshow(fits.open(files[count])[0].data,vmin=vmin,vmax=vmax)
		axes[i,j].set_xlim(x1,x2); axes[i,j].set_ylim(y1,y2); axes[i,j].axis("off")
		axes[i,j].set_title(titles[count],fontsize="x-large")
		count += 1

# Find the left and right corner positions
p0 = axes[1,0].get_position().get_points().flatten()
p1 = axes[1,1].get_position().get_points().flatten()

# Add colorbar at the bottom
ax_cbar = f.add_axes([p0[0], 0.1, p1[2]-p0[0], 0.07])
cbar = plt.colorbar(cax, cax=ax_cbar, orientation='horizontal')
cbar.outline.set_visible(False)
cbar.ax.tick_params(labelsize=12) 
cbar.ax.set_title(r"MJy sr$^{-1}$")
plt.subplots_adjust(wspace=-.15,hspace=-.1,top=1,bottom=.120,right=1,left=0)

if save:
	plt.savefig("Plots/E0102_Regrid_Convolved.png",dpi=200)

######################
# Plot All Colorbars #
######################

fig, axs = plt.subplots(2,2,figsize=(6,5))
x1 = 116; x2 = 140; y1 = 112; y2 = 128
vmin = 0; vmax = 10

# Adjust the color bar
cbar_left = [0.41,0.9,0.4,0.9]
cbar_bottom = [0.54,0.54,0.06,0.06]


count = 0
for i in range(2):
	for j in range(2):
		# Plot image, set limits, plot titles. 
		im = axs[i,j].imshow(fits.open(files[count])[0].data,vmin=vmin,vmax=np.nanmax(fits.open(files[count])[0].data))
		axs[i,j].set_xlim(x1,x2); axs[i,j].set_ylim(y1,y2); axs[i,j].axis("off")
		axs[i,j].set_title(titles[count],fontsize="x-large")

		points = axs[i,j].get_position().get_points().flatten()
		# Adjust ColorBar: Left Bottom Width Height
		ax_cb = fig.add_axes([cbar_left[count],cbar_bottom[count],.035,.39])
		cb = plt.colorbar(im,cax=ax_cb)
		cb.ax.tick_params(labelsize=12) 
		cb.outline.set_visible(False)
		cb.ax.set_title(r"MJy sr$^{-1}$")
		count += 1
# Adjust Plot
plt.subplots_adjust(wspace=-0.05,hspace=0,top=.98,bottom=0,right=.95,left=-0.05)

if save:
	plt.savefig("Plots/E0102_Regrid_Convolved_Individual.png",dpi=200)

# Since we're saving plots we don't always need to open a GUI
if show:
	plt.show()