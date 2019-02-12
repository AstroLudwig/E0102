# -*- Copyright (c) 2019, Bethany Ann Ludwig, All rights reserved. -*-
"""
NAME:
	Data Visualization - Pixel by Pixel SED
PURPOSE:
	Plot of a single data point, 160 microns, 
	with 3 sample seds from two cold "hotspots" and the hotspot at 24 microns.
"""
import sys
sys.path.insert(0, '../SED_Fitting/')

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt 
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
import Eqs

##############
## Switches ##
##############
version_1 = False
version_2 = True
save = True
show = False
sigma_fit = 0

##################
# Data Retrieval #
##################
files = ['../Final_Files/24/24_SNR_Convolve_Regrid_Prune.fits','../Final_Files/70/70_SNR_Convolve_Regrid_Prune.fits',
		'../Final_Files/100/100_SNR_Convolve_Regrid_Prune.fits','../Final_Files/160/160_SNR_Prune.fits']
# http://docs.astropy.org/en/stable/constants/
# Using CGS units

xdim = [120,137]
ydim = [112,128]
title = ["24um","70um","100um","160um"]

lam, k_amc = Eqs.log_extrap('../SED_Fitting/Kappa/kappa_amc.dat')
wv = [24,70,100,160]

# Get SNR data, the physical area of each SNR, and the average pixel intensity
data = []; 
for i in range(4):
	data.append(fits.open(files[i])[0].data)
# Which image to plot 
img = data[3]
# Get Pixel by Pixel Solutions 
Pix_Warm_SED =  np.load("../SED_Fitting/Sols/PixbyPix/Warm_SED.npy")
Pix_Cold_SED =  np.load("../SED_Fitting/Sols/PixbyPix/Cold_SED.npy")

# Select pixels to view 
Px = [133,123,128]; Py = [121,118,116]
##########
# Plot   #
##########
if version_1:
	# Initialize Plot
	f = plt.figure(figsize=(6,7))

	# Initialize SED Plots
	cx = plt.subplot2grid((3,2),(0,1))
	dx = plt.subplot2grid((3,2),(1,1))
	ix = plt.subplot2grid((3,2),(2,1))

	# Initialize 160 um image. 
	ex = plt.subplot2grid((3,2),(0,0),rowspan=3)

	plots = [cx,dx,ix]
	#colors = ["red","aqua","yellow"]
	colors = ["darkviolet","#FF00FF","red"]
	cmap = plt.get_cmap('viridis')

	# SNR Image
	ex.imshow(img,vmin=np.nanmin(img),vmax=np.nanmax(img))
	ex.set_xlim(xdim)
	ex.set_ylim(ydim)
	ex.axis("off")
	ex.scatter(Px,Py,s=60, marker='s',facecolors='none',edgecolors=colors)

	for i in range(3):
		ObsvSed = [data[0][Py[i],Px[i]],data[1][Py[i],Px[i]],data[2][Py[i],Px[i]],data[3][Py[i],Px[i]]]
		ObsvErr = Eqs.error(ObsvSed)
		plots[i].plot(lam,Pix_Warm_SED[Py[i],Px[i]]+Pix_Cold_SED[Py[i],Px[i]],label="Total SED",color="#424186")
		plots[i].plot(lam,Pix_Warm_SED[Py[i],Px[i]],label="Warm SED",color="#84D44B",ls='dashed')
		plots[i].plot(lam,Pix_Cold_SED[Py[i],Px[i]],label="Cold SED",color="#23A883")		

		plots[i].errorbar(wv,ObsvSed,yerr=ObsvErr,marker='o',linestyle='none',c='black')
		plots[i].grid(color='white',linestyle='-')
		plots[i].set_facecolor("#EAEAF2")

		plots[i].spines['bottom'].set_color(colors[i])
		plots[i].spines['top'].set_color(colors[i]) 
		plots[i].spines['right'].set_color(colors[i])
		plots[i].spines['left'].set_color(colors[i])

	plots[0].legend()	
	plots[1].set_ylabel("Spectral Intensity (MJy sr$^{-1}$)")
	plots[2].set_xlabel("Wavelength ($\mu m$)")

	plt.subplots_adjust(top=.99, bottom=0.05, left=0.0, right=.99, hspace=0.2, wspace=0.09)
if version_2:
	# Initialize Plot
	f = plt.figure(figsize=(7,6))

	# Initialize SED Plots
	cx = plt.subplot2grid((2,2),(0,1))
	dx = plt.subplot2grid((2,2),(0,0))
	ix = plt.subplot2grid((2,2),(1,1))

	# Initialize 160 um image. 
	ex = plt.subplot2grid((2,2),(1,0),rowspan=1)

	plots = [cx,dx,ix]
	#colors = ["red","aqua","yellow"]
	colors = ["darkviolet","#FF00FF","red"]
	cmap = plt.get_cmap('viridis')

	# SNR Image
	ex.imshow(img,vmin=np.nanmin(img),vmax=np.nanmax(img))
	ex.set_xlim(xdim)
	ex.set_ylim(ydim)
	ex.axis("off")
	ex.scatter(Px,Py,s=80, marker='s',facecolors='none',edgecolors=colors)

	for i in range(3):
		ObsvSed = [data[0][Py[i],Px[i]],data[1][Py[i],Px[i]],data[2][Py[i],Px[i]],data[3][Py[i],Px[i]]]
		ObsvErr = Eqs.error(ObsvSed)
		plots[i].plot(lam,Pix_Warm_SED[Py[i],Px[i]]+Pix_Cold_SED[Py[i],Px[i]],label="Total SED",color="#424186")
		plots[i].plot(lam,Pix_Warm_SED[Py[i],Px[i]],label="Warm SED",color="#84D44B",ls='dashed')
		plots[i].plot(lam,Pix_Cold_SED[Py[i],Px[i]],label="Cold SED",color="#23A883")		

		plots[i].errorbar(wv,ObsvSed,yerr=ObsvErr,marker='o',linestyle='none',c='black')
		plots[i].grid(color='white',linestyle='-')
		plots[i].set_facecolor("#EAEAF2")

		plots[i].spines['bottom'].set_color(colors[i])
		plots[i].spines['top'].set_color(colors[i]) 
		plots[i].spines['right'].set_color(colors[i])
		plots[i].spines['left'].set_color(colors[i])

	plots[0].legend()	
	plots[1].set_ylabel("Spectral Intensity (MJy sr$^{-1}$)")
	plots[2].set_xlabel("Wavelength ($\mu m$)")

	plt.subplots_adjust(top=.99, bottom=0.12, left=0.09, right=.99, hspace=0.2, wspace=0.12)
if save:
	plt.savefig("Plots/PixelByPixel.png",dpi=200)
if show:
	plt.show()
