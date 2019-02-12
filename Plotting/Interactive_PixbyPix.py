# -*- Copyright (c) 2019, Bethany Ann Ludwig, All rights reserved. -*-
"""
NAME:
	Data Visualization - Interactive Pixel by Pixel Plot
PURPOSE:
	Click and plot an SED
"""
import sys
sys.path.insert(0, '../SED_Fitting/')
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt 
import Eqs

###########
## Files ##
###########

files = ['../Final_Files/24/24_SNR_Convolve_Regrid_Prune.fits','../Final_Files/70/70_SNR_Convolve_Regrid_Prune.fits',
		'../Final_Files/100/100_SNR_Convolve_Regrid_Prune.fits','../Final_Files/160/160_SNR_Prune.fits']
data = []; 
for i in range(4):
	data.append(fits.open(files[i])[0].data)

#############
# Constants #
#############

xdim = [116,140]; ydim = [110,130]
wv = [24,70,100,160]

title = ["24um","70um","100um","160um"]
lam, k_amc = Eqs.log_extrap('../SED_Fitting/Kappa/kappa_amc.dat')


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

##########
# Plot   #
##########
f, (ax,bx) = plt.subplots(1,2)

inObsvSed = [data[0][0,0],data[1][0,0],data[2][0,0],data[3][0,0]]
inObsvErr = Eqs.error(inObsvSed)

ax.imshow(data[3])
ax.scatter(0,0,c='r',marker='s',s=40)
bx.plot(lam,Pix_Warm_SED[0,0]+Pix_Cold_SED[0,0],label="Total SED")
bx.plot(lam,Pix_Warm_SED[0,0],label="Warm SED")
bx.plot(lam,Pix_Cold_SED[0,0],label="Cold SED")
bx.errorbar(wv,inObsvSed,yerr=inObsvErr,marker='o',linestyle='none',c='purple')

ax.set_xlim(xdim[0],xdim[1])
ax.set_ylim(ydim[0],ydim[1])
bx.set_ylim(0, np.max(Pix_Warm_SED[0,0]+Pix_Cold_SED[0,0])+2)
bx.grid()

plt.legend()
plt.title(("Temperature: {:2d} K,  Cold Mass: {:.2f} $M_\odot$,  Warm Mass: {:.2E} $M_\odot$,  $\chi^2$: {:.2f}").format(int(Pix_Temp[0,0]),Pix_Cold_Mass[0,0],Pix_Warm_Mass[0,0],Pix_chisqrd[0,0]))
bx.set_xlabel("Wavelength ($\mu$m)")
bx.set_ylabel("Spectral Intensity (MJy sr$^{-1}$)")
ax.axis("off")

def onclick(event):
	i = int(event.ydata)
	j = int(event.xdata)

	ObsvSed = [data[0][i,j],data[1][i,j],data[2][i,j],data[3][i,j]]
	ObsvErr = Eqs.error(ObsvSed)        
	ax.clear()
	bx.clear()

	ax.imshow(data[3])
	ax.scatter(j,i,c='r',marker='s',s=40)
	ax.set_xlim(xdim[0],xdim[1])
	ax.set_ylim(ydim[0],ydim[1])
	ax.axis("off")

	bx.plot(lam,Pix_Warm_SED[i,j]+Pix_Cold_SED[i,j],label="Total SED")
	bx.plot(lam,Pix_Warm_SED[i,j],label="Warm SED")
	bx.plot(lam,Pix_Cold_SED[i,j],label="Cold SED")
	
	bx.set_ylim(0, np.max(Pix_Warm_SED[i,j]+Pix_Cold_SED[i,j])+2)
	bx.errorbar(wv,ObsvSed,yerr=ObsvErr,marker='o',linestyle='none',c='purple')
	bx.grid()
	plt.legend()
	plt.title(("Temperature: {:2d} K,  Cold Mass: {:.2E} $M_\odot$,  Warm Mass: {:.2E} $M_\odot$,  $\chi^2$: {:.2f}").format(int(Pix_Temp[i,j]),Pix_Cold_Mass[i,j],Pix_Warm_Mass[i,j],Pix_chisqrd[i,j]))
	bx.set_xlabel("Wavelength ($\mu$m)")
	bx.set_ylabel("Spectral Intensity (MJy sr$^{-1}$)")

	plt.draw()

cid = f.canvas.mpl_connect('button_press_event', onclick)

plt.show()
