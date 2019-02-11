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

save = True
# ~~~ Plotting ~~~ #
# Plot 4 Background Removed and Convolved/Regridded to 160 images.
plot_imgs = False


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

##################
# Plotting Range #
##################

lam, k_amc = Eqs.log_extrap('../SED_Fitting/Kappa/kappa_amc.dat')
lam, k_mg2 = Eqs.log_extrap('../SED_Fitting/Kappa/kappa_mg2sio4.dat')
wv = [24,70,100,160]
pix_area = Eqs.pixarea(files[0])

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


##########
# Plot   #
##########

# Single Colorbar
f, axes = plt.subplots(2,2,figsize=(10,15))
x1 = 116; x2 = 140; y1 = 112; y2 = 128
vmin = 0; vmax = 10
axes[0,0].imshow(fits.open(files[0])[0].data,vmin=vmin,vmax=vmax)
axes[0,1].imshow(fits.open(files[1])[0].data,vmin=vmin,vmax=vmax)
axes[1,0].imshow(fits.open(files[2])[0].data,vmin=vmin,vmax=vmax)
cax = axes[1,1].imshow(fits.open(files[3])[0].data,vmin=vmin,vmax=vmax)
axes[0,0].set_xlim(x1,x2); axes[0,0].set_ylim(y1,y2); axes[0,0].axis("off")
axes[0,1].set_xlim(x1,x2); axes[0,1].set_ylim(y1,y2); axes[0,1].axis("off")
axes[1,0].set_xlim(x1,x2); axes[1,0].set_ylim(y1,y2); axes[1,0].axis("off")
axes[1,1].set_xlim(x1,x2); axes[1,1].set_ylim(y1,y2); axes[1,1].axis("off") 
 
axes[0,0].set_title("24 $\mu$m",fontsize="x-large"); axes[0,1].set_title("70 $\mu$m",fontsize="x-large")
axes[1,0].set_title("100 $\mu$m",fontsize="x-large"); axes[1,1].set_title("160 $\mu$m",fontsize="x-large")
f.set_size_inches(5, 5)


p0 = axes[1,0].get_position().get_points().flatten()
p1 = axes[1,1].get_position().get_points().flatten()

ax_cbar = f.add_axes([p0[0], 0.1, p1[2]-p0[0], 0.07])
cbar = plt.colorbar(cax, cax=ax_cbar, orientation='horizontal')
cbar.outline.set_visible(False)
cbar.ax.tick_params(labelsize=12) 
plt.subplots_adjust(wspace=-.15,hspace=-.1,top=1,bottom=.120,right=1,left=0)

if save:
	plt.savefig("Plots/E0102_Regrid_Convolved.png",dpi=200, bbox_inches='tight', pad_inches = 0 )

# Individual Colorbar 

fig, axs = plt.subplots(2,2,figsize=(10,15))
x1 = 116; x2 = 140; y1 = 112; y2 = 128
vmin = 0; vmax = 10

titles = ["24 $\mu$m","70 $\mu$m","100 $\mu$m","160 $\mu$m"] 
cbar_pad = [0.47,0.93,0.47,0.93]
#fig.set_size_inches(5, 5)
#ims = [im1,im2,im3,im4]

count = 0
for i in range(2):
	for j in range(2):
		im = axs[i,j].imshow(fits.open(files[count])[0].data,vmin=vmin,vmax=np.nanmax(fits.open(files[count])[0].data))
		axs[i,j].set_xlim(x1,x2); axs[i,j].set_ylim(y1,y2); axs[i,j].axis("off")
		axs[i,j].set_title(titles[count],fontsize="x-large")

		points = axs[i,j].get_position().get_points().flatten()
		ax_cb = fig.add_axes([cbar_pad[count],points[1]+0.07,.035,.35])
		print(points)
		cb = plt.colorbar(im,cax=ax_cb)
		cb.ax.tick_params(labelsize=12) 
		cb.outline.set_visible(False)
		count += 1
plt.subplots_adjust(wspace=-.15,hspace=-.1,top=1,bottom=.120,right=1,left=0)

if save:
	plt.savefig("Plots/E0102_Regrid_Convolved_Individual.png",dpi=200, bbox_inches='tight', pad_inches = 0 )

plt.show()