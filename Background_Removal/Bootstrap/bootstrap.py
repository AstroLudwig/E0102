# -*- Copyright (c) 2018, Bethany Ann Ludwig, All rights reserved. -*-
"""
NAME:
    Boot Strap
PURPOSE:
    Solve for a scale factor between 70 micron modeled background and both 100 and 160 micron images.
    Create a new modeled background and subtract from remnant.
"""
import sys
sys.path.insert(0, '../../Convolve_Regrid')

from astropy.io import fits
import numpy as np 
from astropy import wcs
import matplotlib.pyplot as plt 
from scipy import stats

import Convolve
import Regrid
from astropy.coordinates import SkyCoord
from astropy import units as u
import prune
####################
#     SWITCHES     #
####################

# Plot Switches For Linear Fit. 
plt_linearFit = True
plt_bkgdSubtraction = True
plt_Histogram = True 
plt_FitRegion = True
plt_subAnnulus = True  
plt_compareHist = True

# Convolve and regrid 70 um modeled background
# to whichever file you're using.
# Doesn't need to be done more than once.
convolve_regrid = True

# Fitting Switch
ForceSlope = False
clip = False # Changes the range you're fitting. 

# Annulus Size defined as 22 arcseconds from 
# E0102 center out to some thickness.
thickness = 40 #arcseconds

# Image Switches 
im100 = True
im160 = False

# Save File
save = False
##################################
# File Handling and Convolutions #
##################################
bkgd = '../Diffusion/im70/bkgd70/70um_diff_3000_steps_bkgd.fits'
print("Background Model is 70 Micron with Mask extended by 20 Pixels Diffused to Time Step: 3000.")

if im100:
	print("Image is 100 Microns")
	kern = '../../Convolve_Regrid/Kernels/Kernel_LoRes_PACS_70_to_PACS_100.fits'
	im =  '../../Sky_Remove/Median_Removed/100um_medianRemoved.fits'
	conv_sname = '../../Convolve_Regrid/Convolutions/70Bkgd_ext_20_tstep_3000_to_100umRes.fits'
	regrid_sname = '../../Convolve_Regrid/Convolutions/bkgd70_ext_20_tstep_3000_to_100umRes_Resampled.fits'

	vmin=0;vmax=30 #bkgd
	vmin_=-10; vmax_=100 #img
	vmin_sub = 4; vmax_sub= 30 # subtraction
	
	xdim = [170,230]; ydim = [160,220]
	fname = '100um_70modeled'
	sname = '100um_70modeled_ForceSlope'
	sq_row = [200,200,200,201,201,201,202,202,202]
	sq_col = [198,199,200,198,199,200,198,199,200]

if im160:
	print("Image is 160 Microns")
	kern = '../../Convolve_Regrid/Kernels/Kernel_LoRes_PACS_70_to_PACS_160.fits'	
	im = '../../Sky_Remove/Median_Removed/160um_medianRemoved.fits'
	conv_sname = '../../Convolve_Regrid/Convolutions/70Bkgd_ext_20_tstep_3000_to_160umRes.fits'
	regrid_sname = '../../Convolve_Regrid/Convolutions/bkgd70_ext_20_tstep_3000_to_160umRes_Resampled.fits'

	vmin=5;vmax=50 #bkgd
	vmin_=10; vmax_=50 #img
	vmin_sub = 4; vmax_sub= 30 # subtraction

	xdim = [110,150]; ydim = [100,135]
	fname = '160um_70modeled'
	sname = '160um_70modeled_ForceSlope'
	sq_row = [124,124,124,125,125,125,126,126,126]
	sq_col = [124,125,126,124,125,126,124,125,126]

if convolve_regrid:
	#Convolve.master_convolve(kern,bkgd,conv_sname)
	Regrid.resample(conv_sname,im,regrid_sname)	

bkgd = fits.open(regrid_sname)[0].data
img = fits.open(im)[0].data
hdr = fits.open(im)[0].header
w = wcs.WCS(hdr)


##################
# Fitting Region #
##################

# Create Annulus of data in both the image and the background.

imgShell = prune.Annulus(np.copy(img),"Annulus")		
bkgdShell = prune.Annulus(np.copy(bkgd),"Annulus")

Ratio = np.copy(bkgdShell)/np.copy(imgShell)

###############
#  Clipping   #
###############

def clip(img,ratio,level):
	nData = np.copy(img)
	ind_r, ind_c = np.shape(ratio)
	med = np.nanmedian(ratio); std = np.nanstd(ratio)
	print(("Fitting Region. Median: {} Std: {}").format(med,std))

	clipRange = [med - level*std, med  + level*std]
	print(("Clipping Range: {} to {}").format(clipRange[0],clipRange[1]))				

	for i in range(ind_r):
		for j in range(ind_c):
			if ratio[i,j] < clipRange[0]: #or clipImg[i,j] > clipRange[1]:
				nData[i,j] = np.nan
			elif ratio[i,j] > clipRange[1]:
				nData[i,j] = np.nan
	return nData 

if clip:
	level = 2	# Select  Level of Clipping

	imgShell = clip(imgShell,Ratio,level)
	bkgdShell = clip(bkgdShell,Ratio,level)

	newRatio = np.copy(bkgdShell)/np.copy(imgShell)

	med = np.nanmedian(newRatio)

	afterClip = np.copy(bkgdShell/imgShell)

in_img = img[np.where(np.isfinite(imgShell))]
in_bkgd = bkgd[np.where(np.isfinite(imgShell))]

in_bkgd = in_bkgd.flatten()
in_img = in_img.flatten()

#######################################################
# Fitting elements in Annulus with Linear Regression  #
#######################################################

# Create x values
xvals = np.arange(np.nanmin(in_bkgd),np.nanmax(in_bkgd))

# Uses the median as a slope.
line = xvals * 1/med

if ForceSlope:
	# Select a region in img 

	# Starting place 
	med = np.nanmedian(np.copy(bkgdShell)/np.copy(imgShell))
	# Create a starting image
	line = xvals * 1/med
	newbkgd = np.copy(bkgd) * 1/med
	
	newsub = np.copy(img) - newbkgd
	# This is the region we want to be 0
	sq = newsub[sq_row,sq_col]	
	# range of other slopes to test things out on 
	medRange = np.arange(med-5,med+5,.01)
	# create new subs
	testMedians = []

	for i in range(len(medRange)):

		newBkgd = np.copy(bkgd) * 1/medRange[i]
		newSub = np.copy(img)-newBkgd
		testRegion = newSub[sq_row,sq_col]
		testMedians.append(np.median(testRegion))

	find = np.where(np.isclose(testMedians,0.,atol=0.1))
	print(("New Slope: {} ").format(medRange[find]))

	med = medRange[find]

######################
# Create Subtraction #
######################

newbkgd = np.copy(bkgd) * 1/med
newsub = np.copy(img) - newbkgd
newsubSNR = prune(np.copy(newsub),"snr")
plt.figure()

plt.imshow(newsubSNR,vmin = np.nanmin(newsubSNR),vmax=np.nanmax(newsubSNR))
plt.xlim(xdim); plt.ylim(ydim)
plt.scatter(sq_col,sq_row,c='r',s=0.5)
histsub = prune(np.copy(newsub),"snr")
histdata = histsub[np.where(np.isfinite(histsub))].flatten()

subSum = prune(newsub,"snr")
print(("Sum of Subtraction: {}").format(np.nansum(subSum)))	
subSnr = subSum[np.where(np.isfinite(newsub))]
print(("Sum of Subtraction Negatives: {}").format(np.nansum(subSnr[np.where(subSnr<0)])))	

########################
# Plotting and Testing #
########################
# Check that everythings imported and fine. 

if plt_linearFit:
	plt.figure()
	plt.scatter(in_bkgd,in_img,s=.1)
	plt.plot(xvals,line)
	plt.xlabel('Background Pixel Intensities')
	plt.ylabel('Original Image Pixel Intensities')
	plt.title("Annulus Intensity Ratio")
if plt_bkgdSubtraction:
	g, (ex,fx) = plt.subplots(1,2,sharey=True)
	#ex.imshow(newsub,vmin=vmin_sub,vmax=vmax_sub)
	ex.imshow(newsub,vmin=vmin,vmax=vmax)
	cfx=fx.imshow(newbkgd,vmin=vmin,vmax=vmax)
	#cfx=fx.imshow(newbkgd,vmin=vmin_sub,vmax=vmax_sub)
	ex.set_xlim(xdim); ex.set_ylim(ydim)
	fx.set_xlim(xdim); fx.set_ylim(ydim)
	ex.set_title("Subtracted")
	fx.set_title("Background")
	cbar = g.colorbar(cfx)
if plt_Histogram:
	k, (lx,kx) = plt.subplots(1,2)
	lx.hist(histdata)
	kx.imshow(histsub)
	kx.set_xlim(xdim); kx.set_ylim(ydim)
	kx.set_title("Subtraction. Outer Radius Removed")
	lx.set_title("Pixel Intensities in Subtracted Remnant")
	#negs
	plt.figure()
	a = np.where(histdata<0)
	plt.hist(histdata[a])
if plt_FitRegion:
#	x,y = w.all_world2pix(ra,dec,1)
	p, (qx,rx) = plt.subplots(1,2)
	qx.imshow(bkgdShell,vmin=vmin,vmax=vmax)
	qx.scatter(x,y,s=.1)
	rx.imshow(imgShell,vmin=vmin_,vmax=vmax_)
	rx.scatter(x,y,s=.1)
	qx.set_xlim(xdim); qx.set_ylim(ydim)
	rx.set_xlim(xdim); rx.set_ylim(ydim)
	qx.set_title("Background Annulus")
	rx.set_title("Image Annulus")
if plt_subAnnulus:
	f2, (ax2,bx2) = plt.subplots(1,2)
	subannulus = prune(np.copy(newsub),"Annulus")
	if clip:
		subannulus = clip(subannulus)
	ax2.imshow(subannulus)
	ax2.set_xlim(xdim); ax2.set_ylim(ydim)
	print(("Subtracted Annulus: Max {} Min {}").format(np.nanmax(subannulus),np.nanmin(subannulus)))
	bx2.hist(subannulus[np.where(np.isfinite(subannulus))].flatten())
	ax2.set_title("Linear Fit")
if plt_compareHist:
	plt.figure(7)
	beforeClipHist = beforeClip[np.where(np.isfinite(beforeClip))]
	afterClipHist = afterClip[np.where(np.isfinite(afterClip))]
	plt.hist(beforeClipHist,bins=10,color='blue')
	plt.hist(afterClipHist,bins=10,color='red')

##############
# Save Game? #
##############

if save:

	fits.writeto(fname+'/'+sname+'_snr.fits',newsub,hdr,overwrite=True)
	fits.writeto(fname+'/'+sname+'_bkgd.fits',newbkgd,hdr,overwrite=True)

plt.show()