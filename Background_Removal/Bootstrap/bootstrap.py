# -*- Copyright (c) 2018, Bethany Ann Ludwig, All rights reserved. -*-
"""
NAME:
    Boot Strap
PURPOSE:
    Solve for a scale factor between 70 micron modeled background and both 100 and 160 micron images.
    Create a new modeled background and subtract from remnant.
KNOWN PROBLEM:
	Welcome to the joys of cross platform development. reproject_exact blew up my windows computer but 
	appears to run fine on ubuntu. If getting an insane infite loop when running this code then regrid 
	is being a problem. 
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

# Plot pixel intensities in the original image
# as a function of pixel intensities in the modeled background.
# Fit a line to this data to use as the scale factor. 
plt_linearFit = True
# Show the results of the SNR and the new modeled background
plt_result = True
plt_Histogram = False
plt_FitRegion = False
plt_subAnnulus = True
plt_compareHist = False
plt_sloppy = False

# Convolve and regrid 70 um modeled background
# to whichever file you're using.
# Doesn't need to be done more than once.
convolve_regrid = False

# Fitting Switch
Detection_Limited = True
OverSubtracted_Region = False
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
	regrid_sname = '../../Convolve_Regrid/Convolutions/70Bkgd_ext_20_tstep_3000_to_100umRes_Resampled.fits'
	sigma  = np.loadtxt('../../Sky_Remove/Sigma.txt')[2]
	fname = '100um_70modeled'
	sname = '100um_70modeled_ForceSlope'

	# Arbitrary stuff for plotting
	vmin=0;vmax=30 #bkgd
	vmin_=-10; vmax_=100 #img
	vmin_sub = 4; vmax_sub= 30 # subtraction
	xdim = [170,230]; ydim = [160,220]

	# Region that is commonly over subtracted
	sq_row = [200,200,200,201,201,201,202,202,202]
	sq_col = [198,199,200,198,199,200,198,199,200]

if im160:
	print("Image is 160 Microns")
	kern = '../../Convolve_Regrid/Kernels/Kernel_LoRes_PACS_70_to_PACS_160.fits'	
	im = '../../Sky_Remove/Median_Removed/160um_medianRemoved.fits'
	conv_sname = '../../Convolve_Regrid/Convolutions/70Bkgd_ext_20_tstep_3000_to_160umRes.fits'
	regrid_sname = '../../Convolve_Regrid/Convolutions/70Bkgd_ext_20_tstep_3000_to_160umRes_Resampled.fits'
	sigma  = np.loadtxt('../../Sky_Remove/Sigma.txt')[3]
	fname = '160um_70modeled'
	sname = '160um_70modeled_ForceSlope'
	
	# Arbitrary stuff for plotting
	vmin=5;vmax=50 #bkgd
	vmin_=10; vmax_=50 #img
	vmin_sub = 4; vmax_sub= 30 # subtraction
	xdim = [110,150]; ydim = [100,135]
	
	# Region that is commonly over subtracted
	sq_row = [124,124,124,125,125,125,126,126,126]
	sq_col = [124,125,126,124,125,126,124,125,126]

if convolve_regrid: # Note that this gave me a lot of problems on windows but not on linux
	Convolve.master_convolve(kern,bkgd,conv_sname)
	Regrid.resample(conv_sname,im,regrid_sname)	

bkgd = fits.open(regrid_sname)[0].data
bkgdhdr = fits.open(regrid_sname)[0].header
img = fits.open(im)[0].data
hdr = fits.open(im)[0].header

##################
# Fitting Region #
##################
# Create Annulus of data in both the image and the background.

imgShell = prune.Annulus(np.copy(img),hdr,thickness)		
bkgdShell = prune.Annulus(np.copy(bkgd),bkgdhdr,thickness)

Ratio = np.copy(bkgdShell)/np.copy(imgShell)
med = np.nanmedian(Ratio)

# Get 1D arrays of pixel intensities. 
in_img = img[np.where(np.isfinite(imgShell))].flatten()
in_bkgd = bkgd[np.where(np.isfinite(imgShell))].flatten()

########################
#  Detection Limited   #
########################
if Detection_Limited:
	def FindLimitedPixels(data,level): 
		# Returns sensitivity limited data, anything below X sigma is included.
		# in the original image we are modeling
		# Level is to adjust sensitivity (I.e. 2 sigma, 3 sigma, etc.)
		row = []; column = [];
		# We only care about pixels within the SNR
		nImg = prune.Prune(data,hdr)
		for i in range(np.shape(nImg)[0]):
			for j in range(np.shape(nImg)[1]):
				if np.isfinite(nImg[i,j]) and nImg[i,j] < sigma * level:
					row.append(i); column.append(j)			
		return data[row, column]
	# We want to pick a slope that doesn't over compensate
	# for extreme intensity values due to the nebula.

	# Our test for a slope, will be one that gets the detection 
	# limited pixels closest to zero.

	# Get an initial slope without adjusting for sensitivity limited pixels.
	slope = 1/np.nanmedian(np.copy(bkgdShell)/np.copy(imgShell))
	xvals = np.arange(np.nanmin(in_bkgd),np.nanmax(in_bkgd))
	line = xvals * slope
	print("The initial slope based on the inverse median is: "+str(slope))
	# Create a range of those slopes
	slopeRange = np.arange(slope-5,slope+5,.1)

	# Find the subtraction that gets the detection limited pixels closest to 0.
	testResults = []; 
	for i in range(len(slopeRange)):
		newBkgd = np.copy(bkgd) * slopeRange[i]
		newSnr = np.copy(img) - newBkgd
		testResults.append(np.median(FindLimitedPixels(newSnr,2)))
	find = np.where(np.isclose(testResults,0.,atol=0.1))	
	slope = slopeRange[find]
	line = xvals * slope
	print("Adjusting slope to account for intensity extremes. New slope is: "+str(slope))

########################
#  Detection Limited   #
########################
if OverSubtracted_Region:
	# We know there's a dark region of the SNR that is commonly oversubtracted in all images.
	# This hand defines that region as sq_row sq_column and then adjusts the slope to best
	# make that region 0. 

	# Starting slope
	slope = 1/np.nanmedian(np.copy(bkgdShell)/np.copy(imgShell))

	# Range of slopes based on starting slope
	slopeRange = np.arange(slope-5,slope+5,.01)
	# create new subs
	testMedians = []

	for i in range(len(slopeRange)):
		# Create new SNR based on an adjusted slope.
		newBkgd = np.copy(bkgd) * slopeRange[i]
		newSub = np.copy(img)- newBkgd
		# Get the region that is usually very negative.
		testRegion = newSub[sq_row,sq_col]
		# Take the median 
		testMedians.append(np.median(testRegion))

	# Find which slope makes that median closest to 0 and then use that slope.	
	find = np.where(np.isclose(testMedians,0.,atol=0.1))
	
	slope = np.mean(slopeRange[find])
	print(("New Slope: {} ").format(slope))
	
	xvals = np.arange(np.nanmin(in_bkgd),np.nanmax(in_bkgd))
	line = xvals * slope

################################
#  Create Background and SNR   #
################################

newbkgd = np.copy(bkgd) * slope
newsub = np.copy(img) - newbkgd
snr = prune.Prune(np.copy(newsub),hdr)


# Create data for histogram
histData = snr[np.where(np.isfinite(snr))].flatten()

# Print some results to quantify oversubtraction
print(("Supernova Total: {}").format(np.nansum(snr)))	
print(("Oversubtraction total: {}").format(np.nansum(histData[np.where(histData<0)])))	

########################
# Plotting and Testing #
########################

if plt_linearFit:
	# Plot image as a function of background.
	# Fit a line using median as slope
	plt.figure()
	plt.scatter(in_bkgd,in_img,s=.1)
	plt.plot(xvals,line)
	plt.xlabel('Background Pixel Intensities')
	plt.ylabel('Original Image Pixel Intensities')
	plt.title("Annulus Intensity Ratio")
if plt_result:
	# Show the result
	# subtraction and background
	g, (ex,fx) = plt.subplots(1,2,sharey=True)
	ex.imshow(newsub,vmin=vmin,vmax=vmax)
	cfx=fx.imshow(newbkgd,vmin=vmin,vmax=vmax)
	ex.set_xlim(xdim); ex.set_ylim(ydim)
	fx.set_xlim(xdim); fx.set_ylim(ydim)
	ex.set_title("Subtracted")
	fx.set_title("Background")
	cbar = g.colorbar(cfx)
if plt_Histogram:
	# Show the pruned supernova remnant
	# Look at how much of the pixels are negative
	# Look closer at how the negatives are distributed
	k, (lx,kx,jx) = plt.subplots(1,3)
	lx.hist(histData)
	kx.imshow(snr)
	jx.hist(histData[np.where(histData<0)])
	kx.set_xlim(xdim); kx.set_ylim(ydim)
	kx.set_title("Pruned Subtraction. ")
	lx.set_title("Pixel Intensities in Subtracted Remnant")
	jx.set_title("Negative Pixel Values Only")
	
if plt_FitRegion:
	# Look at annulus in background and image.
	# Take clipping into account.
	p, (qx,rx) = plt.subplots(1,2)
	qx.imshow(bkgdShell,vmin=vmin,vmax=vmax)
	rx.imshow(imgShell,vmin=vmin_,vmax=vmax_)

	qx.set_xlim(xdim); qx.set_ylim(ydim)
	rx.set_xlim(xdim); rx.set_ylim(ydim)
	qx.set_title("Background Annulus")
	rx.set_title("Image Annulus")
if plt_subAnnulus:
	f2, (ax2,bx2) = plt.subplots(1,2)
	subannulus = prune.Annulus(np.copy(newsub),hdr,thickness)
	#if clip:
	#	subannulus = clip(subannulus,Ratio,level)
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