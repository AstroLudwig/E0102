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
from astropy.wcs import WCS
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

## For Testing 

# Plot pixel intensities in the original image
# as a function of pixel intensities in the modeled background.
# Fit a line to this data to use as the scale factor. 
plt_linearFit = False
# View region we use to adjust the scale factor.
plt_MedianAdjustRegionSelection = False
# Show the results of the SNR and the new modeled background
plt_SnrBkgd = False
# Show Histogram of intensities in entire SNR and the negative intensities
plt_Histogram = False
# Show Annulus in Image and Background before Subtraction
plt_FitRegion = False
# Show subtracted Annulus and histogram of values in it.
plt_subAnnulus = False
# Convolve and regrid 70 um modeled background
# to whichever file you're using.
# Doesn't need to be done more than once.
convolve_regrid = False

## For Code

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
	sname = '100um/100um_70modeled'

	# Arbitrary stuff for plotting
	vmin=0;vmax=30 #bkgd
	vmin_=-10; vmax_=100 #img
	vmin_sub = 4; vmax_sub= 30 # subtraction
	xdim = [170,230]; ydim = [160,220]

if im160:
	print("Image is 160 Microns")
	kern = '../../Convolve_Regrid/Kernels/Kernel_LoRes_PACS_70_to_PACS_160.fits'	
	im = '../../Sky_Remove/Median_Removed/160um_medianRemoved.fits'
	conv_sname = '../../Convolve_Regrid/Convolutions/70Bkgd_ext_20_tstep_3000_to_160umRes.fits'
	regrid_sname = '../../Convolve_Regrid/Convolutions/70Bkgd_ext_20_tstep_3000_to_160umRes_Resampled.fits'
	sigma  = np.loadtxt('../../Sky_Remove/Sigma.txt')[3]
	sname = '160um/160um_70modeled'
	
	# Arbitrary stuff for plotting
	vmin=5;vmax=50 #bkgd
	vmin_=10; vmax_=50 #img
	vmin_sub = 0; vmax_sub= 10 # subtraction
	xdim = [110,150]; ydim = [100,135]

if convolve_regrid: # Note that this gave me a lot of problems on windows but not on linux
	Convolve.master_convolve(kern,bkgd,conv_sname)
	Regrid.resample(conv_sname,im,regrid_sname)	

bkgd = fits.open(regrid_sname)[0].data
bkgdhdr = fits.open(regrid_sname)[0].header
img = fits.open(im)[0].data
hdr = fits.open(im)[0].header
w = WCS(hdr)
# Region that is commonly over subtracted
# Used to just adjust median slope to a value
# That makes this region have a median of 0. 
ra = 16.014763979322012 
dec = -72.02709660389223 
c = SkyCoord(ra=ra*u.degree,dec=dec*u.degree,frame='icrs')
print("Null Region Center Coordinate: "+str(c.ra.hms)+str(c.dec))
SubRadius = 3.; # ArcSeconds

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

#######################
# Adjust Scale Factor #
#######################

# We know there's a dark region of the SNR that is commonly oversubtracted in all images.
# This hand defines that region as sq_row sq_column and then adjusts the slope to best
# make that region 0. 

# Starting slope
slope = np.nanmedian(np.copy(imgShell)/np.copy(bkgdShell))

# Range of slopes based on starting slope
slopeRange = np.arange(slope-5,slope+5,.01)
# create new subs
testMedians = []

# Save some time by getting coordinates for the test region outside of the loop.
newBkgd = np.copy(bkgd) * slopeRange[0]
newSub = np.copy(img)- newBkgd

if im100:
	# This method works best in 100. 
	# Save coordinates to be transformed for 160.
	testRegion_r, testRegion_c = np.where(np.isfinite(prune.SelectRegionalData(newSub,hdr,ra,dec,SubRadius)))
	np.savetxt("NullRegionCoordinates.txt",w.all_pix2world(testRegion_c,testRegion_r,1))
if im160:
	# Transform coordinates used for null region in 100
	Ra,Dec = np.loadtxt("NullRegionCoordinates.txt")
	testRegion_c, testRegion_r = w.all_world2pix(Ra,Dec,1)
	testRegion_r = np.round(testRegion_r).astype(int); testRegion_c = np.round(testRegion_c).astype(int)

for i in range(len(slopeRange)):
	# Create new SNR based on an adjusted slope.
	newBkgd = np.copy(bkgd) * slopeRange[i]
	newSub = np.copy(img)- newBkgd
	# Get the region that is usually very negative.
	testRegion = newSub[testRegion_r,testRegion_c]
	# Take the median 
	testMedians.append(np.median(testRegion))

# Find which slope makes that median closest to 0 and then use that slope.	
find = np.where(np.isclose(testMedians,0.,atol=0.1))
newSlope = np.mean(slopeRange[find])
print(("Adjust slope by {}%").format(100*(slope-newSlope)/slope))
slope = newSlope
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
if plt_MedianAdjustRegionSelection:
	plt.figure()
	plt.imshow(snr)
	plt.ylim(ydim)
	plt.xlim(xdim)
	plt.scatter(testRegion_c,testRegion_r,c='red',s=0.1)
if plt_SnrBkgd:
	# Show the result
	# subtraction and background
	g, (ex,fx) = plt.subplots(1,2,sharey=True)
	ex.imshow(newsub,vmin=vmin_sub,vmax=vmax_sub)
	cfx=fx.imshow(newbkgd,vmin=vmin,vmax=vmax)
	ex.set_xlim(xdim); ex.set_ylim(ydim)
	fx.set_xlim(xdim); fx.set_ylim(ydim)
	ex.set_title("Subtracted")
	fx.set_title("Background")
	cbar = g.colorbar(cfx)
if plt_Histogram:
	# Show the pruned supernova remnant
	# Look at how many pixels are negative
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
	# Look at annulus in background and image before doing anything.
	p, (qx,rx) = plt.subplots(1,2)
	qx.imshow(bkgdShell,vmin=vmin,vmax=vmax)
	rx.imshow(imgShell,vmin=vmin_,vmax=vmax_)

	qx.set_xlim(xdim); qx.set_ylim(ydim)
	rx.set_xlim(xdim); rx.set_ylim(ydim)
	qx.set_title("Background Annulus")
	rx.set_title("Image Annulus")
if plt_subAnnulus:
	# Annulus in subtracted snr image
	f2, (ax2,bx2) = plt.subplots(1,2)
	subannulus = prune.Annulus(np.copy(newsub),hdr,thickness)
	ax2.imshow(subannulus)
	ax2.set_xlim(xdim); ax2.set_ylim(ydim)
	print(("Subtracted Annulus: Max {} Min {}").format(np.nanmax(subannulus),np.nanmin(subannulus)))
	bx2.hist(subannulus[np.where(np.isfinite(subannulus))].flatten())
	ax2.set_title("Linear Fit")
##############
# Save Game? #
##############

if save:
	fits.writeto(sname+'_snr.fits',newsub,hdr,overwrite=True)
	fits.writeto(sname+'_bkgd.fits',newbkgd,hdr,overwrite=True)
	print("Files saved.")
plt.show()