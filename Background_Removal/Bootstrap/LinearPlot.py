from astropy.io import fits
import numpy as np 
from astropy import wcs
import matplotlib.pyplot as plt 
from scipy import stats
import sys
sys.path.insert(0, '../../Convolve_Regrid')
import Convolve
import regrid
from astropy.coordinates import SkyCoord
from astropy import units as u
# Plot Switches For Linear Fit. 
# Turn linear fit on before using these
plt_checkMask = False
plt_linearFit = True
plt_bkgdSubtraction = False
plt_Histogram = True # 
plt_FitRegion = False
plt_subAnnulus = True# 
plt_compareHist = False
# Convolve Switch. 
# Use for each new file.
convolve_regrid = False
# Fitting Switch
slopeInt = False
NoOffset = False # Uses the median as a slope. 
ForceSlope = True
clip = False # Changes the range you're fitting. 
# Need a file switch and image switch turned on at the same time.
# File Switches
image70 = False
image70_ext_20 = True
image100 = False
# Image Switches 
im70_to_100 = True
im70_to_160 = False
im100_to_160 = False
# Save File
save = False
# How far out you want the fitting region to extend. Doesnt have anything to do with the mask. 
# Has to do with how far out the annulus goes. 
extend = 40 #arcseconds
# Which time step file you want to use.
step = 3000

##################################
# File Handling and Convolutions #
##################################
"""
If you want to add another file. You'll need to turn all the switches off 
except convolve and the image switch and run that before using the file.
Then turn convolve off and everything else should work.
"""
# Select which 70 file to use.
if image70: 
	# For Convolution
	bkgd_conv = '../Diffusion/im70/bkgd70/70um_diff_'+str(step)+'_steps_bkgd.fits'
	conv_sname = '../../Convolve_Regrid/Convolutions/bkgd70_ext_0_tstep_'+str(step)+'_to_'
	regrid_sname = '../../Convolve_Regrid/Conv_Regrids/bkgd70_ext_0_tstep_'+str(step)+'_to_'
	# After Convolution
	if convolve_regrid == False:
		bkgd_70to100 = fits.open(regrid_sname+'100_CR.fits')[0].data
		bkgd_70to160 = fits.open(regrid_sname+'160_CR.fits')[0].data	
if image70_ext_20: # At the 1740th step. Can change that later.
	print("Background Model is 70 micron with mask extended by 20 pixels.")
	# For Convolution
	bkgd_conv = '../Diffusion/im70_ext_20/bkgd70/70um_diff_'+str(step)+'_steps_bkgd.fits'
	
	conv_sname = '../../Convolve_Regrid/Convolutions/bkgd70_ext_20_tstep_'+str(step)+'_to_'
	
	regrid_sname = '../../Convolve_Regrid/Conv_Regrids/bkgd70_ext_20_tstep_'+str(step)+'_to_'
	# I just did these below because I wanted to test something with the original sky removed
	# 70 micron img convolved and regridded to 100. Sorry for making a little mess here.
	#bkgd_conv = '../../Sky_Remove/Median_Removed/70um_medianRemoved.fits'
	#conv_sname = '../../Convolve_Regrid/Convolutions/SkyRemovd70_to_'
	#regrid_sname = '../../Convolve_Regrid/Conv_Regrids/SkyRemovd70_to_'
	if convolve_regrid == False:
		bkgd_70to100 = fits.open(regrid_sname+'100_CR.fits')[0].data
		bkgd_70to160 = fits.open(regrid_sname+'160_CR.fits')[0].data		
if image100: 
	print("Background Model is 100 micron bootstrapped by 70 with mask extended by 20 pixels.")
	# For Convolution
	bkgd_conv = '../Bootstrap/100um_70modeled/100um_70modeled_bkgd.fits'
	conv_sname = '../../Convolve_Regrid/Convolutions/bkgd100_modled_by_70_ext_20_tstep_'+str(step)+'to_'
	regrid_sname = '../../Convolve_Regrid/Conv_Regrids/bkgd100_modled_by_70_ext_20_tstep_'+str(step)+'_to_'

	# Being messy again sorry
	#bkgd_conv = '../../Sky_Remove/Median_Removed/100um_medianRemoved.fits'
	#conv_sname = '../../Convolve_Regrid/Convolutions/SkyRemovd100_to_'
	#regrid_sname = '../../Convolve_Regrid/Conv_Regrids/SkyRemovd100_to_'	
	# After Convolution
	if convolve_regrid == False:
		bkgd_100to160 = fits.open(regrid_sname+'160_CR.fits')[0].data
# Turn everything off to convolve and then you can switch this off forever.
if convolve_regrid:
	k70_100 = '../../Convolve_Regrid/Kernels/Kernel_LoRes_PACS_70_to_PACS_100.fits'
	k70_160 = '../../Convolve_Regrid/Kernels/Kernel_LoRes_PACS_70_to_PACS_160.fits'
	k100_160 = '../../Convolve_Regrid/Kernels/Kernel_LoRes_PACS_100_to_PACS_160.fits'
	im100 =  '../../Original_Files/Infrared/e0102_pacs100_new.fits'
	im160 = '../../Original_Files/Infrared/e0102_pacs160_new.fits'

	if im70_to_100:
		Convolve.master_convolve(k70_100,bkgd_conv,conv_sname+'100.fits')
		regrid.resample(conv_sname+'100.fits',im100,regrid_sname+'100_CR.fits')

	if im70_to_160:
		Convolve.master_convolve(k70_160,bkgd_conv,conv_sname+'160.fits')
		regrid.resample(conv_sname+'160.fits',im160,regrid_sname+'160_CR.fits')

	if im100_to_160:
		print("Convolving 100 to 160")
		Convolve.master_convolve(k100_160,bkgd_conv,conv_sname+'160.fits')
		regrid.resample(conv_sname+'160.fits',im160,regrid_sname+'160_CR.fits')

# Shouldn't have to change anything else in here if you're only changing the
# extension on the 70 mask. To add 24 just copy this stuff and change what's appropriate.
if im70_to_100:
	print("Image is 100 Microns")
	bkgd = bkgd_70to100
	img = fits.open('../../Sky_Remove/Median_Removed/100um_medianRemoved.fits')[0].data
	hdr = fits.open('../../Sky_Remove/Median_Removed/100um_medianRemoved.fits')[0].header
	ra, dec = np.loadtxt('../../Masks/Xray/Mask_World_100.txt')
	w = wcs.WCS(hdr)
	x,y = w.all_world2pix(ra,dec,1)
	vmin=0;vmax=30 #bkgd
	vmin_=-10; vmax_=100 #img
	vmin_sub = 4; vmax_sub= 30 # subtraction
	xdim = [170,230]; ydim = [160,220]
	m1 = 1; m2 = 64;
	fname = '100um_70modeled'
	sname = '100um_70modeled_ForceSlope'
	sq_row = [200,200,200,201,201,201,202,202,202]
	sq_col = [198,199,200,198,199,200,198,199,200]
if im70_to_160:
	print("Image is 160 Microns")
	bkgd = bkgd_70to160
	img = fits.open('../../Sky_Remove/Median_Removed/160um_medianRemoved.fits')[0].data
	hdr = fits.open('../../Sky_Remove/Median_Removed/160um_medianRemoved.fits')[0].header
	ra, dec = np.loadtxt('../../Masks/Xray/Mask_World_160.txt')
	w = wcs.WCS(hdr)
	x,y = w.all_world2pix(ra,dec,1)
	vmin=5;vmax=50 #bkgd
	vmin_=10; vmax_=50 #img
	vmin_sub = 4; vmax_sub= 30 # subtraction
	xdim = [110,150]; ydim = [100,135]
	m1 = 1; m2 = 64;
	fname = '160um_70modeled'
	sname = '160um_70modeled_ForceSlope'
	sq_row = [124,124,124,125,125,125,126,126,126]
	sq_col = [124,125,126,124,125,126,124,125,126]
if im100_to_160:
	print("Image is 160 Microns")
	bkgd = bkgd_100to160
	img = fits.open('../../Sky_Remove/Median_Removed/160um_medianRemoved.fits')[0].data
	hdr = fits.open('../../Sky_Remove/Median_Removed/160um_medianRemoved.fits')[0].header
	ra, dec = np.loadtxt('../../Masks/Xray/Mask_World_160.txt')
	w = wcs.WCS(hdr)
	x,y = w.all_world2pix(ra,dec,1)
	vmin=5;vmax=50 #bkgd
	vmin_=10; vmax_=50 #img
	vmin_sub = 4; vmax_sub= 30 # subtraction
	xdim = [110,150]; ydim = [100,135]
	m1 = 1; m2 = 64;
	fname = '160um_100modeled'	
	sname = '160um_70modeled_ForceSlope'
	sq_row = [124,124,124,125,125,125,126,126,126]
	sq_col = [124,125,126,124,125,126,124,125,126]
###################
# Get Area to Fit #
###################
# Create new array to meddle with

imgShell = np.copy(img)
bkgdShell = np.copy(bkgd)

def prune(img,option):
	# Get flattened array of all the coordinates in the image.
	a,b = np.shape(img)
	row = np.arange(0,a); col = np.arange(0,b)
	row,col=np.meshgrid(row,col)
	row=row.flatten(); col=col.flatten()
	# Turn those coordinates into ra/dec
	all_ra, all_dec = w.all_pix2world(col,row,1)
	# c1 is the origin as defined by karins paper. I got the distance from wikipedia
	# probably not super accurate, but I don't think it needs to be right now.
	c1 = SkyCoord(16.00875*u.deg, -72.03125*u.deg, distance=58*u.kpc, frame='icrs')	
	c2 = SkyCoord(all_ra*u.deg, all_dec*u.deg, distance=58*u.kpc, frame='icrs')
	# Get array for how far each pixel is from the center
	sep = c1.separation_3d(c2)
	# Got the 22 arcsecond from karin's paper to create a radius for the SNR. 
	radius = 58*u.kpc*22/206265
	# This is what sets the outer circle.
	radiusOuter = 58*u.kpc*extend/206265
	# Get rid of all the stuff you dont care about
	look = np.where(sep < radius) # Everything inside snr
	lookAgain = np.where(sep>radiusOuter) #Everything outside outer radius
	waldo = np.where(sep > radius) # Everything outside snr
	if option == "Annulus":
		img[row[look],col[look]] = np.nan
		img[row[lookAgain],col[lookAgain]] = np.nan
	if option == "snr":
		img[row[waldo],col[waldo]] = np.nan
	return img

imgShell = prune(imgShell,"Annulus")		
bkgdShell = prune(bkgdShell,"Annulus")

beforeClip = np.copy(bkgdShell/imgShell)

clipImg = np.copy(bkgdShell/imgShell)  # Create Ratio Img 

med = np.nanmedian(clipImg)
std = np.nanstd(clipImg)

if clip:
	level = 7	# Select Sigma Level 

	ind_r,ind_c = np.shape(clipImg)
	print(("Full Range: {} to {}").format(np.nanmin(clipImg),np.nanmax(clipImg)))

	clipRange = [med - level*std,med  + level*std]

	def clip(data):
		nData = np.copy(data)
		for i in range(ind_r):
			for j in range(ind_c):
				if clipImg[i,j] < clipRange[0]: #or clipImg[i,j] > clipRange[1]:
					nData[i,j] = np.nan
				elif clipImg[i,j] > clipRange[1]:
					nData[i,j] = np.nan
		return nData 

	imgShell = clip(imgShell)
	bkgdShell = clip(bkgdShell)

	med = np.nanmedian(np.copy(bkgdShell/imgShell))
	print(("Fitting Region. Median: {} Std: {}").format(med,std))				
	
	print(("Clipping Range: {} to {}").format(clipRange[0],clipRange[1]))				

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

if slopeInt:
	# Fits data, uses slope and intercept
	slope, intercept, r_value, p_value, std_err = stats.linregress(in_bkgd,in_img)
	line = xvals * slope + intercept
	newbkgd = np.copy(bkgd) * slope + intercept
	print(("For Linear Fit Slope: {} Intercept: {}").format(slope,intercept))
if NoOffset:
	# Uses the median as a slope.
	line = xvals * 1/med
	newbkgd = np.copy(bkgd) * 1/med
	print(("For Linear Fit Median Slope: {} ").format(1/med))
if ForceSlope: # THIS IS NOT DYNAMIC. DONT USE ON 100
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

newbkgd = np.copy(bkgd) * 1/med
newsub = np.copy(img) - newbkgd
newsubSNR = prune(np.copy(newsub),"snr")
plt.figure()

plt.imshow(newsubSNR,vmin = np.nanmin(newsubSNR),vmax=np.nanmax(newsubSNR))
plt.xlim(xdim); plt.ylim(ydim)
plt.scatter(sq_col,sq_row,c='r',s=0.5)
histsub = prune(np.copy(newsub),"snr")
histdata = histsub[np.where(np.isfinite(histsub))].flatten()




##############
# Save Game? #
##############

if save:

	fits.writeto(fname+'/'+sname+'_snr.fits',newsub,hdr,overwrite=True)
	fits.writeto(fname+'/'+sname+'_bkgd.fits',newbkgd,hdr,overwrite=True)

########################
# Plotting and Testing #
########################
# Check that everythings imported and fine. 
if plt_checkMask: 
	f, (ax,bx) = plt.subplots(1,2)
	ax.imshow(bkgd,vmin=vmin,vmax=vmax)
	ax.scatter(x,y,s=.1)
	bx.imshow(img,vmin=vmin_,vmax=vmax_)
	bx.scatter(x,y,s=.1)
	ax.set_xlim(xdim); ax.set_ylim(ydim)
	bx.set_xlim(xdim); bx.set_ylim(ydim)

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
subSum = prune(newsub,"snr")
print(("Sum of Subtraction: {}").format(np.nansum(subSum)))	
subSnr = subSum[np.where(np.isfinite(newsub))]
print(("Sum of Subtraction Negatives: {}").format(np.nansum(subSnr[np.where(subSnr<0)])))	
plt.show()