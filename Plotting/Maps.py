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
# Plot Integrated SED
plot_integrated = False
# Plot Chi Squared Confidence Intervals for Integrated SED 
plot_integrated_ChiSquaredConfidence = False
# Plot Pixel by Pixel SED. Interactive
plot_interactive_pixbypix = False
# Plot 4 Background Removed and Convolved/Regridded to 160 images.
plot_imgs = False
# Plot Temperature, Warm and Cold dust mass maps.
plot_maps = False
# Plot Error from Fit Maps for each image. 
plot_sigma = False
# Plot Templates
plot_template = False

# Plot SEDs For Warm Mass Peak and two blobs. 
plot_SelectSED = False
plot_FewPixels = True


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

#########################
# Load integrated files #
#########################

cold_sed = np.loadtxt("../SED_Fitting/Sols/Integrated/cold_sed.txt")
warm_sed = np.loadtxt("../SED_Fitting/Sols/Integrated/warm_sed.txt")
chi_squared_cube = np.load("../SED_Fitting/Sols/Integrated/chi_squared_cube.npy")
temp,cold_mass,warm_mass,chi_squared = np.loadtxt("../SED_Fitting/Sols/Integrated/temp_coldmass_warmmass_chisqrd.txt")
total_sed = cold_sed + warm_sed

#################
# Load Template #
#################	
template_0n_nan = np.loadtxt("../SED_Fitting/Sols/Templates/Template_0_with_NaNs.txt")
template_2n_nan = np.loadtxt("../SED_Fitting/Sols/Templates/Template_2_with_NaNs.txt")
template_3n_nan = np.loadtxt("../SED_Fitting/Sols/Templates/Template_3_with_NaNs.txt")
template_3n_nanplus = np.loadtxt("../SED_Fitting/Sols/Templates/Template_3_with_NaNs_Plus.txt")
template_5n_nanplus = np.loadtxt("../SED_Fitting/Sols/Templates/Template_5_with_NaNs.txt")

template = template_3n_nanplus
##########
# Plot   #
##########
if plot_integrated:
	fig = plt.figure(figsize=(11,8))
	plot = fig.add_subplot(111)
	# Plot SEDs
	plot.plot(lam,total_sed,color="#424186")
	plot.plot(lam,warm_sed,color="#84D44B",ls='dashed') 
	plot.plot(lam,cold_sed,color="#23A883")
	# Plot Measured Values
	plot.errorbar([24,70,100,160],AverageIntensities,yerr=AverageError,marker='o',linestyle='none',color="black")
	# Plot Labels
	plot.set_xlabel("Wavelength ($\mu m$)",size=18)
	plot.set_ylabel("Spectral Intensity (MJy sr$^{-1}$)",size=18)
	plot.set_title("Integrated Spectral Energy Distribution",size=20)
	plot.legend(("Total SED","Warm SED","Cold SED"),prop={'size':14})
	
	plot.tick_params(axis='both', which='major', labelsize=16)
	plot.tick_params(axis='both', which='minor', labelsize=14)
	
	plot.grid(color='white',linestyle='-')
	plot.set_facecolor("#EAEAF2")

	print(("Temp {} Cold Mass {} Warm Mass {} Total Mass {} Chi Squared {} ").format(temp,cold_mass,warm_mass,cold_mass+warm_mass,chi_squared))
	if save:
		plt.savefig("Plots/IntegratedSED.png")

if plot_integrated_ChiSquaredConfidence:

	def ChiSquaredMap(Map,interval):
		Map = np.copy(Map)
		# Change everything outside the interval to NaN.
		R,C = np.where(Map > (chi_squared + interval))
		Map[R,C] = np.nan 
		return np.copy(Map)

	f, axes = plt.subplots(3,2)
	
	# What the elements actually represent on an axis
	physical_X = [np.where(np.isclose(ColdMass,cold_mass))[0][0],
					np.where(np.isclose(ColdTemp,temp))[0][0],
						np.where(np.isclose(ColdMass,cold_mass))[0][0]]

	physical_Y = [np.where(np.isclose(WarmMass,warm_mass))[0][0],
					np.where(np.isclose(WarmMass,warm_mass))[0][0],
						np.where(np.isclose(ColdTemp,temp))[0][0]]
	
	titles = ["68.3% Confidence","95.4% Confidence"]
	ylabels = ["Warm Dust Mass M$_\odot$","Warm Dust Mass M$_\odot$","Temperature K",]
	xlabels = ["Cold Dust Mass M$_\odot$", "Temperature K", "Cold Dust Mass M$_\odot$"]
	
	# Temperature Slice, Cold Mass Slice, Warm Mass Slice at solutions
	chi_squared_slices = [np.copy(chi_squared_cube)[:,physical_X[1],:],
							np.copy(chi_squared_cube)[:,:,physical_X[0]],
								np.copy(chi_squared_cube)[physical_Y[0],:,:]]
	
	for i in range(3):
		for j in range(2):
			img = ChiSquaredMap(chi_squared_slices[i],intervals[j])
			axes[i,j].imshow(img)
			axes[i,j].scatter(physical_X[i],physical_Y[i],s=5,c='r')
			axes[i,j].set_ylim(physical_Y[i]-10,physical_Y[i]+10); axes[i,j].set_xlim(physical_X[i]-10,physical_X[i]+10);
			axes[i,j].set_xlabel(xlabels[i]); axes[i,j].set_ylabel(ylabels[i])
   
	axes[0,0].set_title(titles[0])
	axes[0,1].set_title(titles[1])  
if plot_interactive_pixbypix:

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
if plot_imgs:    
	f, axes = plt.subplots(2,2,figsize=(10,15))
	x1 = 116; x2 = 140; y1 = 112; y2 = 128
	vmin = 0; vmax = 10
	axes[0,0].imshow(fits.open(files[0])[0].data,vmin=vmin,vmax=vmax)
	axes[0,1].imshow(fits.open(files[1])[0].data,vmin=vmin,vmax=vmax)
	axes[1,0].imshow(fits.open(files[2])[0].data,vmin=vmin,vmax=vmax)
	axes[1,1].imshow(fits.open(files[3])[0].data,vmin=vmin,vmax=vmax)
	axes[0,0].set_xlim(x1,x2); axes[0,0].set_ylim(y1,y2); axes[0,0].axis("off")
	axes[0,1].set_xlim(x1,x2); axes[0,1].set_ylim(y1,y2); axes[0,1].axis("off")
	axes[1,0].set_xlim(x1,x2); axes[1,0].set_ylim(y1,y2); axes[1,0].axis("off")
	axes[1,1].set_xlim(x1,x2); axes[1,1].set_ylim(y1,y2); axes[1,1].axis("off") 
	 
	axes[0,0].set_title("24 $\mu$m",fontsize="x-large"); axes[0,1].set_title("70 $\mu$m",fontsize="x-large")
	axes[1,0].set_title("100 $\mu$m",fontsize="x-large"); axes[1,1].set_title("160 $\mu$m",fontsize="x-large")
	f.set_size_inches(5, 5)
	plt.subplots_adjust(wspace=-.15,hspace=-.1,top=1,bottom=0,right=1,left=0)

	if save:
		plt.savefig("Plots/E0102_Regrid_Convolved.png",dpi=200, bbox_inches='tight', pad_inches = 0 )



if plot_maps: # Temperature Map
	# Prep Temp Map for Mean Measurement by only including relevant pixels
	a,b = np.shape(Pix_Temp)
	temp = np.copy(Pix_Temp)
	for i in range(a):
		for j in range(b):
			# If a value is the temperature map is exactly 0
			# (i.e. falls out of the radius) then make it a nan
			if Pix_Temp[i,j] == 0.:
				temp[i,j] = np.nan
	temp = temp * template			
	g, (ax,bx,cx) = plt.subplots(1,3)
	tim = ax.imshow(Pix_Temp * template ,vmin = np.nanmin(temp), vmax = np.nanmax(temp))
	ax.set_title(("Temperature Map"))
	ax.set_ylim(110,130)
	ax.set_xlim(115,140)
	ax.axis("off")
	g.colorbar(tim, ax=ax)

	wim = bx.imshow(Pix_Warm_Mass * template,vmin = np.nanmin(Pix_Warm_Mass*template), vmax = np.nanmax(Pix_Warm_Mass*template))
	bx.set_title(("Warm Mass Map"))
	bx.set_ylim(110,130)
	bx.set_xlim(115,140)
	bx.axis("off")
	g.colorbar(wim, ax=bx)

#	Pix_Cold_Mass[126,133] = np.nan; Pix_Cold_Mass[124,135] = np.nan;
#	Pix_Cold_Mass[123,135] = np.nan; Pix_Cold_Mass[123,134] = np.nan; Pix_Cold_Mass[122,133] = np.nan;
#	cim = cx.imshow(Pix_Cold_Mass *template ,vmin = np.min(Pix_Cold_Mass), vmax = np.max(Pix_Cold_Mass))
	cim = cx.imshow(Pix_Cold_Mass  * template,vmin = .00025, vmax = .002)
	cx.set_title(("Cold Mass Map"))
	cx.set_ylim(110,130)
	cx.set_xlim(115,140)
	cx.axis("off")
	g.colorbar(cim, ax=cx)
	print(np.nanmin(Pix_Cold_Mass  * template))
	print(np.nanmax(Pix_Cold_Mass  * template))




if plot_SelectSED:
	# Initialize Plot
	f = plt.figure()

	cx = plt.subplot2grid((3,3),(0,2))#,colspan=2)

	dx = plt.subplot2grid((3,3),(1,0))#,colspan=2)
	ex = plt.subplot2grid((3,3),(1,1))

	ix = plt.subplot2grid((3,3),(2,2))#,colspan=2)

	plots = [cx,dx,ix]

	# Data 

	img = data[3]
	
	Px = [128,133,123]; Py = [116,121,118]
	# Plot 
	cmap = plt.get_cmap('viridis')
	
	# SNR Image
	ex.imshow(img,vmin=np.nanmin(img),vmax=np.nanmax(img))
	ex.set_xlim(xdim)
	ex.set_ylim(ydim)
	ex.axis("off")
	ex.scatter(Px,Py,s=40, marker='s',facecolors='none',edgecolors='r')

	for i in range(3):
		ObsvSed = [data[0][Py[i],Px[i]],data[1][Py[i],Px[i]],data[2][Py[i],Px[i]],data[3][Py[i],Px[i]]]
		ObsvErr = Eqs.error(ObsvSed)
		plots[i].plot(lam,Pix_Warm_SED[Py[i],Px[i]]+Pix_Cold_SED[Py[i],Px[i]],label="Total SED",color="#424186")
		plots[i].plot(lam,Pix_Warm_SED[Py[i],Px[i]],label="Warm SED",color="#84D44B",ls='dashed')
		plots[i].plot(lam,Pix_Cold_SED[Py[i],Px[i]],label="Cold SED",color="#23A883")
		
		plots[i].set_xlabel("Wavelength ($\mu m$)")
		plots[i].set_ylabel("Spectral Intensity (MJy sr$^{-1}$)")
	
		plots[i].errorbar(wv,ObsvSed,yerr=ObsvErr,marker='o',linestyle='none',c='black')
		plots[i].grid(color='white',linestyle='-')
		plots[i].set_facecolor("#EAEAF2")
		plt.legend()

if plot_FewPixels:
	# Initialize Plot
	f = plt.figure()

	cx = plt.subplot2grid((3,2),(0,1))#,colspan=2)
	dx = plt.subplot2grid((3,2),(1,1))#,colspan=2)
	ix = plt.subplot2grid((3,2),(2,1))#,colspan=2)

	# Image
	ex = plt.subplot2grid((3,2),(0,0),rowspan=3)

	plots = [cx,dx,ix]
	colors = ["red","aqua","yellow"]
	# Data 

	img = data[3]
	
	#Px = [128,133,123]; Py = [116,121,118]
	Px = [133,123,128]; Py = [121,118,116]
	# Plot 
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
	plt.subplots_adjust(top=0.88, bottom=0.11, left=0.11, right=0.9, hspace=0.2, wspace=0.2)

if plot_sigma:
	a,b = np.shape(data[0])
	# Set up maps
	sigma24 = np.zeros(a*b).reshape(a,b); sigma70 = np.copy(sigma24); sigma100 = np.copy(sigma24); sigma160 = np.copy(sigma24); 
	# Loop over real values, find error, plug in error value into map.
	for i in range(a):
		for j in range(b):

			if np.isfinite(data[0][i,j]) and np.isfinite(data[1][i,j]) and np.isfinite(data[2][i,j]) and np.isfinite(data[3][i,j]):

				arr = [data[0][i,j],data[1][i,j],data[2][i,j],data[3][i,j]]
				err = Eqs.error(arr)

				sigma24[i,j] = err[0]; sigma70[i,j] = err[1];
				sigma100[i,j] = err[2]; sigma160[i,j] = err[3];
	# Plot
	f, (ax,bx,cx,dx) = plt.subplots(1,4)				
	vmin = 0; vmax = 4
	# Show Plots
	ax.imshow(sigma24,vmin=vmin,vmax=vmax)
	bx.imshow(sigma70,vmin=vmin,vmax=vmax)
	cx.imshow(sigma100,vmin=vmin,vmax=vmax)
	cb = dx.imshow(sigma160,vmin=vmin,vmax=vmax)
	# Set Limits
	ax.set_xlim(xdim[0],xdim[1])
	ax.set_ylim(ydim[0],ydim[1])
	bx.set_xlim(xdim[0],xdim[1])
	bx.set_ylim(ydim[0],ydim[1])
	cx.set_xlim(xdim[0],xdim[1])
	cx.set_ylim(ydim[0],ydim[1])
	dx.set_xlim(xdim[0],xdim[1])
	dx.set_ylim(ydim[0],ydim[1])
	# Set colorbar
	cax = f.add_axes([0.27, 0.8, 0.5, 0.05])
	f.colorbar(cb,cax=cax,orientation='horizontal')

	plt.figure(2)
	plt.imshow(sigma24+sigma70+sigma100+sigma160,vmin=vmin,vmax=8)
	plt.xlim(xdim)
	plt.ylim(ydim)
	plt.colorbar()
if plot_template: 

	box_r, box_c = np.where(template == 0)	
	
	boxes = []
	for x,y in zip(box_c, box_r): 

		rect = Rectangle((x-.5,y-.5),1,1)
		boxes.append(rect)

	#pc = PatchCollection(boxes,facecolor="none",edgecolor="red")
	
	# What does the histogram of errors look like without the threshhold? 
	sigma_fit = 1;# num_bins = 25
	# Load Stuff
	Temperature_Confidence = np.load("Sols/PixbyPix/Temperature_Confidence.npy")[:,:,sigma_fit]
	ColdMass_Confidence = np.load("Sols/PixbyPix/Cold_Mass_Confidence.npy")[:,:,sigma_fit] 
	WarmMass_Confidence = np.load("Sols/PixbyPix/Warm_Mass_Confidence.npy")[:,:,sigma_fit] 

	Temp_err = Temperature_Confidence[(Temperature_Confidence != 0) & (np.isfinite(Temperature_Confidence))]
	Cold_err = ColdMass_Confidence[(ColdMass_Confidence != 0) & (np.isfinite(ColdMass_Confidence))]
	Warm_err = WarmMass_Confidence[(WarmMass_Confidence != 0) & (np.isfinite(WarmMass_Confidence))]

	Temp_err_clip = (Temperature_Confidence * template_3n_nan)[(Temperature_Confidence != 0) & (np.isfinite(Temperature_Confidence))]
	Cold_err_clip = (ColdMass_Confidence * template_3n_nan)[(ColdMass_Confidence != 0) & (np.isfinite(ColdMass_Confidence))]
	Warm_err_clip = (WarmMass_Confidence * template_3n_nan)[(WarmMass_Confidence != 0) & (np.isfinite(WarmMass_Confidence))]
	
	#########
	# Plots #
	#########

	p = [data[3],Temperature_Confidence,ColdMass_Confidence,WarmMass_Confidence]
	# Which pixels are being removed by the threshhold?
	f, axs = plt.subplots(1,4)
	t = ["24 Microns","Temperature Error Intervals","Cold Mass Error Intervals","Warm Mass Error Intervals"]
	v = [0.5,1,8e-5,8e-9]; v_= [9,17,.01,8e-8]
	for i in range(4):
		pc = PatchCollection(boxes,facecolor="none",edgecolor="red")
		p[i][p[i] == 0] = np.nan
		axs[i].imshow(p[i],vmin=v[i],vmax=v_[i]); axs[i].set_xlim(118,138); axs[i].set_ylim(110,130)
		axs[i].add_collection(pc); axs[i].axis('off'); axs[i].set_title(t[i])

	f.suptitle("2 Sigma Error Distributions, Clipped by 3 * $\sigma_{sky}$")
	g, axes = plt.subplots(3,2)

	plots = [Temp_err,Temp_err_clip,Cold_err,Cold_err_clip,Warm_err,Warm_err_clip ]
	axes[1,0].set_ylabel("Cold Mass"); axes[2,0].set_ylabel("Warm Mass"); axes[0,0].set_ylabel("Temperature")

	g.suptitle("Chi Squared Confidence Intervals", fontsize=16)
	axes[0,0].set_title("Before Clipping")
	axes[0,1].set_title("After Clipping")


	def histplot(array,row,col,num_bins):
		# if row == 0 and col == 1:
		# 	num_bins = 9
		# if row == 1 and col == 0:
		# 	num_bins = 15
		# elif row == 2 and col == 0: 
		# 	num_bins = 6
		# else:
		# 	# Freedman Diaconis rule for number of bins
		# 	num_bins = int((np.max(array) - np.min(array)) / (2 * iqr(array) * len(array)**(-1/3)))
		print(num_bins)
		# Plot histogram
		n, bins, patches = axes[row,col].hist(array,bins=num_bins,density=True)
		# Fit histogram
		(mu, sigma) = norm.fit(array)
		hist_y = mlab.normpdf(bins,mu,sigma)
		# Plot Fit
		axes[row,col].plot(bins, hist_y, '--')

	count = 0 ; nbins = [11,7,15,15,6,8]
	for i in range(3):
		for j in range(2):
			histplot(plots[count],i,j,nbins[count])
			count += 1
plt.show()
