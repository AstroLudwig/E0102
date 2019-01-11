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

save = False
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
# Plot Temperature Map for Detection Limit. 
# Can also turn on use_limit as look at the other plots.
plot_Selecttmap = False
# Plot Pixel by Pixel Chi Squared Contours. Interactive
plot_chi = False
# Plot both Pixel by Pixel SED and Chi Contours. Interactive.
plot_all = False
# Plot Sigma Maps for each image. 
plot_sigma = False
# Plot Contour Maps for Selected Pixels. 
plot_SelectContour = False
# Generate signal to noise maps without calibration or background error
generate_SN = False
# Plot SEDs For Warm Mass Peak and two blobs. 
plot_SelectSED = False
# Plot Templates
plot_template = True

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

lam, k_amc = Eqs.log_extrap('Kappa/kappa_amc.dat')
lam, k_mg2 = Eqs.log_extrap('Kappa/kappa_mg2sio4.dat')
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

Pix_Warm_SED =  np.load("Sols/PixbyPix/Warm_SED.npy")
Pix_Cold_SED =  np.load("Sols/PixbyPix/Cold_SED.npy")

# Temperature                                         
Pix_Temp = np.loadtxt("Sols/PixbyPix/Temperature.txt")

# Mass 
Pix_Warm_Mass = np.loadtxt("Sols/PixbyPix/WarmMass.txt")
Pix_Cold_Mass = np.loadtxt("Sols/PixbyPix/ColdMass.txt")   

# Chi Squared Maps                            
Pix_chisqrd = np.loadtxt("Sols/PixbyPix/ChiSquared.txt")

#########################
# Load integrated files #
#########################

cold_sed = np.loadtxt("Sols/Integrated/cold_sed.txt")
warm_sed = np.loadtxt("Sols/Integrated/warm_sed.txt")
chi_squared_cube = np.load("Sols/Integrated/chi_squared_cube.npy")
temp,cold_mass,warm_mass,chi_squared = np.loadtxt("Sols/Integrated/temp_coldmass_warmmass_chisqrd.txt")
total_sed = cold_sed + warm_sed

#################
# Load Template #
#################	

template_3n_nan = np.loadtxt("Sols/Templates/Template_3_with_NaNs.txt")

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
	plot.set_ylabel("Spectral Intensity (Mjy sr$^{-1}$)",size=18)
	plot.set_title("Average Spectral Energy Distribution",size=20)
	plot.legend(("Total SED","Warm SED","Cold SED"),prop={'size':14})
	
	plot.tick_params(axis='both', which='major', labelsize=16)
	plot.tick_params(axis='both', which='minor', labelsize=14)
	
	plot.grid(color='white',linestyle='-')
	plot.set_facecolor("#EAEAF2")

	print(("Temp {} Cold Mass {} Warm Mass {} Total Mass {} Chi Squared {} ").format(temp,cold_mass,warm_mass,cold_mass+warm_mass,chi_squared))
	if save:
		plt.savefig("AverageSED.png")

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
	bx.plot(lam,Pix_Total_SED[0,0],label="Total SED")
	bx.plot(lam,Pix_Warm_SED[0,0],label="Warm SED")
	bx.plot(lam,Pix_Cold_SED[0,0],label="Cold SED")
	bx.errorbar(wv,inObsvSed,yerr=inObsvErr,marker='o',linestyle='none',c='purple')
	ax.set_xlim(xdim[0],xdim[1])
	ax.set_ylim(ydim[0],ydim[1])
	bx.set_ylim(0, np.max(Pix_Total_SED[0,0])+2)
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

		bx.plot(lam,Pix_Total_SED[i,j],label="Total SED")
		bx.plot(lam,Pix_Warm_SED[i,j],label="Warm SED")
		bx.plot(lam,Pix_Cold_SED[i,j],label="Cold SED")
		
		bx.set_ylim(0, np.max(Pix_Total_SED[i,j])+2)
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
				
	g, (ax,bx,cx) = plt.subplots(1,3)
	tim = ax.imshow(Pix_Temp ,vmin = np.nanmin(temp), vmax = np.nanmax(temp))
	ax.set_title(("Temperature Map.\n Average is {}").format(int(np.nanmean(temp))))
	ax.set_ylim(110,130)
	ax.set_xlim(115,140)
	ax.axis("off")
	g.colorbar(tim, ax=ax)

	wim = bx.imshow(Pix_Warm_Mass ,vmin = np.min(Pix_Warm_Mass), vmax = np.max(Pix_Warm_Mass))
	bx.set_title(("Warm Mass Map.\n Sum is {:2f}").format(np.sum(Pix_Warm_Mass)))
	bx.set_ylim(110,130)
	bx.set_xlim(115,140)
	bx.axis("off")
	g.colorbar(wim, ax=bx)

	cim = cx.imshow(Pix_Cold_Mass ,vmin = np.min(Pix_Cold_Mass), vmax = np.max(Pix_Cold_Mass))
	cx.set_title(("Cold Mass Map.\n Sum is {:2f}").format(np.sum(Pix_Cold_Mass)))
	cx.set_ylim(110,130)
	cx.set_xlim(115,140)
	cx.axis("off")
	g.colorbar(cim, ax=cx)


if plot_Selecttmap: # Selected Temperature Map
	SelTmap = np.loadtxt("Sols/2SigmaTemperatureMap.txt")
	plt.figure(4)
	plt.imshow(SelTmap,vmin = 25, vmax = 50)
	plt.title("Temperature Map")
	plt.ylim(110,130)
	plt.xlim(115,140)

if plot_chi:

	T = np.arange(2,70,1) # Kelvin
	# Fractional Mass of the Sun using Dex
	# Cold Mass
	M = 10**np.arange(-4,0.1,.1)
	T,M = np.meshgrid(T,M)
	chiMap = np.load("Sols/Parameter_Chi_Map.npy")

	

	f, (ax,bx) = plt.subplots(1,2)
	
	ax.imshow(data[3])
	ax.scatter(0,0,c='r',marker='s',s=40)
	ax.axis("off")
	ax.set_xlim(xdim[0],xdim[1])
	ax.set_ylim(ydim[0],ydim[1])
	ax.set_title("Chi Contours for a Fixed Warm Value Mass of 1.97E-7")
	
	cmap = plt.cm.get_cmap("winter")

	cs = bx.contourf(M,T,chiMap[0,0].transpose(),cmap=cmap)
	bx.set_xlabel("Mass Range in Solar Fractional Mass")
	bx.set_ylabel("Temperature in Kelvin")
	bx.grid()


	bx.set_title(("Temperature: {:2d} K,  Cold Mass: {:.2E} $M_\odot$,  Warm Mass: {:.2E} $M_\odot$,  $\chi^2$: {:.2f}").format(int(Pix_Temp[0,0]),Pix_Cold_Mass[0,0],Pix_Warm_Mass[0,0],Pix_chisqrd[0,0]))
	


	def onclick(event):
		i = int(event.ydata)
		j = int(event.xdata)

		ObsvSed = [data[0][i,j],data[1][i,j],data[2][i,j],data[3][i,j]]
		ObsvErr = Eqs.error(ObsvSed)        
		ax.clear()
		bx.clear()
		#f.delaxes(f.axes[2]) 

		ax.imshow(data[3])
		ax.set_title("Chi Contours for a Fixed Warm Value Mass of 1.97E-7")
		ax.scatter(j,i,c='r',marker='s',s=40)
		ax.axis("off")
		ax.set_xlim(xdim[0],xdim[1])
		ax.set_ylim(ydim[0],ydim[1])
	
		minChi_y, minChi_x = np.where(chiMap[i,j]==np.min(chiMap[i,j]))	
		cs = bx.contourf(M,T,np.log(chiMap[i,j].transpose()),cmap=cmap)#,levels=(0,1,2))
		print("Min Chi in Contour Map: " + str(chiMap[minChi_y,minChi_x]))
		#cb = f.colorbar(cs)
		#bx.imshow(chiMap[i,j],vmin=0,vmax=100)
		#bx.scatter(minChi_x,minChi_y,s=20,c='r')
		#bx.set_xlim(M_[0],M_[len(M_)-1])
		#bx.grid()
		#bx.set_ylim(2,70)
		

		bx.set_xlabel("Mass Range in Solar Fractional Mass")
		bx.set_ylabel("Temperature in Kelvin")
		bx.set_title(("Temperature: {:2d} K,  Cold Mass: {:.2E} $M_\odot$,  $\chi^2$: {:.2f}").format(int(Pix_Temp[i,j]),Pix_Cold_Mass[i,j],Pix_chisqrd[i,j]))
		#cx.imshow(np.log(chiMap[i,j]))
		#cb.remove()

		plt.draw()

		
	cid = f.canvas.mpl_connect('button_press_event', onclick)	

if plot_all:
	if use_all or use_noBkgdErr:
		img = data[3]
	if use_limit:
		img = fits.open("Sols/DetectionLimited_160umImg.fits")[0].data
	
	
	# Parameters
	T_ = np.arange(2,70,1) # Kelvin
	# Fractional Mass of the Sun using Dex
	# Cold Mass
	M_ = 10**np.arange(-4,0.1,.1)
	T,M = np.meshgrid(T_,M_)
	
	
	
	# Contour Levels
	lvls = np.arange(5.9,6.5,.01)

	# Set up plot. 
	cmap = plt.get_cmap('viridis')
	f = plt.figure()
	ax = plt.subplot2grid((2,2),(0,0),rowspan=2)
	bx = plt.subplot2grid((2,2),(0,1),colspan=1)
	cx = plt.subplot2grid((2,2),(1,1),colspan=3)
	
	# Get 4 Intensity Values
	chiMap = np.load("Sols/Parameter_Chi_Map.npy")
	inObsvSed = [data[0][0,0],data[1][0,0],data[2][0,0],data[3][0,0]]
	inObsvErr = Eqs.error(inObsvSed)
	
	# SNR Remnant
	ax.imshow(img)
	ax.scatter(0,0,c='r',marker='s',s=40)
	ax.axis("off")
	ax.set_xlim(xdim[0],xdim[1])
	ax.set_ylim(ydim[0],ydim[1])
	ax.set_title(("Temperature: {:2d} K,  Cold Mass: {:.2f} $M_\odot$, \n Warm Mass: {:.2E} $M_\odot$,  $\chi^2$: {:.2f}").format(int(Pix_Temp[0,0]),Pix_Cold_Mass[0,0],Pix_Warm_Mass[0,0],Pix_chisqrd[0,0]))

	# SED Plot 
	bx.plot(lam,Pix_Total_SED[0,0],label="Total SED")
	bx.plot(lam,Pix_Warm_SED[0,0],label="Warm SED")
	bx.plot(lam,Pix_Cold_SED[0,0],label="Cold SED")
	bx.errorbar(wv,inObsvSed,yerr=inObsvErr,marker='o',linestyle='none',c='purple')
	bx.set_ylim(0, np.max(Pix_Total_SED[0,0])+2)
	bx.grid()
	bx.set_xlabel("Wavelength ($\mu$m)")
	bx.set_ylabel("Spectral Intensity (MJy sr$^{-1}$)")
	
	chiPlot = np.log(chiMap[0,0].transpose())
	#chiPlot = chiMap[0,0].transpose()
	cs = cx.contourf(M,T,chiPlot,cmap=cmap,levels=lvls)
	cb = f.colorbar(cs)
	cx.set_xlabel("Mass Range in Solar Fractional Mass")
	cx.set_ylabel("Temperature in Kelvin")
	xr = [M_[0],Pix_Cold_Mass[0,0]+2*Pix_Cold_Mass[0,0]]
	yr = [Pix_Temp[0,0]-.5*Pix_Temp[0,0],Pix_Temp[0,0]+.5*Pix_Temp[0,0]]
	cx.set_xlim(xr)
	cx.set_ylim(yr)
	#cx.grid()


	

	def onclick(event):
		i = int(event.ydata)
		j = int(event.xdata)
		print(("Coordinates {},{}").format(j,i))
		chiPlot = np.log(chiMap[i,j].transpose())
		ObsvSed = [data[0][i,j],data[1][i,j],data[2][i,j],data[3][i,j]]
		ObsvErr = Eqs.error(ObsvSed)        
		ax.clear()
		bx.clear()
		cx.clear()

		ax.imshow(img)
		ax.scatter(j,i,c='r',marker='s',s=40)
		ax.set_xlim(xdim[0],xdim[1])
		ax.set_ylim(ydim[0],ydim[1])
		ax.axis("off")

		bx.plot(lam,Pix_Total_SED[i,j],label="Total SED")
		bx.plot(lam,Pix_Warm_SED[i,j],label="Warm SED")
		bx.plot(lam,Pix_Cold_SED[i,j],label="Cold SED")
		
		bx.set_ylim(0, np.max(Pix_Total_SED[i,j])+2)
		bx.errorbar(wv,ObsvSed,yerr=ObsvErr,marker='o',linestyle='none',c='purple')
		bx.grid()
		#plt.legend()
		ax.set_title(("Temperature: {:2d} K,  Cold Mass: {:.2E} $M_\odot$, \n Warm Mass: {:.2E} $M_\odot$,  $\chi^2$: {:.2f}").format(int(Pix_Temp[i,j]),Pix_Cold_Mass[i,j],Pix_Warm_Mass[i,j],Pix_chisqrd[i,j]))
		bx.set_xlabel("Wavelength ($\mu$m)")
		bx.set_ylabel("Spectral Intensity (MJy sr$^{-1}$)")

		
		#chiPlot = chiMap[i,j].transpose()
		cs = cx.contourf(M,T,chiPlot,cmap=cmap,levels=lvls)
		#cb = f.colorbar(cs)
		cx.set_xlabel("Mass Range in Solar Fractional Mass")
		cx.set_ylabel("Temperature in Kelvin")
		#cx.grid()
		xr = [M_[0],2*Pix_Cold_Mass[i,j]]
		yr = [Pix_Temp[i,j]-.5*Pix_Temp[i,j],Pix_Temp[i,j]+.5*Pix_Temp[i,j]]
		cx.set_xlim(xr)
		cx.set_ylim(yr)
		cx.ticklabel_format(style='sci',scilimits=(-3,4),axis='x')
		plt.draw()

	cid = f.canvas.mpl_connect('button_press_event', onclick)
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

if plot_SelectContour:
	# Initialize Plot
	f = plt.figure()
	#ax = plt.subplot2grid((3,3),(0,0))
	#bx = plt.subplot2grid((3,3),(0,1))
	cx = plt.subplot2grid((3,3),(0,2))

	dx = plt.subplot2grid((3,3),(1,0))
	ex = plt.subplot2grid((3,3),(1,1))
	#fx = plt.subplot2grid((3,3),(1,2))


	#gx = plt.subplot2grid((3,3),(2,0))
	#hx = plt.subplot2grid((3,3),(2,1))
	ix = plt.subplot2grid((3,3),(2,2))

	#plots = [ax,bx,cx,dx,fx,gx,hx,ix]
	plots = [cx,dx,ex,ix]
	# Data 

	img = data[3]
	chiMap = np.load("Sols/Parameter_Chi_Map.npy")
	stats = np.loadtxt("Sols/SelectPixelStats.txt")
	t_pix = stats[:,0]; cm_pix = stats[:,1]; wm_pix= stats[:,2]; chi_pix= stats[:,3]
	def cp (x,y):
		return chiMap[y,x].transpose() #np.log(chiMap[y,x].transpose())

	# Parameters
	
	T_ = np.arange(2,70,1) # Kelvin
	# Fractional Mass of the Sun using Dex
	# Cold Mass
	M_ = 10**np.arange(-4,0.1,.1)
	T,M = np.meshgrid(T_,M_)
	
	# Pixels I'm selecting arbitrarily. 
	#Px = [124,128,133,121,136,123,128,133]; Py = [126,127,126,120,121,115,113,115]
	Px = [128,133,123]; Py = [116,121,118]
	# Plot 
	cmap = plt.get_cmap('viridis')
	
	# SNR Image
	ex.imshow(img,vmin=np.nanmin(img),vmax=np.nanmax(img))
	ex.set_xlim(xdim)
	ex.set_ylim(ydim)
	ex.axis("off")

	#xDim_ = [1e-4,4e-4,1.8e-3,4e-4,  1e-3,2e-4,  1.5e-4,6e-4]
	#xDim =  [9e-4,3e-3,9e-3,  1.4e-3,9e-3,1.5e-3,8e-4,2.2e-3]
	#yDim_ = [2,21.4,19.7,30.5,19.5,32,32.8,27]
	#yDim = [20,31,26,38.2,25.5,41.2,42,33.5]
	
	#lvls = np.arange(0,9,1)
	# Contour Maps
	factor = 1
	for i in range(3):
		# Paramter range
		#plots[i].set_xlim(xDim_[i],xDim[i])
		#plots[i].set_ylim(yDim_[i],yDim[i])
		# Colorbar 
		ConMap = cp(Px[i],Py[i])
		lvls = np.arange(np.min(ConMap),np.min(ConMap)+factor,.1)
		cplot = plots[i].contourf(M,T,ConMap,cmap=cmap ,levels=lvls)
		#plots[i].set_title(np.min(chiMap[Py[i],Px[i]]))
		if plots[i] != ex:
			f.colorbar(cplot,ax=plots[i])
		plots[i].ticklabel_format(style='sci',scilimits=(-3,4),axis='x')
		plots[i].scatter(cm_pix[i],t_pix[i],c='r',marker='*',s=5)

	WarmMassMap = np.loadtxt("Sols/WarmMassMap.txt")
	PixIntArr = np.zeros(4*8).reshape(8,4); stats = np.copy(PixIntArr)
	for i in range(3):
		dataArr = np.zeros(4)
		for j in range(4):
			dataArr[j] = data[j][Py[i],Px[i]]
		PixIntArr[i] = dataArr
		#total_Sed, warm_Sed, cold_Sed, best_cold_temp, best_cold_mass, best_warm_mass, best_error
		#ts,ws,cs,bct,bcm,bwm,bChi = Eqs.CalculateBestSed(T_,145,M_,WarmMassMap[Py[i],Px[i]],k_amc,k_mg2,lam,pix_area,dataArr)
		#print(("For pixel in x: {},y:{} ... Chi value is: {}").format(Px[i],Py[i],bChi))
		#stats[i] = [bct,bcm,bwm,bChi]
	#np.savetxt("Sols/SelectPixelStats.txt",stats)
	
	# Turn selected pixel outline red.
	#p = Rectangle((124, 126), 5, 5, fill=False)#, #transform=ax.transAxes, clip_on=False)
	
	#ax.add_patch(p)
	ex.scatter(Px,Py,s=40, marker='s',facecolors='none',edgecolors='r')
if generate_SN:
	
	Noise = np.loadtxt('../Sky_Remove/Sigma.txt')

	# Create binary maps based on these contours. 
	a,b = np.shape(data[0])
	BiMap = []
	for k in range(4):
		shell = np.zeros(a*b).reshape(a,b)
		for i in range(a):
			for j in range(b):
				if data[k][i,j] > 2 * Noise[k]:
					shell[i,j] = 1
				else:
					shell[i,j] = 0
		BiMap.append(shell)

	# Get SED maps based on these binary ones. 
	#AMC_Total_Mass = np.loadtxt("Sols/BkgdErrRemovd_TotalMassMap.txt"); 
	#Pix_Warm_Mass = np.loadtxt("Sols/BkgdErrRemovd_WarmMassMap.txt");
	#Pix_Cold_Mass = np.loadtxt("Sols/BkgdErrRemovd_ColdMassMap.txt");   

	TotalMass = np.copy(AMC_Total_Mass)
	WarmMass = np.copy(Pix_Warm_Mass)
	ColdMass = np.copy(Pix_Cold_Mass)
	TemP = np.copy(Pix_Temp)
	for i in range(a):
		for j in range(b):
			if BiMap[0][i,j] == 0 or BiMap[1][i,j] ==0 or BiMap[2][i,j] == 0 or BiMap[3][i,j] == 0:
				TotalMass[i,j] = np.nan
				WarmMass[i,j] = np.nan
				ColdMass[i,j] = np.nan
				TemP[i,j] = np.nan
	
	
	h, (zx,yx,xx) = plt.subplots(1,3,figsize=(11,4))
	
	LeMap = [TotalMass,WarmMass,ColdMass]
	Titles = ["Temperature Map","Warm Mass Map","Cold Mass Map"]
	
	zPlts = [zx,yx,xx]
	#tcks = [[1.59e-4,.5e-3,1e-3,1.5e-3,2e-3,2.5e-3],[7.94e-8,1e-7,1.5e-7,2e-7,2.5e-7,3e-7,3.5e-7,3.98e-7],[1e-4,5e-4,1e-3,1.5e-3,2e-3,2.5e-3]]
	
	for i in range(3):
		Jim = zPlts[i].imshow(LeMap[i])
		diff = (np.nanmax(LeMap[i])-np.nanmin(LeMap[i]))/5
		tcks = np.arange(np.nanmin(LeMap[i]),np.nanmax(LeMap[i])+diff,diff)
		zPlts[i].set_xlim(xdim)
		zPlts[i].set_ylim(ydim)
		zPlts[i].set_title(Titles[i])
		zPlts[i].axis("off")
		h.colorbar(Jim, ax=zPlts[i],fraction=.04,ticks=tcks,format='%.1e')#,format='%.0e',ticks=range(6))
	
	plt.subplots_adjust(left=0.0, right=0.9, top=0.983,bottom=0.053,hspace=0.155,wspace=0.089)	

	plt.savefig("Maps.png")	
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
		plots[i].plot(lam,Pix_Total_SED[Py[i],Px[i]],label="Total SED",color="#424186")
		plots[i].plot(lam,Pix_Warm_SED[Py[i],Px[i]],label="Warm SED",color="#84D44B",ls='dashed')
		plots[i].plot(lam,Pix_Cold_SED[Py[i],Px[i]],label="Cold SED",color="#23A883")
		
		plots[i].set_xlabel("Wavelength ($\mu m$)")
		plots[i].set_ylabel("Spectral Intensity (MJy sr$^{-1}$)")
	
		plots[i].errorbar(wv,ObsvSed,yerr=ObsvErr,marker='o',linestyle='none',c='black')
		plots[i].grid(color='white',linestyle='-')
		plots[i].set_facecolor("#EAEAF2")
		plt.legend()
if plot_template: 

	# Which pixels are being removed by the threshhold?
	f, ax = plt.subplots(1)

	box_r, box_c = np.where(template_3n_nan == 0)	
	
	boxes = []
	for x,y in zip(box_c, box_r): 

		rect = Rectangle((x-.5,y-.5),1,1)
		boxes.append(rect)

	pc = PatchCollection(boxes,facecolor="none",edgecolor="red")


	stacked_data = np.nansum(np.copy(data),axis=0)
	ax.imshow(data[3],vmin=-1,vmax=8); ax.set_xlim(118,138); ax.set_ylim(110,130)
	ax.add_collection(pc)
	
	# What does the histogram of errors look like without the threshhold? 
	sigma_fit = 0;# num_bins = 25
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
	
	# Plots 
	g, axes = plt.subplots(3,2)

	plots = [Temp_err,Temp_err_clip,Cold_err,Cold_err_clip,Warm_err,Warm_err_clip ]
	axes[1,0].set_ylabel("Cold Mass"); axes[2,0].set_ylabel("Warm Mass"); axes[0,0].set_ylabel("Temperature")

	g.suptitle("Chi Squared Confidence Intervals", fontsize=16)
	axes[0,0].set_title("Before Clipping")
	axes[0,1].set_title("After Clipping")


	def histplot(array,row,col):
		if row == 1 and col == 0:
			num_bins = 13
		elif row == 2 and col == 0: 
			num_bins = 10
		else:
			# Freedman Diaconis rule for number of bins
			num_bins = int((np.max(array) - np.min(array)) / (2 * iqr(array) * len(array)**(-1/3)))
		print(num_bins)
		# Plot histogram
		n, bins, patches = axes[row,col].hist(array,bins=num_bins,density=True)
		# Fit histogram
		(mu, sigma) = norm.fit(array)
		hist_y = mlab.normpdf(bins,mu,sigma)
		# Plot Fit
		axes[row,col].plot(bins, hist_y, '--')

	count = 0 
	for i in range(3):
		for j in range(2):
			histplot(plots[count],i,j)
			count += 1
plt.show()
