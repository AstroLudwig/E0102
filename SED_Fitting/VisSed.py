from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt 
from matplotlib.patches import Rectangle
import Eqs
import seaborn as sns
##############
## Switches ##
##############
# Chose whether to use both cold amc and cold mg2 
# or the detection limited cold amc sed images.
use_all = False
use_limit = True
use_noBkgdErr = False
use_PixLimit = False # Meant to be used with noBkgdErr.
save = False
# ~~~ Plotting ~~~ #
# Plot General SED from Average of Each image
plot_general = False
# Plot Pixel by Pixel SED. Interactive
plot_pix = False
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
plot_sigma = True
# Plot Contour Maps for Selected Pixels. 
plot_SelectContour = False
# Generate signal to noise maps without calibration or background error
generate_SN = False
# Plot SEDs For Warm Mass Peak and two blobs. 
plot_SelectSED = False
###########
## Files ##
###########

files = ['../Final_Files/24/Final_24_SNR_CR_Prnd.fits','../Final_Files/70/Final_70_SNR_CR_Prnd.fits','../Final_Files/100/Final_100_SNR_CR_Prnd.fits','../Final_Files/160/160um_70modeledFS_prnd_snr.fits']

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

lam, k_amc = Eqs.log_extrap('../Kappa/kappa_amc.dat')
lam, k_mg2 = Eqs.log_extrap('../Kappa/kappa_mg2sio4.dat')
wv = [24,70,100,160]
pix_area = Eqs.pixarea(files[0])
##############
# Get Data   #
##############
data = []; sed_means = []
for i in range(4):
	sm, ar = Eqs.sed_avg(files[i])
	data.append(fits.open(files[i])[0].data)
	sed_means.append(sm)

sigma = Eqs.error(sed_means)
if use_all:
	# SEDS
	# AMC Maps                                         # MG2 Maps
	AMC_Total_SED = np.load("Sols/Total_SED_Map.npy"); MG2_Total_SED = np.load("Sols/Total_SED_Map_MG2.npy")
	AMC_Warm_SED = np.load("Sols/Warm_SED_Map.npy");   MG2_Warm_SED = np.load("Sols/Warm_SED_Map_MG2.npy")
	AMC_Cold_SED = np.load("Sols/Cold_SED_Map.npy");   MG2_Cold_SED = np.load("Sols/Cold_SED_Map_MG2.npy")
	# Temperature
	# AMC Maps                                         # MG2 Maps
	AMC_Temp = np.loadtxt("Sols/TemperatureMap.txt");  MG2_Temp = np.loadtxt("Sols/MG2_TemperatureMap.txt")
	# Mass 
	# AMC Maps                                            # MG2 Maps
	AMC_Total_Mass = np.loadtxt("Sols/TotalMassMap.txt"); MG2_Total_Mass = np.loadtxt("Sols/MG2_TotalMassMap.txt")
	AMC_Warm_Mass = np.loadtxt("Sols/WarmMassMap.txt");   MG2_Warm_Mass = np.loadtxt("Sols/MG2_WarmMassMap.txt")
	AMC_Cold_Mass = np.loadtxt("Sols/ColdMassMap.txt");   MG2_Cold_Mass = np.loadtxt("Sols/MG2_ColdMassMap.txt")
	# Chi Squared Maps
	# AMC Maps                                   # MG2 Maps
	AMC_chi = np.loadtxt("Sols/ChiSqrdMap.txt"); MG2_chi = np.loadtxt("Sols/MG2_ChiSqrdMap.txt")
	# General SEDs 
	Gen_TSed = np.loadtxt("Sols/General_TotalSED.txt");Gen_WSed = np.loadtxt("Sols/General_WarmSED.txt")
	Gen_CSed = np.loadtxt("Sols/General_ColdSED.txt");Gen_Stats = np.loadtxt("Sols/General_TempMassWarmMassChi.txt")
	Gen_Mg2TSed = np.loadtxt("Sols/MG2_General_TotalSED.txt"); Gen_Mg2Stats = np.loadtxt("Sols/MG2_General_TempMassWarmMassChi.txt")
if use_limit:
	# SEDS
	# AMC Maps                                         
	AMC_Total_SED = np.load("Sols/DetectionLimited_Total_SED_Map.npy");
	AMC_Warm_SED = np.load("Sols/DetectionLimited_Warm_SED_Map.npy")
	AMC_Cold_SED = np.load("Sols/DetectionLimited_Cold_SED_Map.npy")
	# Temperature
	# AMC Maps                                         
	AMC_Temp = np.loadtxt("Sols/DetectionLimited_TemperatureMap.txt"); 
	# Mass 
	# AMC Maps                                           
	AMC_Total_Mass = np.loadtxt("Sols/DetectionLimited_TotalMassMap.txt"); 
	AMC_Warm_Mass = np.loadtxt("Sols/DetectionLimited_WarmMassMap.txt");
	AMC_Cold_Mass = np.loadtxt("Sols/DetectionLimited_ColdMassMap.txt");   
	# Chi Squared Maps
	# AMC Maps                                   
	AMC_chi = np.loadtxt("Sols/DetectionLimited_ChiSqrdMap.txt"); 

if use_noBkgdErr:
	# SEDS #BkgdErrRemovd_
	# AMC Maps                                         
	AMC_Total_SED = np.load("Sols/BkgdErrRemovd_Total_SED_Map.npy");
	AMC_Warm_SED =  np.load("Sols/BkgdErrRemovd_Warm_SED_Map.npy")
	AMC_Cold_SED =  np.load("Sols/BkgdErrRemovd_Cold_SED_Map.npy")
	# Temperature
	# AMC Maps                                         
	AMC_Temp = np.loadtxt("Sols/BkgdErrRemovd_TemperatureMap.txt"); 
	# Mass 
	# AMC Maps                                           
	AMC_Total_Mass = np.loadtxt("Sols/BkgdErrRemovd_TotalMassMap.txt"); 
	AMC_Warm_Mass = np.loadtxt("Sols/BkgdErrRemovd_WarmMassMap.txt");
	AMC_Cold_Mass = np.loadtxt("Sols/BkgdErrRemovd_ColdMassMap.txt");   
	# Chi Squared Maps
	# AMC Maps                                   
	AMC_chi = np.loadtxt("Sols/BkgdErrRemovd_ChiSqrdMap.txt"); 
	# General SEDs 
	Gen_TSed = np.loadtxt("Sols/BkgdErrRemovd_General_TotalSED.txt");Gen_WSed = np.loadtxt("Sols/BkgdErrRemovd_General_WarmSED.txt")
	Gen_CSed = np.loadtxt("Sols/BkgdErrRemovd_General_ColdSED.txt");Gen_Stats = np.loadtxt("Sols/BkgdErrRemovd_General_TempMassWarmMassChi.txt")
	
	if use_PixLimit:
		# Get a map of 1's and 0's that tells you which pixels to use. 
		template = np.copy(data[0])
		for i in range(np.shape(data[0])[0]):
			for j in range(np.shape(data[0])[1]):
				arr_int = np.asarray([data[0][i,j],data[1][i,j],data[2][i,j],data[3][i,j]])
				# If the mean of the intensities in each image is less than twice the mean error than remove it. 
				if np.isfinite(arr_int[0]) and np.isfinite(arr_int[1]) and np.isfinite(arr_int[2]) and np.isfinite(arr_int[3]):
					if np.nanmean(arr_int) < 3 * np.nanmean(Eqs.errorNoBkgd(arr_int)):
						template[i,j] = 0
					else: 
						template[i,j] = 1
		plt.imshow(template,vmin=0,vmax=1)
		plt.ylim(110,130); plt.xlim(118,140)
		plt.title("Template Image")
		# Multiply each SED by this template. 
		for i in range(np.shape(AMC_Total_SED)[2]):
			AMC_Total_SED[:,:,i] = AMC_Total_SED[:,:,i] * template
			AMC_Warm_SED[:,:,i] = AMC_Warm_SED[:,:,i] * template  
			AMC_Cold_SED[:,:,i] = AMC_Cold_SED[:,:,i] * template    

##########
# Plot   #
##########

if plot_general:

	fig = plt.figure(figsize=(11,8))

	plot = fig.add_subplot(111)
	#plt.style.use('bmh')
	#plot.style.use('seaborn-darkgrid')
	#plt.style.use('seaborn')
	
	plot.plot(lam,Gen_TSed,color="#424186")#"#440154")
	plot.plot(lam,Gen_WSed,color="#84D44B",ls='dashed') ##fde724
	plot.plot(lam,Gen_CSed,color="#23A883")
	if use_all:
		plot.errorbar(wv,sed_means,yerr=Eqs.error(sed_means),marker='o',linestyle='none',color="black")#"#84D44B")
	elif use_noBkgdErr:	
		plot.errorbar(wv,sed_means,yerr=Eqs.errorNoBkgd(sed_means),marker='o',linestyle='none',color="black")#"#84D44B")
	plot.set_xlabel("Wavelength ($\mu m$)",size=18)
	plot.set_ylabel("Spectral Intensity (Mjy sr$^{-1}$)",size=18)
	plot.set_title(("Average Spectral Energy Distribution (SED)\n").format(int(Gen_Stats[0]),Gen_Stats[1],Gen_Stats[2],Gen_Stats[3]),size=20)
	plot.legend(("Total SED","Warm SED","Cold SED"),prop={'size':14})
	
	plot.tick_params(axis='both', which='major', labelsize=16)
	plot.tick_params(axis='both', which='minor', labelsize=14)
	
	plot.grid(color='white',linestyle='-')
	plot.set_facecolor("#EAEAF2")

	#[t,Cm,wm,chi]
	print(("Temp {} Cold Mass {} Warm Mass {} Total Mass {} ChiSqrd {} ").format(int(Gen_Stats[0]),Gen_Stats[1],Gen_Stats[2],Gen_Stats[1]+Gen_Stats[2],Gen_Stats[3]))
	if save:
		plt.savefig("AverageSED.png")
if plot_pix:

	f, (ax,bx) = plt.subplots(1,2)
	inObsvSed = [data[0][0,0],data[1][0,0],data[2][0,0],data[3][0,0]]
	inObsvErr = Eqs.error(inObsvSed)
	ax.imshow(data[3])
	ax.scatter(0,0,c='r',marker='s',s=40)
	bx.plot(lam,AMC_Total_SED[0,0],label="Total SED")
	bx.plot(lam,AMC_Warm_SED[0,0],label="Warm SED")
	bx.plot(lam,AMC_Cold_SED[0,0],label="Cold SED")
	bx.errorbar(wv,inObsvSed,yerr=inObsvErr,marker='o',linestyle='none',c='purple')
	ax.set_xlim(xdim[0],xdim[1])
	ax.set_ylim(ydim[0],ydim[1])
	bx.set_ylim(0, np.max(AMC_Total_SED[0,0])+2)
	bx.grid()
	plt.legend()
	plt.title(("Temperature: {:2d} K,  Cold Mass: {:.2f} $M_\odot$,  Warm Mass: {:.2E} $M_\odot$,  $\chi^2$: {:.2f}").format(int(AMC_Temp[0,0]),AMC_Cold_Mass[0,0],AMC_Warm_Mass[0,0],AMC_chi[0,0]))
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

		bx.plot(lam,AMC_Total_SED[i,j],label="Total SED")
		bx.plot(lam,AMC_Warm_SED[i,j],label="Warm SED")
		bx.plot(lam,AMC_Cold_SED[i,j],label="Cold SED")
		
		bx.set_ylim(0, np.max(AMC_Total_SED[i,j])+2)
		bx.errorbar(wv,ObsvSed,yerr=ObsvErr,marker='o',linestyle='none',c='purple')
		bx.grid()
		plt.legend()
		plt.title(("Temperature: {:2d} K,  Cold Mass: {:.2E} $M_\odot$,  Warm Mass: {:.2E} $M_\odot$,  $\chi^2$: {:.2f}").format(int(AMC_Temp[i,j]),AMC_Cold_Mass[i,j],AMC_Warm_Mass[i,j],AMC_chi[i,j]))
		bx.set_xlabel("Wavelength ($\mu$m)")
		bx.set_ylabel("Spectral Intensity (MJy sr$^{-1}$)")

		plt.draw()

	cid = f.canvas.mpl_connect('button_press_event', onclick)

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


	bx.set_title(("Temperature: {:2d} K,  Cold Mass: {:.2E} $M_\odot$,  Warm Mass: {:.2E} $M_\odot$,  $\chi^2$: {:.2f}").format(int(AMC_Temp[0,0]),AMC_Cold_Mass[0,0],AMC_Warm_Mass[0,0],AMC_chi[0,0]))
	


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
		bx.set_title(("Temperature: {:2d} K,  Cold Mass: {:.2E} $M_\odot$,  $\chi^2$: {:.2f}").format(int(AMC_Temp[i,j]),AMC_Cold_Mass[i,j],AMC_chi[i,j]))
		#cx.imshow(np.log(chiMap[i,j]))
		#cb.remove()

		plt.draw()

		
	cid = f.canvas.mpl_connect('button_press_event', onclick)	

if plot_imgs:    
	f, axes = plt.subplots(2,2,figsize=(10,15))
	x1 = 115; x2 = 140; y1 = 110; y2 = 130
	vmin = 0; vmax = 10
	axes[0,0].imshow(fits.open(files[0])[0].data,vmin=vmin,vmax=vmax)
	axes[0,1].imshow(fits.open(files[1])[0].data,vmin=vmin,vmax=vmax)
	axes[1,0].imshow(fits.open(files[2])[0].data,vmin=vmin,vmax=vmax)
	axes[1,1].imshow(fits.open(files[3])[0].data,vmin=vmin,vmax=vmax)
	axes[0,0].set_xlim(x1,x2)
	axes[0,1].set_xlim(x1,x2)
	axes[1,0].set_xlim(x1,x2)
	axes[1,1].set_xlim(x1,x2)
	axes[0,0].set_ylim(y1,y2)
	axes[0,1].set_ylim(y1,y2)
	axes[1,0].set_ylim(y1,y2)
	axes[1,1].set_ylim(y1,y2)  

if plot_maps: # Temperature Map
	# Prep Temp Map for Mean Measurement
	a,b = np.shape(AMC_Temp)
	temp = np.copy(AMC_Temp)
	for i in range(a):
		for j in range(b):
			if AMC_Temp[i,j] == 0.:
				temp[i,j] = np.nan
				
	g, (ax,bx,cx) = plt.subplots(1,3)
	tim = ax.imshow(AMC_Temp ,vmin = np.nanmin(temp), vmax = np.nanmax(temp))
	ax.set_title(("Temperature Map.\n Average is {}").format(int(np.nanmean(temp))))
	ax.set_ylim(110,130)
	ax.set_xlim(115,140)
	ax.axis("off")
	g.colorbar(tim, ax=ax)

	wim = bx.imshow(AMC_Warm_Mass ,vmin = np.min(AMC_Warm_Mass), vmax = np.max(AMC_Warm_Mass))
	bx.set_title(("Warm Mass Map.\n Sum is {:2f}").format(np.sum(AMC_Warm_Mass)))
	bx.set_ylim(110,130)
	bx.set_xlim(115,140)
	bx.axis("off")
	g.colorbar(wim, ax=bx)

	cim = cx.imshow(AMC_Cold_Mass ,vmin = np.min(AMC_Cold_Mass), vmax = np.max(AMC_Cold_Mass))
	cx.set_title(("Cold Mass Map.\n Sum is {:2f}").format(np.sum(AMC_Cold_Mass)))
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
	ax.set_title(("Temperature: {:2d} K,  Cold Mass: {:.2f} $M_\odot$, \n Warm Mass: {:.2E} $M_\odot$,  $\chi^2$: {:.2f}").format(int(AMC_Temp[0,0]),AMC_Cold_Mass[0,0],AMC_Warm_Mass[0,0],AMC_chi[0,0]))

	# SED Plot 
	bx.plot(lam,AMC_Total_SED[0,0],label="Total SED")
	bx.plot(lam,AMC_Warm_SED[0,0],label="Warm SED")
	bx.plot(lam,AMC_Cold_SED[0,0],label="Cold SED")
	bx.errorbar(wv,inObsvSed,yerr=inObsvErr,marker='o',linestyle='none',c='purple')
	bx.set_ylim(0, np.max(AMC_Total_SED[0,0])+2)
	bx.grid()
	bx.set_xlabel("Wavelength ($\mu$m)")
	bx.set_ylabel("Spectral Intensity (MJy sr$^{-1}$)")
	
	chiPlot = np.log(chiMap[0,0].transpose())
	#chiPlot = chiMap[0,0].transpose()
	cs = cx.contourf(M,T,chiPlot,cmap=cmap,levels=lvls)
	cb = f.colorbar(cs)
	cx.set_xlabel("Mass Range in Solar Fractional Mass")
	cx.set_ylabel("Temperature in Kelvin")
	xr = [M_[0],AMC_Cold_Mass[0,0]+2*AMC_Cold_Mass[0,0]]
	yr = [AMC_Temp[0,0]-.5*AMC_Temp[0,0],AMC_Temp[0,0]+.5*AMC_Temp[0,0]]
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

		bx.plot(lam,AMC_Total_SED[i,j],label="Total SED")
		bx.plot(lam,AMC_Warm_SED[i,j],label="Warm SED")
		bx.plot(lam,AMC_Cold_SED[i,j],label="Cold SED")
		
		bx.set_ylim(0, np.max(AMC_Total_SED[i,j])+2)
		bx.errorbar(wv,ObsvSed,yerr=ObsvErr,marker='o',linestyle='none',c='purple')
		bx.grid()
		#plt.legend()
		ax.set_title(("Temperature: {:2d} K,  Cold Mass: {:.2E} $M_\odot$, \n Warm Mass: {:.2E} $M_\odot$,  $\chi^2$: {:.2f}").format(int(AMC_Temp[i,j]),AMC_Cold_Mass[i,j],AMC_Warm_Mass[i,j],AMC_chi[i,j]))
		bx.set_xlabel("Wavelength ($\mu$m)")
		bx.set_ylabel("Spectral Intensity (MJy sr$^{-1}$)")

		
		#chiPlot = chiMap[i,j].transpose()
		cs = cx.contourf(M,T,chiPlot,cmap=cmap,levels=lvls)
		#cb = f.colorbar(cs)
		cx.set_xlabel("Mass Range in Solar Fractional Mass")
		cx.set_ylabel("Temperature in Kelvin")
		#cx.grid()
		xr = [M_[0],2*AMC_Cold_Mass[i,j]]
		yr = [AMC_Temp[i,j]-.5*AMC_Temp[i,j],AMC_Temp[i,j]+.5*AMC_Temp[i,j]]
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
	#AMC_Warm_Mass = np.loadtxt("Sols/BkgdErrRemovd_WarmMassMap.txt");
	#AMC_Cold_Mass = np.loadtxt("Sols/BkgdErrRemovd_ColdMassMap.txt");   

	TotalMass = np.copy(AMC_Total_Mass)
	WarmMass = np.copy(AMC_Warm_Mass)
	ColdMass = np.copy(AMC_Cold_Mass)
	TemP = np.copy(AMC_Temp)
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
		plots[i].plot(lam,AMC_Total_SED[Py[i],Px[i]],label="Total SED",color="#424186")
		plots[i].plot(lam,AMC_Warm_SED[Py[i],Px[i]],label="Warm SED",color="#84D44B",ls='dashed')
		plots[i].plot(lam,AMC_Cold_SED[Py[i],Px[i]],label="Cold SED",color="#23A883")
		
		plots[i].set_xlabel("Wavelength ($\mu m$)")
		plots[i].set_ylabel("Spectral Intensity (MJy sr$^{-1}$)")
	
		plots[i].errorbar(wv,ObsvSed,yerr=ObsvErr,marker='o',linestyle='none',c='black')
		plots[i].grid(color='white',linestyle='-')
		plots[i].set_facecolor("#EAEAF2")
		plt.legend()
	
	

	
	

plt.show()
