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
plot_integrated = True
# Plot Chi Squared Confidence Intervals for Integrated SED 
plot_integrated_ChiSquaredConfidence = True


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
plt.show()
