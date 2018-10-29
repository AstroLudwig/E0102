# -*- Copyright (c) 2018, Bethany Ann Ludwig, All rights reserved. -*-
"""
NAME:
    Averaged SED
PURPOSE:
    Average the pixel intensities, use the full area of the remnant, fit an SED to the data.
"""
from timeit import default_timer as timer
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt 
import Eqs


##############
## Switches ##
##############
# Plot Figures
plot = False
plot_ChiSquaredConfidence = False
# Save Figures
save = False
###########
## Files ##
###########
# Insert files to get seds.
files = ['../Final_Files/24/24_SNR_Convolve_Regrid_Prune.fits',
         '../Final_Files/70/70_SNR_Convolve_Regrid_Prune.fits',
         '../Final_Files/100/100_SNR_Convolve_Regrid_Prune.fits',
         '../Final_Files/160/160_SNR_Prune.fits']

################
# Kappa Ranges #
################
# Extrapolate for entire kappa length.
kappa = ['Kappa/kappa_amc.dat','Kappa/kappa_mg2sio4.dat']
lam, k_amc = Eqs.log_extrap('Kappa/kappa_amc.dat')
lam, k_mg2 = Eqs.log_extrap('Kappa/kappa_mg2sio4.dat')

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

#################
#    Fitting    #
#################
# Get SNR data, the physical area of each SNR, and the average pixel intensity
data = []; Areas = []; AverageIntensities = np.zeros(4)
for i in range(4):
    AverageIntensity, Area = Eqs.AverageSED(files[i])
    data.append(fits.open(files[i])[0].data)
    Areas.append(Area)
    AverageIntensities[i] = AverageIntensity
# Get the error associated with each intensity
AverageError = Eqs.error(AverageIntensities)

 # Do the fit
total_sed, warm_sed, cold_sed, temp, cold_mass, warm_mass, chi_squared = Eqs.CalculateBestSed(ColdTemp,ColdMass,WarmMass,kappa[0],Areas[0],AverageIntensities)
# Get the chi map
chiMap = Eqs.CalculateChi(ColdTemp, ColdMass, warm_mass,kappa[0],Areas[0],AverageIntensities)
##################
#    Plotting    #
##################
if plot:
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

if plot_ChiSquaredConfidence:
    # This plots intervals while holding the warm mass solution constant.
    # Can create chi intervals without this requirement, but is trickier to plot.
    f, axes = plt.subplots(1,3)
    intervals = [6.63,2.71,2.3]; titles = ["99% Confidence","90% Confidence","68% Confidence"]
    # 99 % -> 6.63, 90% -> 2.71, 68% -> 2.3
    for i in range(3):
        axes[i].imshow(chiMap,vmin=0,vmax=chi_squared+intervals[i])
        axes[i].scatter(np.where(np.isclose(ColdMass,cold_mass))[0][0],np.where(np.isclose(ColdTemp,temp))[0][0],s=5,c='r')

        axes[i].set_xlabel("Cold Dust Mass")
        axes[i].set_ylabel("Temperature")
        axes[i].set_title(titles[i])
        # Fix some wonkiness with the axis labels
        Yticks = [0,5,10,15,20,25,30,35,40,45,50,55,60,65]
        Xticks = [0,5,10,15,20,25,30,35,40]
        axes[i].set_xticks(Xticks)
        axes[i].set_xticklabels(ColdMass[Xticks])
        axes[i].set_yticks(Yticks)
        axes[i].set_yticklabels(ColdTemp[Yticks])
        axes[i].set_ylim(17,47)
        axes[i].set_xlim(20,40)

plt.show()