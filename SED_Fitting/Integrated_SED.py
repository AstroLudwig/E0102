# -*- Copyright (c) 2018, Bethany Ann Ludwig, All rights reserved. -*-
"""
NAME:
    Integrated SED
PURPOSE:
    Fit each pixel intensity, use the pixel area, fit an SED to the data.
"""
from timeit import default_timer as timer
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt 
import Eqs


##############
## Switches ##
##############

# Fit the entire remnant rather than pixel by pixel.
plot_general = True


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

#########################
#    General Fitting    #
#########################
# Get SNR data, the physical area of each SNR, and the average pixel intensity
data = []; Areas = []; AverageIntensities = np.zeros(4)
for i in range(4):
    AverageIntensity, Area = Eqs.AverageSED(files[i])
    data.append(fits.open(files[i])[0].data)
    Areas.append(Area)
    AverageIntensities[i] = AverageIntensity
# Get the error associated with each intensity
AverageError = Eqs.error(AverageIntensities)

  
total_sed, warm_sed, cold_sed, temp, cold_mass, warm_mass, chi_squared = Eqs.CalculateBestSed(ColdTemp,ColdMass,WarmMass,kappa[0],Areas[0],AverageIntensities)
chiMap = Eqs.CalculateChi(ColdTemp, ColdMass, WarmMass,kappa[0],Areas[0],AverageIntensities)

if plot_general:
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

plt.show()