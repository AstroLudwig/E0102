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
plot_ChiSquaredConfidence = True
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

##############################
#  Chi^2 Confidence Interval #
##############################

# Get the chi map holding one of the solutions constant.
chiMap_holdWarmMass = Eqs.CalculateChi(ColdTemp, ColdMass, warm_mass,kappa[0],Areas[0],AverageIntensities)
chiMap_holdColdMass = Eqs.CalculateChi(ColdTemp, cold_mass, WarmMass,kappa[0],Areas[0],AverageIntensities)
chiMap_holdTemp = Eqs.CalculateChi(temp, ColdMass, WarmMass,kappa[0],Areas[0],AverageIntensities)
FixedChiMaps = [chiMap_holdWarmMass,chiMap_holdColdMass,chiMap_holdTemp]
# Allow all parameters to vary, getting the total chi map cube.
chiCube = Eqs.CalculateChi(ColdTemp, ColdMass, WarmMass,kappa[0],Areas[0],AverageIntensities)

# 90% and 63%
intervals = [2.71,2.3]


## Find error intervals for 63% confidence
row, column, Z = np.where(chiCube < chi_squared + intervals[1])
self_width = 1
# Measure chi width's based on how far from the minimum it can be. 
# Max element - Min element + 1 (since there shouldn't be a 0 width, even 1 pixel has some width)
# This is multiplied by the interval in parameter space and divided by 2 to give a plus minus error.
# The interval is arbitrary so the array element values don't matter. Just want the width of a single parameter.
# Temperature interval
Temperature_Confidence = (np.max(row)-np.min(row)+ self_width) * (ColdTemp[50] - ColdTemp[49]) / 2
# Cold Mass Interval
ColdMass_Confidence = (np.max(column)-np.min(column)+ self_width) * (ColdMass[30] - ColdMass[29]) / 2
# Warm Mass Interval
WarmMass_Confidence = (np.max(Z)-np.min(Z)+ self_width) * (WarmMass[30] - WarmMass[29]) / 2

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
    f, axes = plt.subplots(3,2)
    titles = ["90% Confidence","68% Confidence"]
    # 90% -> 2.71, 68% -> 2.3
    for i in range(3):
        for j in range(2):
            axes[i,j].imshow(FixedChiMaps[i],vmin=0,vmax=chi_squared+intervals[j])
            axes[i,j].scatter(np.where(np.isclose(ColdMass,cold_mass))[0][0],np.where(np.isclose(ColdTemp,temp))[0][0],s=5,c='r')
"""
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
"""
plt.show()