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
from mpl_colormap_tools import truncate_cmap

##############
## Switches ##
##############
# Plot Figures
plot = True
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
chiMap_holdColdMass = Eqs.CalculateChi(ColdTemp, cold_mass, WarmMass,kappa[0],Areas[0],AverageIntensities).transpose()
chiMap_holdTemp = Eqs.CalculateChi(temp, ColdMass, WarmMass,kappa[0],Areas[0],AverageIntensities)
FixedChiMaps = [chiMap_holdWarmMass,chiMap_holdColdMass,chiMap_holdTemp]

# Allow all parameters to vary, getting the total chi map cube.
chiCube = Eqs.CalculateChi(ColdTemp, ColdMass, WarmMass,kappa[0],Areas[0],AverageIntensities)

# 68.3%, 95.4% and 99.73%, or 1, 2 and 3 sigma 
intervals = [1,4,9]


## Find error intervals for some confidence level
def ConfidenceInterval(interval): 
    row, column, Z = np.where(chiCube < chi_squared + interval)
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

    # Stats 
    print("Temperature Elements in Interval: "+str((np.max(row)-np.min(row)+ self_width)))
    print("Cold Dust Mass Elements in Interval: "+str((np.max(column)-np.min(column)+ self_width)))
    print("Warm Dust Mass Elements in Interval: "+str((np.max(Z)-np.min(Z)+ self_width)))
    print("Temperature interval in K: "+str(Temperature_Confidence))
    print("Cold Dust Mass interval in Solar Mass Fraction: "+str(ColdMass_Confidence))
    print("Warm Dust Mass interval in Solar Mass Fraction: "+str(WarmMass_Confidence))

for i in range(3):
    print("...")
    print(str(np.sqrt(intervals[i]))+" sigma level stats:") 
    ConfidenceInterval(intervals[i])

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

    def ChiSquaredMap(Map,interval):
        # Change everything outside the interval to NaN.
        R,C = np.where(Map > chi_squared + interval)
        Map[R,C] = np.nan 
        return Map

    # This plots intervals while holding one of the solutions constant.
    # The intervals aren't calculated like this, they allow all parameters to vary
    # but plotting is easier this way.
    
    f, axes = plt.subplots(3,2)
    
    physical_X = [np.where(np.isclose(ColdMass,cold_mass))[0][0], np.where(np.isclose(WarmMass,warm_mass))[0][0], np.where(np.isclose(ColdMass,cold_mass))[0][0]]
    physical_Y = [np.where(np.isclose(ColdTemp,temp))[0][0],np.where(np.isclose(ColdTemp,temp))[0][0],np.where(np.isclose(WarmMass,warm_mass))[0][0]]
    
    titles = ["68.3% Confidence","95.4% Confidence"]
    ylabels = ["Temperature K","Warm Dust Mass M$_\odot$","Warm Dust Mass M$_\odot$"]
    xlabels = ["Cold Dust Mass M$_\odot$", "Temperature K", "Cold Dust Mass M$_\odot$"]
    
    for i in range(3):
        for j in range(2):
            axes[i,j].imshow(ChiSquaredMap(FixedChiMaps[i],intervals[j]))
            axes[i,j].scatter(physical_X[i],physical_Y[i],s=5,c='r')
            axes[i,j].set_ylim(physical_Y[i]-10,physical_Y[i]+10); axes[i,j].set_xlim(physical_X[i]-10,physical_X[i]+10);
            axes[i,j].set_xlabel(xlabels[i]); axes[i,j].set_ylabel(ylabels[i])
    axes[0,0].set_title(titles[0])
    axes[0,1].set_title(titles[1])            

plt.show()