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

# Save Figures
save = False
calculate_data = True
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

################################
#    Pixel by Pixel Fitting    #
################################
if calculate_data:
    # Get SNR data, and the physical pixel area in each SNR
    pix_area = Eqs.pixarea(files[0])

    shape_r, shape_c = np.shape(fits.open(files[0])[0].data)
    DataCube = np.zeros(shape_r*shape_c*4).reshape(shape_r,shape_c,4)
    SEDCube = np.zeros(shape_r*shape_c*len(lam)).reshape(shape_r,shape_c,len(lam))
    shell = np.copy(DataCube[:,:,0])
    for i in range(4):
        DataCube[:,:,i] = fits.open(files[i])[0].data 

    # Intiializing arrays to store all the solutions in. Reduces amount of files to be saved.
    allTemps = np.copy(shell); allColdMass = np.copy(shell); allWarmMass = np.copy(shell); allChiSquared = np.copy(shell);
    allTotalSED = np.copy(SEDCube); allWarmSED = np.copy(SEDCube); allColdSED = np.copy(SEDCube)
    allChiSquaredCubes = []; blankChiSqrdCube = np.zeros(len(ColdTemp)*len(ColdMass)*len(WarmMass)).reshape(len(WarmMass),len(ColdTemp),len(ColdMass))
    # Look over each pixel, if it is not a nan, do the fit and save the results.
    for i in range(shape_r):
        for j in range(shape_c):
            if np.isfinite(DataCube[i,j,:]).all():
                total_sed, warm_sed, cold_sed, temp, cold_mass, warm_mass, chi_squared, chi_squared_cube = Eqs.CalculateBestSed(ColdTemp,ColdMass,WarmMass,kappa[0],pix_area,DataCube[i,j,:])
                allTotalSED[i,j,:] = total_sed; allWarmSED[i,j,:] = warm_sed; allColdSED[i,j,:] = cold_sed
                allTemps[i,j] = temp; allColdMass[i,j] = cold_mass; allWarmMass[i,j] = warm_mass; allChiSquared[i,j] = chi_squared
                allChiSquaredCubes.append(chi_squared_cube)
            else:
                 allChiSquaredCubes.append(blankChiSqrdCube)               
    # Save arrays
    np.save("Sols/Total_SED.npy",allTotalSED); np.save("Sols/Warm_SED.npy",allWarmSED); np.save("Sols/Cold_SED.npy",allColdSED)            
    np.savetxt("Sols/Temperature.txt", allTemps); np.savetxt("Sols/ColdMass.txt",allColdMass); np.savetxt("Sols/WarmMass.txt",allWarmMass)
    np.save("Sols/ChiSquaredCubes.npy",allChiSquaredCubes); np.savetxt('Sols/ChiSquared.txt',allChiSquared)    
# # Get the error associated with each intensity
# AverageError = Eqs.error(AverageIntensities)

  
# 
# chiMap = Eqs.CalculateChi(ColdTemp, ColdMass, WarmMass,kappa[0],Areas[0],AverageIntensities)

# if plot_general:
#     fig = plt.figure(figsize=(11,8))
#     plot = fig.add_subplot(111)
#     # Plot SEDs
#     plot.plot(lam,total_sed,color="#424186")
#     plot.plot(lam,warm_sed,color="#84D44B",ls='dashed') 
#     plot.plot(lam,cold_sed,color="#23A883")
#     # Plot Measured Values
#     plot.errorbar([24,70,100,160],AverageIntensities,yerr=AverageError,marker='o',linestyle='none',color="black")
#     # Plot Labels
#     plot.set_xlabel("Wavelength ($\mu m$)",size=18)
#     plot.set_ylabel("Spectral Intensity (Mjy sr$^{-1}$)",size=18)
#     plot.set_title("Average Spectral Energy Distribution",size=20)
#     plot.legend(("Total SED","Warm SED","Cold SED"),prop={'size':14})
    
#     plot.tick_params(axis='both', which='major', labelsize=16)
#     plot.tick_params(axis='both', which='minor', labelsize=14)
    
#     plot.grid(color='white',linestyle='-')
#     plot.set_facecolor("#EAEAF2")

#     print(("Temp {} Cold Mass {} Warm Mass {} Total Mass {} Chi Squared {} ").format(temp,cold_mass,warm_mass,cold_mass+warm_mass,chi_squared))
#     if save:
#         plt.savefig("AverageSED.png")

plt.show()