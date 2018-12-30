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
# Running the pixel by pixel fit takes a while
# Run once and save in text file
calculate_data = False
# Load text files
load_data = True
# Plot SEDS of 24 micron hot spot 
# and the 2 blobs at longer wavelengths
plot_selectedpix = True
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
# Dimensions of Data
shape_r, shape_c = np.shape(fits.open(files[0])[0].data)
# Chi Confidence Intervals
intervals = [1,4,9]

################################
#    Pixel by Pixel Fitting    #
################################
# Empty Data Cubes for Maps 
DataCube = np.zeros(shape_r*shape_c*4).reshape(shape_r,shape_c,4)
ChiSquaredIntervalCube = np.zeros(shape_r*shape_c*2).reshape(shape_r,shape_c,2)

# Get SNR data, and the physical pixel area 
pix_area = Eqs.pixarea(files[0])
SEDCube = np.zeros(shape_r*shape_c*len(lam)).reshape(shape_r,shape_c,len(lam))
shell = np.copy(DataCube[:,:,0])

for i in range(4):
    DataCube[:,:,i] = fits.open(files[i])[0].data 

if calculate_data:

    # Intiializing arrays to store all the solutions in. Reduces amount of files to be saved.
    allTemps = np.copy(shell); allColdMass = np.copy(shell); allWarmMass = np.copy(shell); allChiSquared = np.copy(shell);
    allTotalSED = np.copy(SEDCube); allWarmSED = np.copy(SEDCube); allColdSED = np.copy(SEDCube)
    Temperature_Confidence = np.copy(ChiSquaredIntervalCube) 
    ColdMass_Confidence = np.copy(ChiSquaredIntervalCube) 
    WarmMass_Confidence = np.copy(ChiSquaredIntervalCube) 

    # Look over each pixel, if it is not a nan, do the fit and save the results.
    count = 0; count_r, count_c = np.where(np.isfinite(DataCube[:,:,0]))
    for i in range(shape_r):
        for j in range(shape_c):
            if np.isfinite(DataCube[i,j,:]).all():
                count += 1
                print(("Calculating Pixel: [{},{}],{} Percent Complete").format(i,j,count/len(count_r)*100))
                total_sed, warm_sed, cold_sed, temp, cold_mass, warm_mass, chi_squared, chi_squared_cube = Eqs.CalculateBestSed(ColdTemp,ColdMass,WarmMass,kappa[0],pix_area,DataCube[i,j,:])
                allTotalSED[i,j,:] = total_sed; allWarmSED[i,j,:] = warm_sed; allColdSED[i,j,:] = cold_sed
                allTemps[i,j] = temp; allColdMass[i,j] = cold_mass; allWarmMass[i,j] = warm_mass; allChiSquared[i,j] = chi_squared
                # Saving the chi squared cube individually because of memory errors
                np.save("Sols/ChiSquaredCubes/cube_"+str(i)+"_"+str(j)+"_.npy",chi_squared_cube)
                for k in range(2):
                    T,C,W = Eqs.ConfidenceInterval(intervals[k],chi_squared_cube,chi_squared,ColdTemp,ColdMass,WarmMass)
                    Temperature_Confidence[i,j,k] = T
                    ColdMass_Confidence[i,j,k] = C
                    WarmMass_Confidence[i,j,k] = W
   
    # # Save arrays
    np.save("Sols/Total_SED.npy",allTotalSED); np.save("Sols/Warm_SED.npy",allWarmSED); np.save("Sols/Cold_SED.npy",allColdSED)            
    np.savetxt("Sols/Temperature.txt", allTemps); np.savetxt("Sols/ColdMass.txt",allColdMass); np.savetxt("Sols/WarmMass.txt",allWarmMass)
    np.savetxt('Sols/ChiSquared.txt',allChiSquared)    

    np.save("Sols/Temperature_Confidence.npy",Temperature_Confidence)
    np.save("Sols/Cold_Mass_Confidence.npy",ColdMass_Confidence)
    np.save("Sols/Warm_Mass_Confidence.npy",WarmMass_Confidence)

if load_data:    
    ColdMass_Map = np.loadtxt("Sols/ColdMass.txt"); WarmMass_Map = np.loadtxt("Sols/WarmMass.txt")
    Temperature_Map = np.loadtxt("Sols/Temperature.txt"); Chi_Map = np.loadtxt("Sols/ChiSquared.txt")
    WarmSed_Cube = np.load("Sols/Warm_SED.npy"); ColdSed_Cube = np.load("Sols/Cold_SED.npy");
    Temperature_Confidence = np.load("Sols/Temperature_Confidence.npy"); 
    ColdMass_Confidence = np.load("Sols/Cold_Mass_Confidence.npy",)
    WarmMass_Confidence = np.load("Sols/Warm_Mass_Confidence.npy")

    sigma_fit = 1
    # To get the error I'm calculating the frobenius norm of the interval map
    print("Total Mass "+str(np.sum(ColdMass_Map)+np.sum(WarmMass_Map))
        + " pm " + str(np.linalg.norm(ColdMass_Confidence[:,:,sigma_fit]+WarmMass_Confidence[:,:,sigma_fit],ord='fro')))
    print("Cold Mass "+str(np.sum(ColdMass_Map))
        + " pm " + str(np.linalg.norm(ColdMass_Confidence[:,:,sigma_fit],ord='fro')))
    print("Warm Mass "+str(np.sum(WarmMass_Map))
        + " pm " + str(np.linalg.norm(WarmMass_Confidence[:,:,sigma_fit],ord='fro')))
    count_r, count_c = np.where(np.isfinite(DataCube[:,:,0]))

    print("Average Temperature "+str(np.mean(Temperature_Map[np.where(Temperature_Map > 0)]))
            +" pm "+str(np.mean(Temperature_Confidence[:,:,sigma_fit][np.where(Temperature_Map>0)])))
    print("...")
    print("Detection Limited Stats")
    print("...")
    # Create Detection Limited Template
    n = 2
    Noise = np.loadtxt('../Sky_Remove/Sigma.txt')
    FlatData = np.sum(np.copy(DataCube),axis=2)
    FlatData[np.where(np.isnan(FlatData))] = 0
    Template = np.copy(FlatData)
    Template[np.where(np.copy(FlatData) > n * np.sum(Noise))] = 1
    Template[np.where(np.copy(FlatData) < n * np.sum(Noise))] = 0

    # Multiply it by Maps and then take Stats
    print("Total Mass "+str(np.sum(np.copy(ColdMass_Map)*Template)+np.sum(np.copy(WarmMass_Map)*Template))
            + " pm " + str(np.linalg.norm((np.copy(ColdMass_Confidence[:,:,sigma_fit])+np.copy(WarmMass_Confidence[:,:,sigma_fit]))*Template,ord='fro')))
    print("Cold Mass "+str(np.sum(np.copy(ColdMass_Map)*Template)) 
            + " pm " + str(np.linalg.norm(np.copy(ColdMass_Confidence[:,:,sigma_fit])*Template,ord='fro')))
    print("Warm Mass "+str(np.sum(np.copy(WarmMass_Map)*Template))
            + " pm " + str(np.linalg.norm(np.copy(WarmMass_Confidence[:,:,sigma_fit])*Template,ord='fro')))
    
    print("Average Temperature "+str(np.mean((np.copy(Temperature_Map)*Template)[np.where(Temperature_Map > 0)]))
            +" pm "+str(np.mean((np.copy(Temperature_Confidence[:,:,sigma_fit])*Template)[np.where(Temperature_Map>0)])))
#if plot_selectedpix:
    # First I want to look at what pix we're talking about 

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

#plt.show()