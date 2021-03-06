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

# Running the pixel by pixel fit takes a while
# Run once and save in text file
calculate_data = False
# Load text files, print stats
load_data = True

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
                np.save("Sols/PixbyPix/ChiSquaredCubes/cube_"+str(i)+"_"+str(j)+"_.npy",chi_squared_cube)
                for k in range(2):
                    T,C,W = Eqs.ConfidenceInterval(intervals[k],chi_squared_cube,chi_squared,ColdTemp,ColdMass,WarmMass)
                    Temperature_Confidence[i,j,k] = T
                    ColdMass_Confidence[i,j,k] = C
                    WarmMass_Confidence[i,j,k] = W
   
    # # Save arrays
    np.save("Sols/PixbyPix/Total_SED.npy",allTotalSED); np.save("Sols/PixbyPix/Warm_SED.npy",allWarmSED); np.save("Sols/PixbyPix/Cold_SED.npy",allColdSED)            
    np.savetxt("Sols/PixbyPix/Temperature.txt", allTemps); np.savetxt("Sols/PixbyPix/ColdMass.txt",allColdMass); np.savetxt("Sols/PixbyPix/WarmMass.txt",allWarmMass)
    np.savetxt('Sols/PixbyPix/ChiSquared.txt',allChiSquared)    

    np.save("Sols/PixbyPix/Temperature_Confidence.npy",Temperature_Confidence)
    np.save("Sols/PixbyPix/Cold_Mass_Confidence.npy",ColdMass_Confidence)
    np.save("Sols/PixbyPix/Warm_Mass_Confidence.npy",WarmMass_Confidence)

if load_data:    
    ColdMass_Map = np.loadtxt("Sols/PixbyPix/ColdMass.txt"); WarmMass_Map = np.loadtxt("Sols/PixbyPix/WarmMass.txt")
    Temperature_Map = np.loadtxt("Sols/PixbyPix/Temperature.txt"); Chi_Map = np.loadtxt("Sols/PixbyPix/ChiSquared.txt")
    WarmSed_Cube = np.load("Sols/PixbyPix/Warm_SED.npy"); ColdSed_Cube = np.load("Sols/PixbyPix/Cold_SED.npy");
    Temperature_Confidence = np.load("Sols/PixbyPix/Temperature_Confidence.npy") 
    ColdMass_Confidence = np.load("Sols/PixbyPix/Cold_Mass_Confidence.npy")
    WarmMass_Confidence = np.load("Sols/PixbyPix/Warm_Mass_Confidence.npy")

    sky_noise = np.loadtxt("../Sky_Remove/Sigma.txt")
    res_conversion = np.loadtxt("../Sky_Remove/native_res_to_160_res.txt")

    # Confidence Interval
    # ========================================================
    sigma_fit = 1 # 0 is 1 sigma, 1 is 2 sigma, 2 is 3 sigma. 
    # ========================================================

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
    # ========================================================
    # Sets what multiple of the sky's noise to remove. 
    # ========================================================
    noise_multiplier = 3
    # ========================================================
    # Initialize array of empty 2d matrices the size of the images.
    empty_template = np.zeros(np.shape(DataCube[:,:,0])[0]*np.shape(DataCube[:,:,0])[1]).reshape(np.shape(DataCube[:,:,0])[0],np.shape(DataCube[:,:,0])[1])
    templates = [np.copy(empty_template),np.copy(empty_template),np.copy(empty_template),np.copy(empty_template)]

    for i in range(4):
        for j in range(np.shape(DataCube[:,:,i])[0]):
            for k in range(np.shape(DataCube[:,:,i])[1]):
                if np.isfinite(DataCube[:,:,i][j,k]):
                    # If the data is greater than some multiple of the noise, keep it by setting it to 1. 
                    if DataCube[:,:,i][j,k] > noise_multiplier * sky_noise[i] / res_conversion[i]: 
                        templates[i][j,k] = 1

    # Multiply all 4 templates together to get a final template             
    template = np.ones(np.shape(DataCube[:,:,0])[0]*np.shape(DataCube[:,:,0])[1]).reshape(np.shape(DataCube[:,:,0])[0],np.shape(DataCube[:,:,0])[1])
    for item in templates:
        template *= item    

    # Multiply it by Maps and then take Stats
    print("Total Mass "+str(np.sum(np.copy(ColdMass_Map)*template)+np.sum(np.copy(WarmMass_Map)*template))
            + " pm " + str(np.linalg.norm((np.copy(ColdMass_Confidence[:,:,sigma_fit])+np.copy(WarmMass_Confidence[:,:,sigma_fit]))*template,ord='fro')))
    print("Cold Mass "+str(np.sum(np.copy(ColdMass_Map)*template)) 
            + " pm " + str(np.linalg.norm(np.copy(ColdMass_Confidence[:,:,sigma_fit])*template,ord='fro')))
    print("Warm Mass "+str(np.sum(np.copy(WarmMass_Map)*template))
            + " pm " + str(np.linalg.norm(np.copy(WarmMass_Confidence[:,:,sigma_fit])*template,ord='fro')))
    
    print("Average Temperature "+str(np.mean((np.copy(Temperature_Map)*template)[np.where(Temperature_Map > 0)]))
            +" pm "+str(np.mean((np.copy(Temperature_Confidence[:,:,sigma_fit])*template)[np.where(Temperature_Map>0)])))