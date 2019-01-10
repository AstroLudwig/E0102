# -*- Copyright (c) 2018, Bethany Ann Ludwig, All rights reserved. -*-
"""
NAME:
    Averaged SED
PURPOSE:
    Average the pixel intensities, use the full area of the remnant, fit an SED to the data.
"""
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
>>>>>>> Stashed changes
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
total_sed, warm_sed, cold_sed, temp, cold_mass, warm_mass, chi_squared, chi_squared_cube = Eqs.CalculateBestSed(ColdTemp,ColdMass,WarmMass,kappa[0],Areas[0],AverageIntensities)
print(len(total_sed))
# Print Solutions
print("*~~Solutions~~~*")
print("Temperature: " + str(temp)+ " K")
print("Cold Mass: " + str(cold_mass)+" Solar Mass")
print("Warm Mass: " + str(warm_mass)+" Solar Mass")
print("*~~~~~~~~~~~~~~*")

##############################
#  Chi^2 Confidence Interval #
##############################

# 68.3%, 95.4% and 99.73%, or 1, 2 and 3 sigma 
intervals = [1,4,9]

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

    count_r, count_c = np.where(np.isfinite(DataCube[:,:,0]))


for i in range(3):
    print("...")
    print(str(np.sqrt(intervals[i]))+" sigma level stats:") 
    Eqs.ConfidenceInterval(intervals[i],chi_squared_cube,chi_squared,ColdTemp,ColdMass,WarmMass)

##################
# Save Solutions #
##################

np.savetxt("Sols/Integrated/cold_sed.txt",cold_sed)
np.savetxt("Sols/Integrated/warm_sed.txt",warm_sed)
np.savetxt("Sols/Integrated/temp_coldmass_warmmass_chisqrd.txt",[temp,cold_mass,warm_mass,chi_squared])
np.save("Sols/Integrated/chi_squared_cube.npy",chi_squared_cube)
