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
