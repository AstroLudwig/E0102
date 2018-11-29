# -*- Copyright (c) 2018, Bethany Ann Ludwig, All rights reserved. -*-
"""
NAME:
    Equations
PURPOSE:
    Modules for SED and VisSED: 
        Planck Equation 
        Modified Single Temperature Black Body Equation,
        Extrapolate Kappa 
        Physical Pixel Area
        Error 
        General SED
        Integrated SED
        Chi Squared Maps
"""
import numpy as np
from astropy import units as u
from astropy.io import fits
from astropy import constants as const
from scipy.optimize import curve_fit
import numexpr as ne

test = False

######################
# Extrapolate Kappa  #
######################

# The first array is wavelength in microns (10^-6 meters) 
# The second array is kappa in cm^2/g.
def log_extrap(kappafile):
# Admin
    data = np.loadtxt(kappafile)
    kappa = data[:,1]
    lam = data[:,0]

# Define straight line in log space. 
    def pwrfn(x,a,b):
        return a * np.power(x,b)
# Choose a range to model off of, this is 580th element to the end.   
    popt,pcov = curve_fit(pwrfn,lam[580:len(lam)],kappa[580:len(lam)])

# Extrapolate
    lam_extend = np.arange(103.9,164,0.1)
    newk = np.log(popt[0]) + popt[1] * np.log(lam_extend)
    newk_linear = np.exp(newk)
# Merge Arrays 
    full_lam = np.append(lam,lam_extend)
    full_kappa = np.append(kappa,newk_linear)

    return full_lam,full_kappa    

####################
# Planck Function  #
####################

# Inputs wavelength and temperature without units.
# Outputs spectral intensity in MJy/Sr

def B_nu(wav,temp): 
    # Constants
    h_ = const.h.cgs; c_ = const.c.cgs; kb_ = const.k_B.cgs
    nu = (c_ / (wav*u.micron).to(u.cm))
    temp = temp * u.K
    
    # Plank Equation 
    first = (2 * h_ * nu * nu * nu) / (c_ * c_)
    second = 1 / (ne.evaluate('exp((h_ * nu)/(kb_ * temp))')  -  1)

    # Dividing by steradians and hz*s isn't actually changing any units.
    # It makes converting to Mjy/sr possible. 
    first = first / u.sr / u.Hz / u.s
    ans  = (first * second).to(u.MJy/u.sr)
    
    return ans
#################################
# Modified Black Body Function  #
#################################

# Modified Single Temperature Black Body equation

def ModBB(kappa_arr,area,Temp,Mass):
    M_sun = const.M_sun.cgs
    Solutions = []
    lam, kappa = log_extrap(kappa_arr)
    for i in range(len(lam)):
        AbsorptionCrossSection = kappa[i] * u.cm * u.cm / u.g
        SigmaDust = Mass * M_sun / area
        Planck = B_nu(lam[i], Temp)
        Answer = AbsorptionCrossSection * SigmaDust * Planck 
        Solutions.append((Answer).value)

        # Flag if the units are off
        if Answer.unit != "MJy/sr":
            print(("Warning, ModBB units = {}").format(ans.unit))
         

    return np.asarray(Solutions)

################
# Pixel Area   #
################

#Calculate the pixel area using cdelt from the header and x = theta * d
#where d is the distance to the smc (61 kpc)
#Area is then radian^2 (or steradian) * distance_to_the_smc^2

# Units are in sq parsec

def pixarea(filename): 
    SMC_d = 61 * np.power(10.,3) * u.pc #Parsecs                       
    hdr = fits.open(filename)[0].header
    # Assuming Pixel is Square, if cdelts are different 
    # this method wouldn't work. 
    cdelt1 = np.abs(hdr["CDELT1"]) 
    cdelt2 = np.abs(hdr["CDELT2"])
    cd1 = (cdelt1 * u.deg).to(u.rad)
    cd2 = (cdelt2 * u.deg).to(u.rad)
    sr = (cd1 * cd2) # Steradian
    area = (sr * np.power(SMC_d,2) / u.rad**2).to(u.cm*u.cm)
    
    return area


###########
# Error   #
###########

def error(vals):
    # Standard Deviation in the Region chosen for Sky Removal
    Noise = np.loadtxt('../Sky_Remove/Sigma.txt')
    # Calibration Error, may need to update this?
    Calibr = np.repeat(.1,4)

    # Initiate blank error array   
    err = np.zeros(4)
    # Fill array with appropriate errors.
    for i in range(4):
        X = vals[i]
        quad = (Noise[i])**2 + (Calibr[i]*X)**2 
        err[i] = (np.sqrt(quad))
    return np.asarray(err)


####################################
# Average Pixel Intensity and Area #
####################################

# The measured SED is just an average of the intensities in the SNR
# This returns that average and the physical area of the remnant.  ## sed_avg previous name
def AverageSED(imgfile): 
    img = fits.open(imgfile)[0].data
    look = np.where(np.isfinite(img))
    count = np.shape(look)[1]
    avg = np.nansum(img) / count 
    pix_area = pixarea(imgfile)     
    area = (pix_area * count).to(u.cm*u.cm)
    return avg,area 
########################################
# Locate Observed Wavelength Features  #
########################################
# If you need indexes where observed wavelengths match lambda
# This isn't currently being called but is left in the code in
# Case lambda needs to vary in the future.
def GetObservationIndex(Lambda):
    ObservedWavelengths = [24,70,100,160]; index = []
    for i in range(4):
        index.append(np.where(np.isclose(Lambda,ObservedWavelengths[i]))[0][0])
    return index    
#############
# Fit SED   #
#############
"""
NAME:
   CalculateBestSed
PURPOSE:
   Generate SEDs given solutions for various Masses and Temperatures.
INPUT:
   coldTemp= 1D Array of Parameters
   warmTemp= 145 K in most cases (Karin, 2009)
   coldMass= 1D Array of Parameters
   warmMass= 1D Array of Parameters
   coldKappaFile= String for which kappa composition file you want to fit the cold dust to. 
   area= Either the pixel area or the total area for integrated or averaged solutions respectively
   measured_sed= A list or array of 4 values of averaged or individual pixel intensities at each wavelength
OUTPUT: SED Solutions
"""
def CalculateBestSed(coldTemp, coldMass, warmMass, coldKappaFile, area,measured_sed):
    # Constants
    warmTemp = 145 # K
    warmKappaFile = 'Kappa/kappa_mg2sio4.dat'
    lam, warmKappa = log_extrap(warmKappaFile)
    # Create Parameter Grid
    coldTempv, coldMassv, warmMassv = np.meshgrid(coldTemp, coldMass, warmMass) 
    # Get the error based on input measured sed
    sigma = error(measured_sed)
    # Solve for Cold Dust Component
    coldComponent = np.asarray(ModBB(coldKappaFile,area,coldTempv,coldMassv)).transpose()
    # Fixing warm component based on evidence from Sandstrom 2009
    warmComponent = np.asarray(ModBB(warmKappaFile,area,warmTemp,warmMassv)).transpose()
    # Put warm and cold together 
    totalSED = coldComponent + warmComponent
    # The numbers are locations that match up to 24, 70, 100, 160
    # For information on how they are calculated see GetObservationIndex
    index = GetObservationIndex(lam)
    totalSED = totalSED[:,:,:,index]
    
    # Calculate the Error and Determine where it is Lowest.
    chiSquare = (measured_sed - totalSED)**2 / (sigma)**2
    # Create Chi Squared Cube for measuring confidence intervals
    chiSquareCube = np.sum(chiSquare, axis=3)
    bestError = np.min(chiSquareCube)
    bestIndices = np.where(chiSquareCube == bestError)

    # Get Solutions
    best_cold_temp = coldTemp[bestIndices[1]][0]
    best_cold_mass = coldMass[bestIndices[2]][0]
    best_warm_mass = warmMass[bestIndices[0]][0]

    # Create final SED for all wavelengths
    warm_Sed = np.asarray(ModBB(warmKappaFile,area,warmTemp,best_warm_mass)).transpose()
    cold_Sed  = np.asarray(ModBB(coldKappaFile,area,best_cold_temp,best_cold_mass)).transpose()
    total_Sed = cold_Sed + warm_Sed 
    calc_sed = total_Sed[[index]]
    


    # Print Updates
    print(("Measured SED: {}").format(measured_sed))
    print(("Calculated SED: {}").format(calc_sed))
    print(("Chi Squared Value: {}").format(bestError))
    return total_Sed, warm_Sed, cold_Sed, best_cold_temp, best_cold_mass, best_warm_mass, bestError,chiSquareCube  


def CalculateChi(coldTemp,coldMass, warmMass,coldKappaFile,area,measured_sed):
    # Constants 
    warmTemp = 145 ; warmKappaFile = 'Kappa/kappa_mg2sio4.dat'
    lam, warmKappa = log_extrap(warmKappaFile)
    coldTempv, coldMassv, warmMassv = np.meshgrid(coldTemp, coldMass, warmMass) 

    # Get the error based on input measured sed
    sigma = error(measured_sed)

    # ModBB Parameters: (kappa_arr,area,wav,lam,Temp,Mass)

    coldComponentFinal = np.asarray(ModBB(coldKappaFile,area,coldTempv,coldMassv)).transpose()
    warmComponentFinal = np.asarray(ModBB(warmKappaFile,area,warmTemp,warmMassv)).transpose()

    # For  Cold AMC
    sed_c_final = coldComponentFinal + warmComponentFinal
    index = GetObservationIndex(lam)
    chi_final = (measured_sed - sed_c_final[:,:,:,index])**2 / (sigma)**2
    chi_final_sumd = np.sum(chi_final, axis=3)
    # Generalize so the warm mass could be fixed
    # if type(warmMass) == np.float64: 
    #     chi_map = chi_final_sumd[:][:][0]   
    #     print("Fixed Variable")
    # elif type(coldTemp) == np.int64:
    #     chi_map = chi_final_sumd[:,0,:]
    # elif type(coldMass) == np.float64:        
    #     chi_map = chi_final_sumd[:,:,0]
    # else:
    chi_map = chi_final_sumd  
    return chi_map

###########
# Testing #
###########
if test: 
    # Ensure that this SED code will accurately predict its own measurements.
    FixedWarmMass = 4e-5; FixedWarmTemp = 145
    FixedColdMass = .1  ; FixedColdTemp = 30

    import matplotlib.pyplot as plt

    coldKappaFile = 'Kappa/kappa_amc.dat'
    warmKappaFile = 'Kappa/kappa_mg2sio4.dat'

    lam, k_amc = log_extrap(coldKappaFile) 
    NotRequired, Area = AverageSED('../Final_Files/24/24_SNR_Convolve_Regrid_Prune.fits')

    # Create Parameter Space
    temp = np.arange(2,70,1) # Kelvin
    mass = 10**np.arange(-4,0.1,0.1) 
    wmass = 10**np.arange(-8,-3,0.1)

    # Get Fake Observations
    index = GetObservationIndex(lam)
    warm_Sed = np.asarray(ModBB(warmKappaFile,Area,FixedWarmTemp,FixedWarmMass)).transpose()
    cold_Sed  = np.asarray(ModBB(coldKappaFile,Area,FixedColdTemp,FixedColdMass)).transpose()
    total_Sed = cold_Sed + warm_Sed 
    calc_sed = total_Sed[[index]]

    # Try to Get Back Fixed Input
    sed, wsed,csed,t, m, wm, chi = CalculateBestSed(temp, mass, wmass,'Kappa/kappa_amc.dat',Area, calc_sed)
    
    # Plot 
    plt.figure(1)
    plt.scatter([24,70,100,160], calc_sed)
    plt.plot(lam,sed)
    print(("Temp {} Mass {} WarmMass {} ").format(t,m,wm))
    plt.show()



