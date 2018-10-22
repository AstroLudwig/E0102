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
from timeit import default_timer as timer
import numpy as np
from astropy import units as u
from astropy.io import fits
from astropy import constants as const
from scipy.optimize import curve_fit
import numexpr as ne

test = False

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

def ModBB(kappa_arr,area,wav,lam,Temp,Mass):
    M_sun = const.M_sun.cgs
    # If you want the ModBB solution for multiple wavelengths.
    if type(wav) == list or type(wav) == np.ndarray:
        Ans = []
        for i in range(len(wav)):
            kappa = kappa_arr[np.where(np.isclose(lam,wav[i]))[0][0]] * u.cm * u.cm / u.g
            E = Mass * M_sun / area
            B = B_nu(wav[i], Temp)
            ans = kappa * E * B 
            if ans.unit != "MJy/sr":
                print(("Warning, ModBB units = {}").format(ans.unit))
            Ans.append((ans).value) 

    return np.asarray(Ans  )

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
# This returns that average and the physical area of the remnant. 
def sed_avg(imgfile): 
    img = fits.open(imgfile)[0].data
    look = np.where(np.isfinite(img))
    count = np.shape(look)[1]
    avg = np.nansum(img) / count 
    pix_area = pixarea(imgfile)     
    area = (pix_area * count).to(u.cm*u.cm)
    return avg,area 

#############
# Fit SED   #
#############

def CalculateBestSed(coldTemp, warmTemp, coldMass, warmMass, coldKappa, warmKappa,lam,area,measured_sed):
    startC = timer()


    coldTempv, coldMassv, warmMassv = np.meshgrid(coldTemp, coldMass, warmMass) 

    # Get the error based on input measured sed
    sigma = error(measured_sed)

    # ModBB Parameters: (kappa_arr,area,wav,lam,Temp,Mass)

    coldComponentFinal = np.asarray(ModBB(coldKappa,area,[24,70,100,160],lam,coldTempv,coldMassv)).transpose()
    warmComponentFinal = np.asarray(ModBB(warmKappa,area,[24,70,100,160],lam,warmTemp,warmMassv)).transpose()

    # For  Cold AMC
    sed_c_final = coldComponentFinal + warmComponentFinal
    chi_final = (measured_sed - sed_c_final)**2 / (sigma)**2
    chi_final_sumd = np.sum(chi_final, axis=3)
    best_error = np.min(chi_final_sumd)
    bestIndices = np.where(chi_final_sumd == best_error)

    # bestIndices will be [cold_mass_index, cold_temp_index, warm_mass_index]
    best_cold_temp = coldTemp[bestIndices[1]][0]
    best_cold_mass = coldMass[bestIndices[2]][0]
    if len(warmMass) > 1:
        best_warm_mass = warmMass[bestIndices[0]][0] # change to scalar if using warm mass map
    else:
        best_warm_mass = warmMass

    # Create final SED for all wavelength
    warm_Sed = np.asarray(ModBB(warmKappa,area,lam,lam,warmTemp,best_warm_mass)).transpose()
    cold_Sed  = np.asarray(ModBB(coldKappa,area,lam,lam,best_cold_temp,best_cold_mass)).transpose()
    total_Sed = cold_Sed + warm_Sed 

    # Create Calculated SED 
    warm_sed = np.asarray(ModBB(warmKappa,area,[24,70,100,160],lam,warmTemp,best_warm_mass)).transpose()
    cold_sed  = np.asarray(ModBB(coldKappa,area,[24,70,100,160],lam,best_cold_temp,best_cold_mass)).transpose()
    calc_sed = cold_sed + warm_sed 
    
    # Print Updates
    #print("")
    print(("Measured SED: {}").format(measured_sed))
    print(("Calculated SED: {}").format(calc_sed))
    print(("Chi Value: {}").format(best_error))
    endC = timer()
    print(("Loading CalculatedBestSed Packages took: {} s").format(endC - startC))
    return total_Sed, warm_Sed, cold_Sed, best_cold_temp, best_cold_mass, best_warm_mass, best_error

def CalculateBestSedNoBkgd(coldTemp, warmTemp, coldMass, warmMass, coldKappa, warmKappa,lam,area,measured_sed):
    startC = timer()


    coldTempv, coldMassv, warmMassv = np.meshgrid(coldTemp, coldMass, warmMass) 

    # Get the error based on input measured sed
    sigma = errorNoBkgd(measured_sed)

    # ModBB Parameters: (kappa_arr,area,wav,lam,Temp,Mass)

    coldComponentFinal = np.asarray(ModBB(coldKappa,area,[24,70,100,160],lam,coldTempv,coldMassv)).transpose()
    warmComponentFinal = np.asarray(ModBB(warmKappa,area,[24,70,100,160],lam,warmTemp,warmMassv)).transpose()

    # For  Cold AMC
    sed_c_final = coldComponentFinal + warmComponentFinal
    chi_final = (measured_sed - sed_c_final)**2 / (sigma)**2
    chi_final_sumd = np.sum(chi_final, axis=3)
    best_error = np.min(chi_final_sumd)
    bestIndices = np.where(chi_final_sumd == best_error)

    # bestIndices will be [cold_mass_index, cold_temp_index, warm_mass_index]
    best_cold_temp = coldTemp[bestIndices[1]][0]
    best_cold_mass = coldMass[bestIndices[2]][0]
    if len(warmMass) > 1:
        best_warm_mass = warmMass[bestIndices[0]][0] # change to scalar if using warm mass map
    else:
        best_warm_mass = warmMass
    # Create final SED for all wavelength
    warm_Sed = np.asarray(ModBB(warmKappa,area,lam,lam,warmTemp,best_warm_mass)).transpose()
    cold_Sed  = np.asarray(ModBB(coldKappa,area,lam,lam,best_cold_temp,best_cold_mass)).transpose()
    total_Sed = cold_Sed + warm_Sed 

    # Create Calculated SED 
    warm_sed = np.asarray(ModBB(warmKappa,area,[24,70,100,160],lam,warmTemp,best_warm_mass)).transpose()
    cold_sed  = np.asarray(ModBB(coldKappa,area,[24,70,100,160],lam,best_cold_temp,best_cold_mass)).transpose()
    calc_sed = cold_sed + warm_sed 
    
    # Print Updates
    #print("")
    print(("Measured SED: {}").format(measured_sed))
    print(("Calculated SED: {}").format(calc_sed))
    print(("Chi Value: {}").format(best_error))
    endC = timer()
    print(("Loading CalculatedBestSed Packages took: {} s").format(endC - startC))
    return total_Sed, warm_Sed, cold_Sed, best_cold_temp, best_cold_mass, best_warm_mass, best_error    

def CalculateChi(coldTemp, warmTemp, coldMass, warmMass,coldKappa, warmKappa,lam,area,measured_sed):
    startC = timer()

  #  warmMass = np.asarray([1.9717e-7])
    
    coldTempv, coldMassv, warmMassv = np.meshgrid(coldTemp, coldMass, warmMass) 

    # Get the error based on input measured sed
    sigma = error(measured_sed)

    # ModBB Parameters: (kappa_arr,area,wav,lam,Temp,Mass)

    coldComponentFinal = np.asarray(ModBB(coldKappa,area,[24,70,100,160],lam,coldTempv,coldMassv)).transpose()
    warmComponentFinal = np.asarray(ModBB(warmKappa,area,[24,70,100,160],lam,warmTemp,warmMassv)).transpose()

    # For  Cold AMC
    sed_c_final = coldComponentFinal + warmComponentFinal
    chi_final = (measured_sed - sed_c_final)**2 / (sigma)**2
    chi_final_sumd = np.sum(chi_final, axis=3)

    return chi_final_sumd[:][:][0] 

def CalculateChiNoBkgd(coldTemp, warmTemp, coldMass, warmMass,coldKappa, warmKappa,lam,area,measured_sed):
    startC = timer()

  #  warmMass = np.asarray([1.9717e-7])
    
    coldTempv, coldMassv, warmMassv = np.meshgrid(coldTemp, coldMass, warmMass) 

    # Get the error based on input measured sed
    sigma = errorNoBkgd(measured_sed)

    # ModBB Parameters: (kappa_arr,area,wav,lam,Temp,Mass)

    coldComponentFinal = np.asarray(ModBB(coldKappa,area,[24,70,100,160],lam,coldTempv,coldMassv)).transpose()
    warmComponentFinal = np.asarray(ModBB(warmKappa,area,[24,70,100,160],lam,warmTemp,warmMassv)).transpose()

    # For  Cold AMC
    sed_c_final = coldComponentFinal + warmComponentFinal
    chi_final = (measured_sed - sed_c_final)**2 / (sigma)**2
    chi_final_sumd = np.sum(chi_final, axis=3)

    return chi_final_sumd[:][:][0]     
###########
# Testing #
###########
if test: 
    # Ensure that SED code will accurately predict it's own measurements.
    import matplotlib.pyplot as plt
    # From some pixel
    sed_test = [1.380621250859944, 3.7808188648242655, 6.4110121622656484, 4.1580572596900751]

    lam, k_amc = log_extrap('Kappa/kappa_amc.dat') 
    lam, k_mg2 = log_extrap('Kappa/kappa_mg2sio4.dat')

    area = pixarea('../Final_Files/24/Final_24_SNR_CR_Prnd.fits') 

    temp = np.arange(2,100,1) # Kelvin

    # Fractional Mass of the Sun using Dex
    mass = 10**np.arange(-6,-2,0.1) 
    wmass = 10**np.arange(-8,-3,0.1)

    tempv,massv,wmassv = np.meshgrid(temp,mass,wmass)

    planck = {}
    

    sed, wsed,csed,t, m, wm, chi = CalculateBestSed(temp, 145, mass, wmass,k_amc, k_mg2, lam, area, sed_test,planck)
    plt.figure(1)
    plt.scatter([24,70,100,160], sed_test)
    plt.plot(lam,sed)
    print(("Temp {} Mass {} WarmMass {} ").format(t,m,wm))
    plt.show()
