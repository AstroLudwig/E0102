from timeit import default_timer as timer

from astropy.io import fits

import numpy as np

import matplotlib.pyplot as plt 

import Eqs


##############
## Switches ##
##############
# Pix by Pix
# Creates Pixel Sed/Temp/Mass solutions and stores them in text file.
generate = True  
# Choose whether to evalute the cold dust composition as carbon or silicate.
generate_AMC = True
generate_MG2 = False
# Fit the entire remnant rather than pixel by pixel.
generate_general = False
# Fit only pixels that fall within a detection limit. 
# If running for something other than 2 sigma will need to change n. 
generate_limit = False
# Try without background removal error.
generate_nobkgderr = True
# Generate chi map, "Warm_Mass_Map" will need to be changed if running
# something other than the full map with cold AMC.  
generate_chiMap = True

###########
## Files ##
###########
# Insert files to get seds.
files = ['../Final_Files/24/Final_24_SNR_CR_Prnd.fits','../Final_Files/70/Final_70_SNR_CR_Prnd.fits','../Final_Files/100/Final_100_SNR_CR_Prnd.fits','../Final_Files/160/160um_70modeledFS_prnd_snr.fits']

################
# Kappa Ranges #
################
# Extrapolate for entire kappa length.
lam, k_amc = Eqs.log_extrap('../Kappa/kappa_amc.dat')
lam, k_mg2 = Eqs.log_extrap('../Kappa/kappa_mg2sio4.dat')

####################
# Parameter Ranges #
####################
# Create a grid to do a fit
# Temperature
T = np.arange(2,70,1) # Kelvin
# Fractional Mass of the Sun using Dex
# Cold Mass
M = 10**np.arange(-4,0.1,.1)
# Warm Mass
wM = 10**np.arange(-8,-3,.1)
T_,M_,wM_ = np.meshgrid(T,M,wM)


#########################
#    General Fitting    #
#########################


data = []; sums = []; areas = []; sed_means = []
for i in range(4):
    sm, ar = Eqs.sed_avg(files[i])
    data.append(fits.open(files[i])[0].data)
    sums.append(np.nansum(fits.open(files[i])[0].data))
    areas.append(ar.value)
    sed_means.append(sm)

sigma = Eqs.error(sed_means)

if generate_general:   
    j,genArea = Eqs.sed_avg(files[0])

    cold_kappa = k_amc

    sed1, wsed1, csed1, t1, m1, wm1, chi1 = Eqs.CalculateBestSed(T,145,M,wM,cold_kappa,k_mg2,lam,genArea,sed_means)
    print(t1,m1,wm1)
    np.savetxt("Sols/General_TotalSED.txt",sed1); np.savetxt("Sols/General_WarmSED.txt",wsed1); np.savetxt("Sols/General_ColdSED.txt",csed1)
    np.savetxt("Sols/General_TempMassWarmMassChi.txt",[t1,m1,wm1,chi1])

    cold_kappa2 = k_mg2

    sed2, wsed2, csed2, t2, m2, wm2, chi2 = Eqs.CalculateBestSed(T,145,M,wM,cold_kappa2,k_mg2,lam,genArea,sed_means)

    np.savetxt("Sols/MG2_General_TotalSED.txt",sed2); np.savetxt("Sols/MG2_General_WarmSED.txt",wsed2); np.savetxt("Sols/MG2_General_ColdSED.txt",csed2)
    np.savetxt("Sols/MG2_General_TempMassWarmMassChi.txt",[t2,m2,wm2,chi2])

    if generate_nobkgderr:

        sed, wsed, csed, t, m, wm, chi = Eqs.CalculateBestSedNoBkgd(T,145,M,wM,k_amc,k_mg2,lam,genArea,sed_means)
        #(coldTemp, warmTemp, coldMass, warmMass,coldKappa, warmKappa,lam,area,measured_sed)
        chiMappa = Eqs.CalculateChiNoBkgd(T, 145, M, wm,k_amc, k_mg2,lam,genArea,sed_means)
        np.savetxt("Sols/BkgdErrRemovd_General_TotalSED.txt",sed); np.savetxt("Sols/BkgdErrRemovd_General_WarmSED.txt",wsed); np.savetxt("Sols/BkgdErrRemovd_General_ColdSED.txt",csed)
        np.savetxt("Sols/BkgdErrRemovd_General_TempMassWarmMassChi.txt",[t,m,wm,chi])

        a,b = np.where(np.isclose(chiMappa,chi,atol=0.1))
        print("Chi Squared Width is "+str(len(a)))
#########################
#    Pixel by Pixel     #
#########################

pix_int = []; # 4 Value Intensity Array
pix_area = Eqs.pixarea(files[0]) # Any file could work all the areas and dimensions are the same.

# Create empty SED Maps
row, col = np.shape(data[0])
SEDmap = np.zeros(row*col*len(lam)).reshape(row,col,len(lam)) 
TotalSEDmap = np.copy(SEDmap); WarmSEDmap = np.copy(SEDmap); ColdSEDmap = np.copy(SEDmap);

# Create empty Temp/Mass/Chi Maps
Tmap = np.copy(data[0]); Tmap[:] = 0
TotalMassMap = np.copy(Tmap); ColdMassMap = np.copy(Tmap); 
WarmMassMap = np.copy(Tmap); ChiMap = np.copy(Tmap); 

# Parameter Chi Map
ParaChiMap = np.zeros(row*col*len(T)*len(M)).reshape(row,col,len(T),len(M))
WidthMap01 = np.copy(Tmap); WidthMap10 = np.copy(Tmap)
# Initialize Select Temp Maps
SelectTmap = np.copy(data[0]); SelectTmap[:] = 0


# Which cold kappa to model
if generate_AMC:
    cold_kappa = k_amc
if generate_MG2:
    cold_kappa = k_mg2

print("Starting SED Fitting...")
startSED = timer() # To pacify the insanity occured while waiting for things.
                   # Depending on how big the parameter range is this can really take a while.          
# Detection Limited Arrays
if generate_limit:

    for i in range(np.shape(data[0])[0]):
        for j in range(np.shape(data[0])[1]):
            arr_int = [data[0][i,j],data[1][i,j],data[2][i,j],data[3][i,j]]  
            # If the mean of the intensities in each image is less than twice the mean error than remove it. 
            if np.nanmean(arr_int) < 2 * np.nanmean(Eqs.error(arr_int)):
                data[0][i,j] = np.nan; data[1][i,j] = np.nan; data[2][i,j] = np.nan; data[3][i,j] = np.nan
    
    a,b = np.where(np.isfinite(data[0]))
    print("Number of pixels has changed to: "+str(len(a)))
    fits.writeto("Sols/DetectionLimited_160umImg.fits",data[3],overwrite=True)
if generate:
    count = 0 
    for i in range(np.shape(data[0])[0]):
        for j in range(np.shape(data[0])[1]):
            if np.isfinite(data[0][i,j]) and np.isfinite(data[1][i,j]) and np.isfinite(data[2][i,j]) and np.isfinite(data[3][i,j]):
                # If pixel in each image contains an intensinty then create array of 4 intensity values.
                arr_sed = [data[0][i,j],data[1][i,j],data[2][i,j],data[3][i,j]]            
                pix_int.append(arr_sed)
                # Start Timer
                count += 1 
                print(("Working on pixel: {} of 189").format(count))

                # Calculated SED for each pixel and save in text file. Can turn this off after running once.
                if generate:
                    # Parameters:(coldTemp, warmTemp, coldMass, warmMass, coldKappa, warmKappa,lam,area,measured_sed)
                    if generate_nobkgderr:
                        print("Removing Background Error")
                        sed, wsed, csed, t, m, wm, chi = Eqs.CalculateBestSedNoBkgd(T,145,M,wM,cold_kappa,k_mg2,lam,pix_area,arr_sed)
                    else:
                        sed, wsed, csed, t, m, wm, chi = Eqs.CalculateBestSed(T,145,M,wM,cold_kappa,k_mg2,lam,pix_area,arr_sed)
                    
                    if generate_chiMap:    
                        # May need to select a different warm mass map depending on what is running.
                        WarmMassMap = np.loadtxt("Sols/WarmMassMap.txt")
                        ParaChi = Eqs.CalculateChi(T, 145, M, WarmMassMap[i,j], cold_kappa, k_mg2,lam,pix_area,arr_sed)
                        print(("Chi VALS.... chi: {}, ParaChi: {}").format(chi,np.min(ParaChi)))
                        print(("Coords: {},{}").format(j,i))
                        ParaChiMap[i,j] = ParaChi

                        row,col = np.where(np.isclose(ParaChi,chi,atol=0.1))
                        row_,col_ = np.where(np.isclose(ParaChi,chi,atol=1))
                        print("Chi Squared Width within 0.1 is "+str(len(row)))
                        print("Chi Squared Width within 1.0 is "+str(len(row_)))

                        WidthMap01[i,j] = len(row)
                        WidthMap10[i,j] = len(row_)

                    # Save Maps. 
                    Tmap[i,j] = t # Temp
                    ChiMap[i,j] = chi # Chi
                    TotalMassMap[i,j] = m + wm; WarmMassMap[i,j] = wm; ColdMassMap[i,j] = m # Mass
                    
                    TotalSEDmap[i,j] = sed; WarmSEDmap[i,j] = wsed; ColdSEDmap[i,j] = csed # SEDs

    endSED = timer() 
    # For sanity
    print(("Generating SEDs took: {} s").format(endSED - startSED))

# Saving and Loading things from files to speed up plotting. Doing Cold Emissivity Types seperate to speed up also.

if generate_AMC:
# Maps
    np.savetxt("Sols/TemperatureMap.txt",Tmap)
    np.savetxt("Sols/TotalMassMap.txt",TotalMassMap)
    np.savetxt("Sols/WarmMassMap.txt",WarmMassMap)
    np.savetxt("Sols/ColdMassMap.txt",ColdMassMap)
    np.savetxt("Sols/ChiSqrdMap.txt",ChiMap)
    # All 3D data must be saved in npy format.
    np.save("Sols/Total_SED_Map.npy", TotalSEDmap)
    np.save("Sols/Warm_SED_Map.npy", WarmSEDmap)
    np.save("Sols/Cold_SED_Map.npy", ColdSEDmap)
    np.save("Sols/Parameter_Chi_Map.npy",ParaChiMap)
        

if generate_MG2:
    np.savetxt("Sols/MG2_TemperatureMap.txt",Tmap)
    np.savetxt("Sols/MG2_TotalMassMap.txt",TotalMassMap)
    np.savetxt("Sols/MG2_WarmMassMap.txt",WarmMassMap)
    np.savetxt("Sols/MG2_ColdMassMap.txt",ColdMassMap)    
    np.savetxt("Sols/MG2_ChiSqrdMap.txt",ChiMap)
    np.save("Sols/Total_SED_Map_MG2.npy", TotalSEDmap)
    np.save("Sols/Warm_SED_Map_MG2.npy", WarmSEDmap)
    np.save("Sols/Cold_SED_Map_MG2.npy", ColdSEDmap)

if generate_limit:
    n = 2
    # Maps
    np.savetxt("Sols/"+str(n)+"SigmaDetectionLimited_TemperatureMap.txt",Tmap)
    np.savetxt("Sols/"+str(n)+"SigmaDetectionLimited_TotalMassMap.txt",TotalMassMap)
    np.savetxt("Sols/"+str(n)+"SigmaDetectionLimited_WarmMassMap.txt",WarmMassMap)
    np.savetxt("Sols/"+str(n)+"SigmaDetectionLimited_ColdMassMap.txt",ColdMassMap)
    np.savetxt("Sols/"+str(n)+"SigmaDetectionLimited_ChiSqrdMap.txt",ChiMap)
    # All 3D data must be saved in npy format.
    np.save("Sols/"+str(n)+"SigmaDetectionLimited_Total_SED_Map.npy", TotalSEDmap)
    np.save("Sols/"+str(n)+"SigmaDetectionLimited_Warm_SED_Map.npy", WarmSEDmap)
    np.save("Sols/"+str(n)+"SigmaDetectionLimited_Cold_SED_Map.npy", ColdSEDmap)
    np.save("Sols/"+str(n)+"SigmaDetectionLimited_Parameter_Chi_Map.npy",ParaChiMap)

if generate_nobkgderr:
    # Maps
    np.savetxt("Sols/BkgdErrRemovd_TemperatureMap.txt",Tmap)
    np.savetxt("Sols/BkgdErrRemovd_TotalMassMap.txt",TotalMassMap)
    np.savetxt("Sols/BkgdErrRemovd_WarmMassMap.txt",WarmMassMap)
    np.savetxt("Sols/BkgdErrRemovd_ColdMassMap.txt",ColdMassMap)
    np.savetxt("Sols/BkgdErrRemovd_ChiSqrdMap.txt",ChiMap)
    np.savetxt("Sols/BkgdErrRemovd_WidthMapPoint1.txt",WidthMap01)
    np.savetxt("Sols/BkgdErrRemovd_WidthMap1.txt",WidthMap10)
    # All 3D data must be saved in npy format.
    np.save("Sols/BkgdErrRemovd_Total_SED_Map.npy", TotalSEDmap)
    np.save("Sols/BkgdErrRemovd_Warm_SED_Map.npy", WarmSEDmap)
    np.save("Sols/BkgdErrRemovd_Cold_SED_Map.npy", ColdSEDmap)
    np.save("Sols/BkgdErrRemovd_Parameter_Chi_Map.npy",ParaChiMap)    

