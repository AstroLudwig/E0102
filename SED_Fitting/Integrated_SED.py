from timeit import default_timer as timer

from astropy.io import fits

import numpy as np

import matplotlib.pyplot as plt 

import Eqs


##############
## Switches ##
##############

# Fit the entire remnant rather than pixel by pixel.
generate_general = True
plot_general = True

# Pix by Pix
# Creates Pixel Sed/Temp/Mass solutions and stores them in text file.
generate = True  
# Choose whether to evalute the cold dust composition as carbon or silicate.
generate_AMC = True
generate_MG2 = False

# Fit only pixels that fall within a detection limit. 
# If running for something other than 2 sigma will need to change n. 
generate_limit = False
# Try without background removal error.
generate_nobkgderr = True
# Generate chi map, "Warm_Mass_Map" will need to be changed if running
# something other than the full map with cold AMC.  
generate_chiMap = True

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

#########################
#    General Fitting    #
#########################
# Get SNR data, the physical area of each SNR, and the average pixel intensity
data = []; Areas = []; AverageIntensities = np.zeros(4)
for i in range(4):
    AverageIntensity, Area = Eqs.AverageSED(files[i])
    data.append(fits.open(files[i])[0].data)
    Areas.append(Area)
    AverageIntensities[i] = AverageIntensity
# Get the error associated with each intensity
AverageError = Eqs.error(AverageIntensities)

if generate_general:   
    total_sed, warm_sed, cold_sed, temp, cold_mass, warm_mass, chi_squared = Eqs.CalculateBestSed(ColdTemp,ColdMass,WarmMass,kappa[0],Areas[0],AverageIntensities)
    chiMap = Eqs.CalculateChi(ColdTemp, ColdMass, WarmMass,kappa[0],Areas[0],AverageIntensities)

if plot_general:
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
#########################
#    Pixel by Pixel     #
#########################
"""
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
"""
plt.show()