# -*- coding: utf-8 -*-
"""
Created on Tue Mar 27 16:20:44 2018

@author: Betha
"""
homogenized_check = False # Check to see if things are changing
step1 = False # Check that mask is working.
step2 = False # Plot background and subtraction
step3 = True # Save
step4 = False # Save to test file
f_24um = False
f_24um_extend = False
f_70um = False
f_70um_extend = False
f_100um = False
f_160um = False
# Ratios
f70_100 = False
f70_160 = False
f100_160 = True


ext = 0 # If you're using an extended mask, select file extension number here.
vmin = 0; vmax=30
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS


if f_24um:
    # Original, Median Removed File
    f = fits.open('../../Sky_Remove/Median_Removed/24um_medianRemoved.fits')
    # Import mask coordinates
    ra,dec = np.loadtxt('../../Masks/Xray/Mask_World_24_V2.txt')
    um = '24' # for saving file types
    folder = 'im24/' #folder name for saving
    data = f[0].data; hdr = f[0].header; wcs = WCS(hdr)
    x,y = np.round(wcs.all_world2pix(ra,dec,1))
if f_24um_extend:
    print(('Diffusing 24 microns extended by {}').format(ext))
    # Original, Median Removed File
    f = fits.open('../../Sky_Remove/Median_Removed/24um_medianRemoved.fits')
    # Import mask coordinates
    ra,dec = np.loadtxt('../../Masks/Xray/Extend/Mask_World_24_V2_ext_'+str(ext)+'.txt')
    um = '24' # for saving file types
    folder = 'im24_ext_'+str(ext)+'/' #folder name for saving   
    data = f[0].data; hdr = f[0].header; wcs = WCS(hdr)
    x,y = np.round(wcs.all_world2pix(ra,dec,1))
if f_70um:
    # Original, Median Removed File
    f = fits.open('../../Sky_Remove/Median_Removed/70um_medianRemoved.fits')
    # Import mask coordinates
    ra,dec = np.loadtxt('../../Masks/Xray/Mask_World_70.txt')
    um = '70'
    folder = 'im70/' #folder name for saving
if f_70um_extend:
    # Original, Median Removed File
    f = fits.open('../../Sky_Remove/Median_Removed/70um_medianRemoved.fits')
    # Import mask coordinates
    ra,dec = np.loadtxt('../../Masks/Xray/Extend/Mask_World_70_ext_'+str(ext)+'.txt')
    um = '70'
    folder = 'im70_ext_'+str(ext)+'/' #folder name for saving  
    data = f[0].data; hdr = f[0].header; wcs = WCS(hdr)
    x,y = np.round(wcs.all_world2pix(ra,dec,1))
if f70_100:
    # Sky removed image 70 convolved and regridded to 100. SNR NOT removed.
    f70 = fits.open('../../Convolve_Regrid/Conv_Regrids/SkyRemovd70_to_100_CR.fits')
    # Sky removed 100
    f100 = fits.open('../../Sky_Remove/Median_Removed/100um_medianRemoved.fits')
    ra,dec = np.loadtxt('../../Masks/Xray/Mask_World_100.txt')
    um = '100'
    folder = 'ratio_100_to_70/'
    data = f100[0].data/f70[0].data; hdr = f100[0].header; wcs = WCS(hdr)
    x,y = np.round(wcs.all_world2pix(ra,dec,1))
if f70_160:
    # Sky removed image 70 convolved and regridded to 100. SNR NOT removed.
    f70 = fits.open('../../Convolve_Regrid/Conv_Regrids/SkyRemoved70_to_160_CR.fits')
    # Sky removed 100
    f160 = fits.open('../../Sky_Remove/Median_Removed/160um_medianRemoved.fits')
    ra,dec = np.loadtxt('../../Masks/Xray/Mask_World_160.txt')
    um = '160'
    folder = 'ratio_160_to_70/'
    data = f160[0].data/f70[0].data; hdr = f160[0].header; wcs = WCS(hdr)
    x,y = np.round(wcs.all_world2pix(ra,dec,1))    
if f100_160:
    # Sky removed image 100 convolved and regridded to 160. SNR NOT removed.
    f100 = fits.open('../../Convolve_Regrid/Conv_Regrids/SkyRemovd100_to_160_CR.fits')
    # Sky removed 100
    f160 = fits.open('../../Sky_Remove/Median_Removed/160um_medianRemoved.fits')
    ra,dec = np.loadtxt('../../Masks/Xray/Mask_World_160.txt')
    um = '160'
    folder = 'ratio_160_to_100/'
    data = f160[0].data/f100[0].data; hdr = f160[0].header; wcs = WCS(hdr)
    x,y = np.round(wcs.all_world2pix(ra,dec,1))        
if f_100um:
    # Original, Median Removed File
    f = fits.open('../../Sky_Remove/Median_Removed/100um_medianRemoved.fits')
    # Import mask coordinates
    ra,dec = np.loadtxt('../../Masks/Xray/Mask_World_100.txt')
    um = '100'
    folder = 'im100/' #folder name for saving
    data = f[0].data; hdr = f[0].header; wcs = WCS(hdr)
    x,y = np.round(wcs.all_world2pix(ra,dec,1))
if f_160um:
    # Original, Median Removed File
    f = fits.open('../../Sky_Remove/Median_Removed/160um_medianRemoved.fits')
    # Import mask coordinates
    ra,dec = np.loadtxt('../../Masks/Xray/Mask_World_160.txt')
    um = '160'
    folder = 'im160/' #folder name for saving
    data = f[0].data; hdr = f[0].header; wcs = WCS(hdr)
    x,y = np.round(wcs.all_world2pix(ra,dec,1))    
# Checking that mask is correct 
if step1:
    plt.figure()
    plt.imshow(data,vmin=0,vmax=30)
    plt.scatter(x,y,s=0.1,c='r')
    plt.xlim(190,320)
    plt.ylim(190,280)

###############
# Snr Removed #
###############
def removeSnr(data,hdr):
    w = WCS(hdr)
    # For each pixel we need to know it's ra/dec coordinates
    a,b = np.shape(data)
    row = np.arange(0,a); col = np.arange(0,b)
    row,col=np.meshgrid(row,col)
    row=row.flatten(); col=col.flatten()

    all_ra, all_dec = w.all_pix2world(col,row,1)
    # Numbers here are from Karin's paper.
    c1 = SkyCoord(c_ra*u.deg, c_dec*u.deg, distance=d*u.kpc, frame='icrs')
    c2 = SkyCoord(all_ra*u.deg, all_dec*u.deg, distance=d*u.kpc, frame='icrs')

    sep = c1.separation_3d(c2)
    radius = d*u.kpc*arcs/206265

    look =np.where(sep < radius)

    data[row[look],col[look]] = 0
    return data
#############
# Diffusion #
#############

#######################
# General Notes #
#######################
#
# x, y are pixel coordinates of the mask
#

# Find rectangle enclosing all data points in diffusion region.
# (minY, minX) will be the bottom left corner of the region.
minX = int(min(x))
maxX = int(max(x)+1)
minY = int(min(y))
maxY = int(max(y)+1)

# Bundle pixel coordinates spanning the (circular-ish) mask.
xyZip = list(zip(x,y))
xyArray = np.asarray(xyZip)
# - Original_data contains the original pixel intensities.
# - New_data contains the simulated diffusion intensities, and will be
#       a 3d array containing diffusion intensities for each timestep.
# - Index_set contains the indices of the boundary pixels.
# - Timestep is the timestep for which the boundary condition is being applied.
#
# For each index in Index_set, this function sets the corresponding intensity
# value in New_data so that it matches Original_data.
# 
# Intensity_Values is an array of zeros set to be the shape of a rectangle 
# around the Snr.
#
# ApplyBoundaryConditions then fills in the appropriate mask pixels with the 
# corresponding pixel intensity from the original image.
def ApplyBoundaryConditions(original_data, new_data, index_set, timestep) :
    for indexPair in index_set :
        # Converts between index in original image and index in smaller template img. 
        # Diffusion occurs in smaller template image and then gets replaced in the 
        i = int(indexPair[1]) - minY
        j = int(indexPair[0]) - minX
        new_data[i,j,timestep] = original_data[int(indexPair[1]), int(indexPair[0])]
# Apply one step of the diffusion.
def StepForward(data, timestep) :
    xLim, yLim, junk = np.shape(data)
    # Kappa should be small, but is arbitrary as we are using the 
    # homoginzed solution as a final stopping place. 
    alpha = 0.1
    for i in range(1, xLim-1) :
        for j in range(1, yLim-1) :
            # Discretizing heat equation
            # Backward Finite Difference in Time
            # Second Order Central Difference in Space
            data[i,j,timestep] = data[i,j,timestep-1] + alpha * (
                    (data[i+1,j,timestep-1] - 2*data[i,j,timestep-1] +
                    data[i-1,j,timestep-1]) + (data[i,j+1,timestep-1] -
                        2*data[i,j,timestep-1] + data[i,j-1,timestep-1]
                            ))
# Evaluates whether or not a point is inside or on the mask coordinates.            
def InsideSet(index_array, i, j) :
    sameJ = np.where(np.isclose(index_array[:, 1], j))[0]
    sameI = np.where(np.isclose(index_array[:, 0], i))[0]
    if len(sameJ) == 0 or len(sameI) == 0 :
        return False
    possibleIs = index_array[sameJ][:,0]
    possibleJs = index_array[sameI][:,1]
    minI = min(possibleIs)
    minJ = min(possibleJs)
    maxI = max(possibleIs)
    maxJ = max(possibleJs)
    if minI <= i and i <= maxI and minJ <= j and j <= maxJ :
        return True
    return False


saveEveryXFrames = 5
maxNumberOfSteps = 10000
saveCounter = 0
intensity_values = np.zeros([maxY-minY, maxX-minX, maxNumberOfSteps])
ApplyBoundaryConditions(data, intensity_values, xyZip, 0);

# Diffusion Loop
for t in range(1,maxNumberOfSteps+1):
    # Do Diffusion Step
    StepForward(intensity_values, t);
    # Re-plug in the mask intensitity values.
    ApplyBoundaryConditions(data, intensity_values, xyZip, t);

    if homogenized_check:
        previousSum = 0
        sumSoFar = 0
        sumSoFar = np.sum(intensity_values[:,:,t])
        print(sumSoFar-previousSum)
        previousSum = sumSoFar

    if saveCounter != (saveEveryXFrames-1):
        saveCounter = saveCounter + 1
        continue

    saveCounter = 0
    print("Timestep : " + str(t))

    shell = np.copy(data)
    a, b = np.shape(data)
    for i in range(a):
        for j in range(b):
            if InsideSet(xyArray, i, j) :
                shell[j,i] = intensity_values[j-minY, i-minX, t]
    
    
    # Ratio, Special Cases:
    if f100_160:
        shell = np.copy(shell) * fits.open('../../Convolve_Regrid/Conv_Regrids/bkgd100_bootstrappedWith70diffusedbkgd_CRto160.fits')[0].data
        sub = f160[0].data - np.copy(shell)
    if f70_160:
        shell = np.copy(shell)  * fits.open('../../Convolve_Regrid/Conv_Regrids/bkgd70_diffused_6000thstep_20ext_CRto160.fits')[0].data
        sub = f160[0].data - np.copy(shell)        
    if f70_100:  
        shell = np.copy(shell) * fits.open('../../Convolve_Regrid/Conv_Regrids/bkgd70_ext_20_tstep_6000_to_100_CR.fits')[0].data
        sub = f100[0].data - np.copy(shell)
    
    # Most Cases:
    else:
       sub = np.copy(data) - np.copy(shell)

    # Plot each:   
    if step2:
        f, (ax,bx) = plt.subplots(1,2,sharey=False)
        ax.imshow(sub,vmin=vmin,vmax=vmax)
        bx.imshow(shell,vmin=vmin,vmax=vmax)
        ax.set_xlim(220,280)
        ax.set_ylim(210,260)
        bx.set_xlim(190,320)
        bx.set_ylim(190,280)

    # Save
    if step3:
        fits.writeto(folder+'snr'+um+'/'+um+'um_diff_'+str(t)+'_steps_snr.fits',sub,header=hdr,overwrite=True)
        fits.writeto(folder+'bkgd'+um+'/'+um+'um_diff_'+str(t)+'_steps_bkgd.fits',shell,header=hdr,overwrite=True)

    if step4:
        fits.writeto('Sky_Remove/Diff/test/snr/70um_diff_'+str(t)+'_steps_snr.fits',sub,header=hdr,overwrite=True)
        fits.writeto('Sky_Remove/Diff/test/bkgd/70um_diff_'+str(t)+'_steps_bkgd.fits',shell,header=hdr,overwrite=True)

plt.show()
