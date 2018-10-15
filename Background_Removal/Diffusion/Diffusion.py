# -*- Copyright (c) 2018, Bethany Ann Ludwig, All rights reserved. -*-
"""
NAME:
    Diffusion
PURPOSE:
    To Model and Remove the Background of 24 and 70 micron images. 
NOTE: 
    If starting from scratch you will need to create the appropriate folders for this to run.        
"""
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

############
# Switches # 
############
homogenized_check = False # Check to see if things are changing
check_mask = False # Check that mask is working.
check_result = False # Plot background and subtraction
save = True # Save
f_24um = False
f_70um = True
#################
# File Handling #
#################
vmin = 0; vmax=30

if f_24um:
    print('Diffusing 24 microns')
    # Original, Median Removed File
    f = fits.open('../../Sky_Remove/Median_Removed/24um_medianRemoved.fits')
    # Import mask coordinates
    ra,dec = np.loadtxt('../../Masks/Xray/Mask_World_24_ext_20.txt')
    um = '24' # for saving file types
    folder = 'im24/' #folder name for saving   
    data = f[0].data; hdr = f[0].header; wcs = WCS(hdr)
    x,y = np.round(wcs.all_world2pix(ra,dec,1))
if f_70um:
    # Original, Median Removed File
    f = fits.open('../../Sky_Remove/Median_Removed/70um_medianRemoved.fits')
    # Import mask coordinates
    ra,dec = np.loadtxt('../../Masks/Xray/Mask_World_70_ext_20.txt')
    um = '70'
    folder = 'im70/' #folder name for saving  
    data = f[0].data; hdr = f[0].header; wcs = WCS(hdr)
    x,y = np.round(wcs.all_world2pix(ra,dec,1))

# Plot image and mask to make sure files are correct.
if check_mask:
    plt.figure()
    plt.imshow(data,vmin=vmin,vmax=vmax)
    plt.scatter(x,y,s=0.1,c='r')
    plt.xlim(190,320)
    plt.ylim(190,280)

#############
# Diffusion #
#############
#
# x, y are pixel coordinates of the mask
#

# Select region to diffuse over rather than the entire image:

# Define rectangle enclosing all data points in diffusion region.
# (minY, minX) will be the bottom left corner of the region.

# Data within this rectangle but outside the mask is not included
# when subtracting the background. 
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
    # Kappa should be small to avoid converging too quicly, 
    # but is arbitrary as we are using the 
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
# Create empty data cube for just snr region 
# (data outside of the rectangle is extraneous and does not need to be considered)
intensity_values = np.zeros([maxY-minY, maxX-minX, maxNumberOfSteps])
# Fill in the mask coordinates with their original intensity values. 
ApplyBoundaryConditions(data, intensity_values, xyZip, 0);

# Diffusion Loop
for t in range(1,maxNumberOfSteps+1):
    # Do Diffusion Step
    StepForward(intensity_values, t);
    # Re-plug in the mask intensitity values.
    ApplyBoundaryConditions(data, intensity_values, xyZip, t);

    if homogenized_check:
        # Check to see if things have stopped changing.
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
    
    # Remove the background (shell) from the original image.
    sub = np.copy(data) - np.copy(shell)

    # Plot each:   
    if check_result:
        f, (ax,bx) = plt.subplots(1,2,sharey=False)
        ax.imshow(sub,vmin=vmin,vmax=vmax)
        bx.imshow(shell,vmin=vmin,vmax=vmax)
        ax.set_xlim(220,280)
        ax.set_ylim(210,260)
        bx.set_xlim(190,320)
        bx.set_ylim(190,280)

    # Save
    if save:
        fits.writeto(folder+'snr'+um+'/'+um+'um_diff_'+str(t)+'_steps_snr.fits',sub,header=hdr,overwrite=True)
        fits.writeto(folder+'bkgd'+um+'/'+um+'um_diff_'+str(t)+'_steps_bkgd.fits',shell,header=hdr,overwrite=True)

plt.show()
