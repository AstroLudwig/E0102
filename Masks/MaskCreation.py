# -*- Copyright (c) 2018, Bethany Ann Ludwig, All rights reserved. -*-
import sys
sys.path.insert(0, '../Convolve_Regrid')

import numpy as np  
from astropy.io import fits 
import matplotlib.pyplot as plt 
import cv2
from astropy.wcs import WCS
from Convolve import master_convolve
#################################################
## ******************************************  ##
##  ~    Mask Coordinate Creation           ~  ##
## ******************************************  ##
#################################################

"""
Creates a file with the world  coordinates of the mask 
achieved by first finding the pixel coordinates around 
the supernova in x-ray where the edges are most clearly defined. 

Convolves Xray to image in question then gets the mask coordinates.

Sanity Check:
 Convolution has zero padding --------> 
 From https://docs.scipy.org/doc/scipy-0.19.0/reference/generated/scipy.signal.fftconvolve.html
 "Gaussian blur implemented using FFT convolution. 
 Notice the dark borders around the image, due to the zero-padding beyond its boundaries. 
 The convolve2d function allows for other types of image boundaries, but is far slower."

"""
#######################
##     Switches      ##
#######################
# Choose how many pixels to extend by
extend = 20
# Extend the mask isotropically. 
push_all = False
# Extend the mask to the left (away from nebula.)
push_left = True

# If the file doesn't exist, it needs to be convolved.
convolve = True 

# Choose which mask to create
mask24 = True
mask70 = False
mask100 = False
mask160 = False

# Check things are working
plot = True

# Save game
save = True


#######################
##  Loading Data     ##
#######################

# Using Chandra data to define the boundary of the shock
xray_file = '../Original_Files/Xray/e0102_1100-2000eV.fits'
xray_data = fits.open(xray_file)[0].data
xray_hdr = fits.open(xray_file)[0].header

# Kernels
if mask24: 
    kern = '../Convolve_Regrid/Kernels/Kernel_LoRes_Gauss_02.0_to_MIPS_24.fits'
    con_sname = '../Convolve_Regrid/Convolutions/Xray_to_24.fits'
    save_name = '../Masks/Xray/Mask_World_24_ext_'+str(extend)+'.txt'
if mask70:
    kern = '../Convolve_Regrid/Kernels/Kernel_LoRes_Gauss_02.0_to_PACS_70.fits'
    con_sname = '../Convolve_Regrid/Convolutions/Xray_to_70.fits'
    save_name = '../Masks/Xray/Mask_World_70_ext_'+str(extend)+'.txt'
if mask100:
    kern= '../Convolve_Regrid/Kernels/Kernel_LoRes_Gauss_02.0_to_PACS_100.fits'
    con_sname = '../Convolve_Regrid/Convolutions/Xray_to_100.fits'
    save_name = '../Masks/Xray/Mask_World_100_ext_'+str(extend)+'.txt'
if mask160:
    kern = '../Convolve_Regrid/Kernels/Kernel_LoRes_Gauss_02.0_to_PACS_160.fits'
    con_sname = '../Convolve_Regrid/Convolutions/Xray_to_160.fits'
    save_name = '../Masks/Xray/Mask_World_160_ext_'+str(extend)+'.txt'

if convolve:    
    conv = master_convolve(kern, xray_file, con_sname)
else: 
    conv = fits.open(con_sname)[0].data

###############################
##  Find Mask Coordinates    ##
###############################

#Edge Detection               
VisCopy = np.uint8(conv)
edges = cv2.Canny(VisCopy,0,5)   
# Some threshhold for what constitutes as an edge.
points = np.where(edges == 255)
x = points[1]
y = points[0]
# Save coordinates
all_coords = np.asarray(list(zip(x,y)))
    
##################################
##  Reduce Mask Coordinates    ##
#################################

# Create mask arrays. 
maskx_x = []
maskx_y = []

for i in range(np.max(x)): 
    # remove duplicates in X          
    look = np.where(x == i)[0]
    dim = np.shape(look)
    if dim[0] > 0: 
         maskx_y.append(y[look[0]])
         maskx_y.append(y[look[len(look)-1]])
         maskx_x.append(x[look[0]])
         maskx_x.append(x[look[len(look)-1]])
        
masky_x = []
masky_y = []         
for i in range(np.max(y)):  
    # remove duplicates in Y
    look = np.where(y == i)[0]
    dim = np.shape(look)
    if dim[0] > 0: 
         masky_y.append(y[look[0]])
         masky_y.append(y[look[len(look)-1]])
         masky_x.append(x[look[0]])
         masky_x.append(x[look[len(look)-1]])

mask_x = maskx_x + masky_x         
mask_y = maskx_y + masky_y

mask = np.asarray(list(zip(mask_x,mask_y)))

############################
##  Extend Coordinates    ##
############################


mid_x = (np.max(x)+np.min(x))/2; mid_y = (np.max(y)+np.min(y))/2
mask_x = np.asarray(mask_x);     mask_y = np.asarray(mask_y)
push_x = np.copy(mask_x);        push_y = np.copy(mask_y)

if push_all:
    for i in range(len(mask_x)): 
        if mask_x[i] < mid_x:
            push_x[i] = mask_x[i] - extend
        if mask_x[i] > mid_x:
            push_x[i] = mask_x[i] + extend
        if mask_y[i] > mid_y:
            push_y[i] = mask_y[i] + extend
        if mask_y[i] < mid_y:
            push_y[i] = mask_y[i] - extend

    # ----> Top and Bottom 
    fpx = mid_x - 1 - extend
    lpx = mid_x + 1 + extend 
    fpy = np.max(mask_y) + extend
    fppy = np.min(mask_y) - extend
    fyrr = np.repeat(fpy,lpx-fpx)
    fypr = np.repeat(fppy,lpx-fpx)           
    fxrr = np.arange(fpx,lpx)
    push_x = np.append(push_x,fxrr)
    push_y = np.append(push_y,fyrr)
    push_x = np.append(push_x,fxrr)
    push_y = np.append(push_y,fypr)
    # -----> Left and Right
    kpx = np.max(mask_x) + extend 
    kppx = np.min(mask_x) - extend
    kpy = mid_y - 1 - extend
    kppy = mid_y + 1 + extend 
    kxrr = np.repeat(kpx,kppy - kpy)
    kxpr = np.repeat(kppx,kppy - kpy)
    kyrr = np.arange(kpy,kppy)
    push_x = np.append(push_x,kxrr)
    push_y = np.append(push_y,kyrr)
    push_x = np.append(push_x,kxpr)
    push_y = np.append(push_y,kyrr)

if push_left: 
    for i in range(len(mask_x)): 
        if mask_x[i] < mid_x:
            push_x[i] = mask_x[i] - extend
        if mask_y[i] > mid_y:
            push_y[i] = mask_y[i] + extend
        if mask_y[i] < mid_y:
            push_y[i] = mask_y[i] - extend
    # ----> Top and Bottom 
    fpx = mid_x - 1 - extend
    lpx = mid_x + 1 
    fpy = np.max(mask_y) + extend
    fppy = np.min(mask_y) - extend
    fyrr = np.repeat(fpy,lpx-fpx)
    fypr = np.repeat(fppy,lpx-fpx)           
    fxrr = np.arange(fpx,lpx)
    push_x = np.append(push_x,fxrr)
    push_y = np.append(push_y,fyrr)
    push_x = np.append(push_x,fxrr)
    push_y = np.append(push_y,fypr)
    # -----> Left and Right
    kpx = np.max(mask_x)  
    kppx = np.min(mask_x) - extend
    kpy = mid_y - 1 - extend
    kppy = mid_y + 1 + extend 
    kxrr = np.repeat(kpx,kppy - kpy)
    kxpr = np.repeat(kppx,kppy - kpy)
    kyrr = np.arange(kpy,kppy)
    push_x = np.append(push_x,kxrr)
    push_y = np.append(push_y,kyrr)
    push_x = np.append(push_x,kxpr)
    push_y = np.append(push_y,kyrr)

mask_x = push_x
mask_y = push_y

##################################
##  Convert Mask Coordinates    ##
#################################
         
          
#Convert From Pixel to World 
pix_x = np.asarray(mask_x)
pix_y = np.asarray(mask_y)  
pix = list(zip(pix_x,pix_y))


w = WCS(xray_hdr)
world = w.all_pix2world(pix, 0) 


##################################
##    Save Mask Coordinates    ##
#################################

#Save World Coordinates of Mask For Use in Next Program
if save:
    np.savetxt(save_name,(world[:,0],world[:,1]))


#######################
##    Plot Result    ##
#######################
if plot:
    # Load Mask Coordinates that were just saved
    ra,dec = np.loadtxt(save_name)
    # Open 24 micron file as an example
    infra_data = fits.open('../Sky_Remove/Median_Removed/24um_medianRemoved.fits')[0].data
    infra_hdr = fits.open('../Sky_Remove/Median_Removed/24um_medianRemoved.fits')[0].header
    # Convert coordinates
    W = WCS(infra_hdr)
    X,Y = W.all_world2pix(ra,dec,1)
    #Show with mask
    plt.figure(2)
    plt.scatter(X,Y, s= .5) 
    plt.imshow(infra_data,vmin=-1,vmax=1,cmap='gray')
    plt.ylim(100,200)
plt.show()




