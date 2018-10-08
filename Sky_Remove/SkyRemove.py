# -*- Copyright (c) 2018, Bethany Ann Ludwig, All rights reserved. -*-
"""
NAME:
    Sky Remove
PURPOSE:
    To understand noise we need to find a region that exists in all of the images.
    The median in this region should be closest to 0. 
NOTES:
    Doing this by eye seems to work best.
    A previous attempt moved through the entire image algorithmically but 
    always chose the edges which isn't a good representation of the noise 
    because it's likely that the edges aren't scanned over by the instrument as well
    as regions away from the center. I could have moved the algorithm away from the
    edges but being that the 160 image is precariously shaped this proved 
    immensely troublesome. When I picked a region myself it seemed to have reasonable
    histograms in all images. 
"""
############
# Switches #
############
# Displays Region Selected in all 4 Images
plot_region = False
# Displays Data in each Region with Associated Histogram.
plot_all = True
stats = True
save = False
###################
# Import Packages #
###################
from astropy.io import fits 
from astropy.wcs import utils, WCS
import numpy as np 
import matplotlib.patches as patches
import matplotlib.pyplot as plt 
from astropy import units as u
import matplotlib.path as mplPath

################
# Import Files #
################
f24 = '../Saturation_Correction/n76crop_smcsage_24_desaturated.fits'
f70 = '../Original_Files/Infrared/e0102_pacs70_new.fits'
f100 = '../Original_Files/Infrared/e0102_pacs100_new.fits'
f160 = '../Original_Files/Infrared/e0102_pacs160_new.fits'

im24 = fits.open(f24)[0].data
im70 = fits.open(f70)[0].data
im100 = fits.open(f100)[0].data
im160 = fits.open(f160)[0].data

#################
# Define Region #
#################
# xn x yn box 
xs = 76; ys = 79 # Starting point of box in 160 image
xn = 8; yn = 8   # Create box with this length and width.  

##############################
# Get Coordinates for Region #
##############################
# Given ra dec coordinates this returns x, y pixel integers
def get_xy(file,ra,dec):
    w = WCS(fits.open(file)[0].header)
    x,y = w.all_world2pix(ra,dec,1)
    x = np.round(x).astype(int)
    y = np.round(y).astype(int)
    return x,y
# Given x,y coordinates this returns ra, dec coordinates
def get_RaDec(file,x,y):
    w = WCS(fits.open(file)[0].header)
    ra, dec = w.all_pix2world(x,y,1)
    return ra, dec
# Some regions in the image may have rotations, this corrects for that.
# Uses a group of points and creates a polygon around it.
def get_path(x,y):
    bbPath = mplPath.Path(np.array([[x[0], y[0]],
                                    [x[1], y[1]],
                                     [x[3], y[3]],
                                      [x[2], y[2]],
                                       [x[0],y[0]]]))   
    return bbPath
# Uses polygon path to get all of the pixels inside of the image. 
# Does this by turning everything outside of the image into a NAN. 
def get_region(img,path):
    a,b = np.shape(img)
    shell = np.copy(img)
    for i in range(a):
        for j in range(b):
            if path.contains_point((j,i)) == False:
                shell[i,j] = np.nan
    return shell

   
 
# Get pixel coordinates for the 4 corners of the region in 160 um image.
sq = [(xs,ys),(xs+xn,ys),(xs,ys+yn),(xs+xn,ys+yn)]
# Get RA/DEC coordinates 
sq_ra, sq_dec = get_RaDec(f160,[x[0] for x in sq],[y[1] for y in sq])
# Convert to pixel Coordinates in the other Images and get the polygonal 
# path around region.
sx24,sy24 = get_xy(f24,sq_ra,sq_dec)
sx70,sy70 = get_xy(f70,sq_ra,sq_dec)
sx100,sy100 = get_xy(f100,sq_ra,sq_dec)
sx160,sy160 = get_xy(f160,sq_ra,sq_dec)

path24 = get_path(sx24,sy24)
path70 = get_path(sx70,sy70)
path100 = get_path(sx100,sy100)
path160 = get_path(sx160,sy160)

# Get all sampled data within the polygonal path
s24= get_region(im24,path24)
s70= get_region(im70,path70)
s100= get_region(im100,path100)
s160= get_region(im160,path160)

samples = [s24,s70,s100,s160]

    
####################
# View All Regions #
####################

if plot_region: 

# Displays Region Selected in all 4 Images
# With Vmin/Vmax set to data within each region. 
# Environment is shown for context.

    f, (ax,bx,cx,dx) = plt.subplots(1,4,sharey=False)
    files = [im24,im70,im100,im160]
    vmin = 0; vmax = 10

    ax.set_title("24 Micron")
    ax.imshow(files[0],vmin=np.nanmin(samples[0]),vmax=np.nanmax(samples[0]))
    ax.scatter(sx24,sy24)
    patch24 = patches.PathPatch(path24, facecolor='none', lw=2)
    ax.add_patch(patch24)             
    ax.set_xlim(120,190); ax.set_ylim(90,140)

    bx.set_title("70 Micron")
    bx.imshow(files[1],vmin=np.nanmin(samples[1]),vmax=np.nanmax(samples[1]))
    bx.scatter(sx70,sy70)
    patch70 = patches.PathPatch(path70, facecolor='none', lw=2)
    bx.add_patch(patch70)
    bx.set_xlim(100,240); bx.set_ylim(110,210)

    cx.set_title("100 Micron")
    cx.imshow(files[2],vmin=np.nanmin(samples[2]),vmax=np.nanmax(samples[2]))
    cx.scatter(sx100,sy100)
    patch100 = patches.PathPatch(path100, facecolor='none', lw=2)
    cx.add_patch(patch100)
    cx.set_xlim(80,180); cx.set_ylim(90,160)

    dx.set_title("160 Micron")
    dx.imshow(files[3],vmin=np.nanmin(samples[3]),vmax=np.nanmax(samples[3]))
    dx.scatter(sx160,sy160)
    patch160 = patches.PathPatch(path160, facecolor='none', lw=2)
    dx.add_patch(patch160)
    dx.set_xlim(60,110); dx.set_ylim(60,100)
    
    
    
    
if plot_all:
    g, ((ax,ex),(bx,fx),(cx,gx),(dx,hx)) = plt.subplots(4,2,sharey=False)
    
    r = samples

    ax.imshow(r[0],vmin=np.nanmin(r[0]),vmax=np.nanmax(r[0]))
    bx.imshow(r[1],vmin=np.nanmin(r[1]),vmax=np.nanmax(r[1]))
    cx.imshow(r[2],vmin=np.nanmin(r[2]),vmax=np.nanmax(r[2]))
    dx.imshow(r[3],vmin=np.nanmin(r[3]),vmax=np.nanmax(r[3]))
    
    
    ax.set_xlim(145,157)
    ax.set_ylim(112,123)
    bx.set_xlim(134,154)
    bx.set_ylim(145,165)
    cx.set_xlim(113,129)
    cx.set_ylim(122,140)
    dx.set_xlim(75,85)
    dx.set_ylim(79,88)
        
    ax.set_title("24 Micron")
    bx.set_title("70 Micron")
    cx.set_title("100 Micron")
    dx.set_title("160 Micron")
    
    ex.hist(r[0][~np.isnan(r[0])])
    fx.hist(r[1][~np.isnan(r[1])])
    gx.hist(r[2][~np.isnan(r[2])])
    hx.hist(r[3][~np.isnan(r[3])])

    
###############
# Sky Removal #
###############
if save:
    # Aquire Statistics of Regions
    medians = []; stds = []

    for samples in samples:
        medians.append(np.nanmedian(samples))
        stds.append(np.nanstd(samples))
        
    print(" ")    
    print("Sky Removals: ")
    print (("24 Micron: Median {} Std {}").format(medians[0],stds[0]))
    print (("70 Micron: Median {} Std {}").format(medians[1],stds[1]))
    print (("100 Micron:Median {} Std {}").format(medians[2],stds[2]))
    print (("160 Micron:Median {} Std {}").format(medians[3],stds[3]))

    pref = 'Median_Removed/'
    strs = ['24','70','100','160']
    fnames = [f24,f70,f100,f160]
    files = [im24,im70,im100,im160]
    med_suffx = 'um_medianRemoved.fits'
    
    # Save Images with Median Subtracted
    for i in range(4):
        hdr = fits.open(fnames[i])[0].header
        fits.writeto(pref+strs[i]+med_suffx, files[i]-medians[i],hdr,overwrite=True)
    
    # Save Error       
    stds = np.asarray(stds)    
    np.savetxt('Sigma.txt',stds.reshape(1,4),delimiter=" ")
    print("Sky Removed")
    print("All Files Saved.")

plt.show()
