# -*- coding: utf-8 -*-
"""
Summary Thoughts:
To understand noise we need to find a region that exists in all of the images.
The median in this region should be closest to 0. 
Doing this by eye seems to work best.
A previous attempt moved through the entire image algorithmically but 
always chose the edges which isn't a good representation of the noise 
because it's likely that the edges aren't scanned over by the instrument as well
as regions away from the center. I could have moved the algorithm away from the
edges but being that the 160 image is precariously shaped this proved 
immensely troublesome. When I picked a region myself it seemed to have reasonable
histograms in all images.

Concerns:
Q: Not sure the area is correct. I measured 43.98 pc2
DS9 says 657.9 pc2 but they dont know the distance to the smc
so I dont know how trustworthy it is. Check that this area makes sense.

Methodology: 
I have defined one pixel coordinate in the 160 image at (76,79) and 
extended out 8 pixels. This gives me a region of 64 pixels with which to sample.
 
"""

################
# Switches #
################
# For testing, turn off parts of code

plot = False
plot_corners = False
plot_all = False
stats = True
save = True
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

#f24 = '../Original_Files/mosaic_24_n76.fits'
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
def create_box(img,xstart,ystart):
    x = np.arange(xstart,xstart+xn)
    y = np.arange(ystart,ystart+yn)
    col,row = np.meshgrid(x,y)
    data = img[row,col]
    return data

# provides ra/dec stats about box, only works for 160 as of right now
# to extend it to other images, would need to do xs, xn method differently 
# couple ways to do that, I'm not sure it's necessary right now. 
def equa_box(fname): 
    w = WCS(fits.open(fname)[0].header)
    # Find Center
    center = utils.pixel_to_skycoord(xs + (xn/2),ys + (ys/2),wcs = w)
    # Find corners,  then Area 
    bL = utils.pixel_to_skycoord(xs,ys,wcs = w)
    tL = utils.pixel_to_skycoord(xs,ys+yn,wcs = w)
    width = tL.separation(bL).radian * 60e3 * u.pc 
    bR = utils.pixel_to_skycoord(xs+xn,ys,wcs = w)
    tR = utils.pixel_to_skycoord(xs+xn,ys+yn,wcs=w)
    length = bR.separation(tR).radian * 60e3 * u.pc
    area = width * length
    return center.to_string('hmsdms'), area
    
dbox = create_box(im160,xs,ys)
# Stats
print("Stats for region considered:")
center, area = equa_box(f160)
print(("Center Coordinate: {}").format(center))
print(("Area Subtended: {}").format(area))
print(("Min {} Max {} Mean {} Std {}").format(np.min(dbox),np.max(dbox),np.mean(dbox),np.std(dbox)))


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

   
 
# Get 4 corners
sq = [(xs,ys),(xs+xn,ys),(xs,ys+yn),(xs+xn,ys+yn)]
# Get coordinates for corners, draw path through corners
sq_ra, sq_dec = get_RaDec(f160,[x[0] for x in sq],[y[1] for y in sq])
sx24,sy24 = get_xy(f24,sq_ra,sq_dec)
path24 = get_path(sx24,sy24)
sx70,sy70 = get_xy(f70,sq_ra,sq_dec)
path70 = get_path(sx70,sy70)
sx100,sy100 = get_xy(f100,sq_ra,sq_dec)
path100 = get_path(sx100,sy100)
sx160,sy160 = get_xy(f160,sq_ra,sq_dec)
path160 = get_path(sx160,sy160)
# Get sampled data within the polygon
s24= get_region(im24,path24)
s70= get_region(im70,path70)
s100= get_region(im100,path100)
s160= get_region(im160,path160)
samples = [s24,s70,s100,s160]

    
####################
# View All Regions #
####################

if plot_corners: 
    f4, (ax,bx,cx,dx) = plt.subplots(1,4,sharey=False)
    files = [im24,im70,im100,im160]
    vmin = 0; vmax = 10
    ax.imshow(files[0],vmin=vmin,vmax=vmax)
    ax.scatter(sx24,sy24)
    patch24 = patches.PathPatch(path24, facecolor='none', lw=2)
    ax.add_patch(patch24)             
    bx.imshow(files[1],vmin=vmin,vmax=vmax)
    bx.scatter(sx70,sy70)
    patch70 = patches.PathPatch(path70, facecolor='none', lw=2)
    bx.add_patch(patch70)
    cx.imshow(files[2],vmin=vmin,vmax=vmax)
    cx.scatter(sx100,sy100)
    patch100 = patches.PathPatch(path100, facecolor='none', lw=2)
    cx.add_patch(patch100)
    dx.imshow(files[3],vmin=vmin,vmax=vmax)
    dx.scatter(sx160,sy160)
    patch160 = patches.PathPatch(path160, facecolor='none', lw=2)
    dx.add_patch(patch160)
    ax.set_title("24 Micron")
    bx.set_title("70 Micron")
    cx.set_title("100 Micron")
    dx.set_title("160 Micron")
    
if plot_all:
    f3, ((ax,ex),(bx,fx),(cx,gx),(dx,hx)) = plt.subplots(4,2,sharey=False)
    # I think corners works, not project. Not ready to throw it out yet though
    
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
    
    ex.hist(r[0][~np.isnan(r[0])])#,bins=bins24)
    fx.hist(r[1][~np.isnan(r[1])])#,bins=bins70)
    gx.hist(r[2][~np.isnan(r[2])])#,bins=bins100)
    hx.hist(r[3][~np.isnan(r[3])])#,bins=bins160)


############
# View Box #
############
if plot:
    f, ((ax,cx)) = plt.subplots(2,1,sharey=False)
    # Vmin Vmax kept constant
    ax.imshow(dbox,vmin=np.min(dbox),vmax = np.max(dbox))
    ax.set_title('Region')
    bins160 = np.arange(0,3,1) 
    cx.hist(dbox,bins=bins160)
    cx.set_title('Histogram')
    
    #View Region
    fig = plt.figure()
    zx = fig.add_subplot(1, 1, 1, projection=WCS(fits.open(f160)[0].header))
    #zx.imshow(im160,vmin = np.min(im160),vmax = np.max(im160))
    zx.imshow(im160,vmin = np.min(dbox),vmax = np.max(dbox))
    zx.add_patch(patches.Rectangle((xs,ys),xn,yn,fill=False,edgecolor="red"))
    zx.set_xlim(xs-15,xs+15)
    zx.set_ylim(ys-15,ys+15)
    zx.set_title("Region in Image")
    
###############
# Sky Removal #
###############
if stats:
    means = []; medians = []; stds = []

    for samples in samples:
        means.append(np.nanmean(samples))
        medians.append(np.nanmedian(samples))
        stds.append(np.nanstd(samples))
        
    print(" ")    
    print("Sky Removals: ")
    print (("24 Micron: Mean {} Median {} Std {}").format(means[0],medians[0],stds[0]))
    print (("70 Micron: Mean {} Median {} Std {}").format(means[1],medians[1],stds[1]))
    print (("100 Micron: Mean {} Median {} Std {}").format(means[2],medians[2],stds[2]))
    print (("160 Micron: Mean {} Median {} Std {}").format(means[3],medians[3],stds[3]))


if save:
    pref = 'Median_Removed/'
    pref_ = 'Mean_Removed/'
    strs = ['24','70','100','160']
    fnames = [f24,f70,f100,f160]
    files = [im24,im70,im100,im160]
    med_suffx = 'um_medianRemoved.fits'
    mea_suffx = 'um_meanRemoved.fits'
    for i in range(4):
        hdr = fits.open(fnames[i])[0].header
        fits.writeto(pref+strs[i]+med_suffx, files[i]-medians[i],hdr,overwrite=True)
        fits.writeto(pref_+strs[i]+mea_suffx, files[i]-means[i],hdr,overwrite=True)
    stds = np.asarray(stds)    
    np.savetxt('Sigma.txt',stds.reshape(1,4),delimiter=" ")
    print("Sky Removed")
    print("All Files Saved.")

plt.show()
