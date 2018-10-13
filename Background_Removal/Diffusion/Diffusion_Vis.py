# -*- Copyright (c) 2018, Bethany Ann Ludwig, All rights reserved. -*-
"""
NAME:
    Diffusion Visualization
PURPOSE:
    Data visualization to see if diffusion is working and makes sense as a method.
"""
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from matplotlib.widgets import Slider
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy import units as u

############
# Switches # 
############

im_24 = False
im_70 = False

# Open a range of fits files
# Opening all of them at once can crash things depending on your system.
# This functionality gives you more control.
get_data = False

# Set up Plot of Background, Subtraction, and Negatives
# Include some interaction 
slider = False

# Save the slider plots incrementally
save_img = False

# Turn those plots into a video 
save_video = False

# Observe over subtraction
hist = True

#################
# File Handling # 
#################

if im_24:
    pref_sub = 'im24/snr24/24um_diff_'
    pref_bkgd = 'im24/bkgd24/24um_diff_'
    original = '../../Sky_Remove/Median_Removed/24um_medianRemoved.fits'
    vmin = 0; vmax = 3
    xdim = [190,230]; ydim = [130,180]
    folder = 'im24'
    save_name = folder+'/imgs/im24'   
    video_name = '24um_diffusion.avi'
    prefix = 'im24/imgs/im24'

if im_70:
    pref_sub = 'im70/snr70/70um_diff_'
    pref_bkgd = 'im70/bkgd70/70um_diff_'
    original = '../../Sky_Remove/Median_Removed/70um_medianRemoved.fits'
    vmin = 0; vmax = 30
    xdim = [210,275]; ydim = [200,260]
    folder = 'im70'
    save_name = folder+'/imgs/im70'
    video_name = '70um_diffusion.avi'
    prefix = 'im70/imgs/im70ext_20_'

# Define length of video
start = 5
stop =  3000 
steps = 5
# save names
app_sub = '_steps_snr.fits'
app_bkgd = '_steps_bkgd.fits'

# Plot details
colormap = 'gnuplot2'

#####################################
# Quantify Possible OverSubtraction #
#####################################

def count_negs(data):
    x = np.arange(xdim[0],xdim[1])
    y = np.arange(ydim[0],ydim[1])
    x,y = np.meshgrid(x,y)
    newdata = data[x,y]

    a,b = np.where(np.isfinite(newdata))
    c,d = np.shape(np.where(newdata[a,b]<0))
    return d

def negative_img(subtraction_img):
    shell = np.copy(subtraction_img)
    a,b = np.where(np.isnan(shell))
    shell[a,b] = 0
    a,b = np.where(shell>0)
    shell[a,b] = 0

    return shell    

# -*- Set up initial plot -*-
   # Get initial data
if get_data:   
    original_neg = count_negs(fits.open(original)[0].data)
    sub_data=[];bkgd_data=[]; sub_negs = [];bkgd_negs = []; neg_imgs = []
    for i in range(start,stop,5):
        with fits.open(pref_sub+str(i)+app_sub) as hdu_sub:
            sub_data.append(hdu_sub[0].data)
            sub_negs.append(count_negs(hdu_sub[0].data))
            neg_imgs.append(negative_img(hdu_sub[0].data))

        with fits.open(pref_bkgd+str(i)+app_bkgd) as hdu_bkgd:
            bkgd_data.append(hdu_bkgd[0].data)
            bkgd_negs.append(count_negs(hdu_bkgd[0].data))

if slider:
     # Initialize Plot
    f, (ax, bx, cx) = plt.subplots(1, 3, sharey=False)
    plt.subplots_adjust(left=0.25, bottom=0.3, top = 0.9 )
    # ax is the background image
    ax.imshow(bkgd_data[0],cmap=colormap,origin='lower',vmin=vmin,vmax=vmax)
    ax.set_xlim(xdim)
    ax.set_ylim(ydim)
    ax.set_title('Background')

    # bx is the subtracted image
    bx.imshow(sub_data[0],cmap=colormap,origin='lower',vmin=vmin,vmax=vmax)
    bx.set_title('Subtraction')
    bx.set_xlim(xdim)
    bx.set_ylim(ydim)

    # cx is where the negatives are
    cx.imshow(neg_imgs[0],cmap=colormap,origin='lower',vmin=np.min(neg_imgs[0]),vmax=np.max(neg_imgs[0]))
    cx.set_title('Negatives')
    cx.set_xlim(xdim)
    cx.set_ylim(ydim)

    # Parameter Box to show stats
    r = np.arange(start,stop,5)
    ax.plot([1,2,3], label= ("N in Original: {:.3f} ").format(original_neg))
    ax.plot([3,2,1], label= ("N in Background: {:.3f}").format(bkgd_negs[0]))
    ax.plot([2,2,2], label= ("N in Subtraction: {:.3f} ").format(sub_negs[0]))
    ax.plot([1,1,1], label=("Time Step: {:.0f}").format(r[0]))
    ax.legend(bbox_to_anchor=(0.05, 1), loc=1, borderaxespad=0.)
    if save_img:
        plt.savefig(folder+'/imgs/im24_0.png') 

    # -*- Update plot everytime slope/intercept value is changed with slider -*-
    def update(val):
        print("Updating...")
        num = int(slide.val)
        print(num)

        # Clear image
        ax.clear(); bx.clear(); cx.clear()
        # Set image shapes
        ax.set_xlim(xdim); ax.set_ylim(ydim)
        bx.set_xlim(xdim); bx.set_ylim(ydim)
        cx.set_xlim(xdim); cx.set_ylim(ydim)

        # Fill images in
        ax.imshow(bkgd_data[num],cmap=colormap,origin='lower',vmin=vmin,vmax=vmax)
        ax.set_title('Background')

        bx.imshow(sub_data[num],cmap=colormap,origin='lower',vmin=vmin,vmax=vmax)
        bx.set_title('Subtraction')

        cx.imshow(neg_imgs[num],cmap=colormap,origin='lower',vmin=np.min(neg_imgs[num]),vmax=np.max(neg_imgs[num]))
        cx.set_title('Negatives')

        # Update parameters
        ax.plot([1,2,3], label= ("N in Original: {:.3f} ").format(original_neg))
        ax.plot([3,2,1], label= ("N in Background: {:.3f}").format(bkgd_negs[num]))
        ax.plot([2,2,2], label= ("N in Subtraction: {:.3f} ").format(sub_negs[num]))
        ax.plot([1,1,1], label=("Time Step: {:.0f}").format(r[num]))
        ax.legend(bbox_to_anchor=(0.05, 1), loc=1, borderaxespad=0.)

        f.canvas.draw_idle()


    # Create Sliders
    axcolor = 'aliceblue'
    length = len(np.arange(1,steps,5)) - 1
    ax_y = plt.axes([0.25, 0.2, 0.65, 0.03], facecolor=axcolor)
    slide = Slider(ax_y, 'File: ', 0, length, valinit=0, valfmt='%0.0f')
    slide.on_changed(update)

######################################## Negative Plot 

    
if save_img:
         # Initialize Plot
    f, (ax, bx,cx) = plt.subplots(1, 3, sharey=False)
    plt.subplots_adjust(left=0.25, bottom=0.3, top = 0.9 )

    for i in np.arange(0,len(bkgd_data)):
        # ax is the background image
        ax.imshow(bkgd_data[i],cmap=colormap,origin='lower',vmin=vmin,vmax=vmax)
        ax.set_xlim(xdim)
        ax.set_ylim(ydim)
        ax.set_title('Background')

        # bx is the subtracted image
        bx.imshow(sub_data[i],cmap=colormap,origin='lower',vmin=vmin,vmax=vmax)
        bx.set_title('Subtraction')
        bx.set_xlim(xdim)
        bx.set_ylim(ydim)

        # cx is where the negatives are
        if neg_imgs[i].any() == 0:
            vMiN = -1000; vMaX = 1
        else:
            vMiN = np.min(neg_imgs[i]); vMaX = np.max(neg_imgs[i])
        cx.imshow(neg_imgs[i],cmap=colormap,origin='lower',vmin=vMiN,vmax=vMaX)
        cx.set_title('Negatives')
        cx.set_xlim(xdim)
        cx.set_ylim(ydim)

        # Parameter Box
        r = np.arange(start,stop,5)
        ax.plot([1,2,3], label= ("N in Original: {:.3f} ").format(original_neg))
        ax.plot([3,2,1], label= ("N in Background: {:.3f}").format(bkgd_negs[i]))
        ax.plot([2,2,2], label= ("N in Subtraction: {:.3f} ").format(sub_negs[i]))
        ax.plot([1,1,1], label=("Time Step: {:.0f}").format(r[i]))
        ax.legend(bbox_to_anchor=(0.05, 1), loc=1, borderaxespad=0.)
      #  plt.show()

        print(("Saving image number: {}").format(r[i]))
        figure = plt.gcf() # get current figure
        figure.set_size_inches(20, 10)
        plt.savefig(save_name+str(r[i])+'.png')

        ax.clear()
        bx.clear()
        cx.clear()

if save_video:
    import Diffusion_Video

    Diffusion_Video.create_video(prefix,video_name,start,stop)

if hist:
    def Prune(fname):
        # Numerical Values from Sandstrom 2009
        c_dec = -72.03125; c_ra = 16.00875
        d = 61 #kpc
        arcs = 22
        # Get snr data and header
        data = fits.open(fname)[0].data
        hdr = fits.open(fname)[0].header
        w = WCS(hdr)
        # For each pixel we need to know it's ra/dec coordinates
        a,b = np.shape(data)
        row = np.arange(0,a); col = np.arange(0,b)
        row,col=np.meshgrid(row,col)
        row=row.flatten(); col=col.flatten()

        all_ra, all_dec = w.all_pix2world(col,row,1)
        c1 = SkyCoord(c_ra*u.deg, c_dec*u.deg, distance=d*u.kpc, frame='icrs')
        c2 = SkyCoord(all_ra*u.deg, all_dec*u.deg, distance=d*u.kpc, frame='icrs')

        sep = c1.separation_3d(c2)
        radius = d*u.kpc*arcs/206265

        look =np.where(sep > radius)

        data[row[look],col[look]] = 0

        return data

    negs = []; sums = []; pruned = [];
    for i in range(steps_start,steps,5):
        hdu = fits.open(pref_sub+str(i)+app_sub) 
        im = hdu[0].data[np.where(np.isfinite(hdu[0].data))]
        count = len(im[im < 0.])
        negs.append(count)
        sums.append(np.nansum(hdu[0].data))
        pruned.append(np.nansum(Prune(pref_sub+str(i)+app_sub)))
        hdu.close()
    def findPer(arr,percent): # percent should be a decimal
        arr = np.asarray(arr)
        atol = 1
        ind = np.where(np.isclose(arr[0]-percent*arr[0],arr,atol=atol))[0] 
        for i in range(50): # iter is arbitrary, it's assuming a max difference          
            if len(ind) == 0:
                atol += 1
                ind = np.where(np.isclose(arr[0]-percent*arr[0],arr,atol=atol))[0]  
            
        return ind[0]     

    x = np.arange(steps_start,steps,5)
    
    plt.figure(1)
    
    plt.plot(x,sums)
    plt.xlabel('Steps')
    plt.ylabel('Total Sum')
    plt.title(folder)
    
    plt.figure(2)
    end = fits.open(pref_sub+str(stop)+app_sub)[0].data
    endSum = np.nansum(end)
    plt.plot(x,(sums-endSum)/endSum*100)
    plt.xlabel('Steps')
    plt.ylabel('Percent Difference from Most Homogenized Point')
    plt.title(folder)

    
    plt.figure(3)
    plt.plot(x,pruned)
    plt.xlabel('Steps')
    plt.ylabel('Sum of Values within a 22 arcsecond Radius')
    plt.title(folder)
    plt.legend(bbox_to_anchor=(.7, .8), loc=2, borderaxespad=0.)
    
    
    plt.figure(4)
    lastSum = pruned[101]
    diffs = np.abs(pruned-lastSum)/lastSum * 100
    plt.plot(x,diffs)
    plt.xlabel('Steps')
    plt.ylabel('Percent Difference from Most Homogenized Within a 22 arcsecond Radius')
    plt.title(folder)
    
    plt.figure(5)
    origSum = np.sum(Prune(original))
    diffs_ = np.abs(pruned-origSum)/origSum * 100
    plt.plot(x,diffs_)
    plt.xlabel('Steps')
    plt.ylabel('Percent Difference from Original Within a 22 arcsecond Radius')
    plt.title(folder)
    
plt.show()
