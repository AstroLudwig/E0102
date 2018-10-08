im_70 = False
im_70_ext_20 = False
im_24 = False
im_24_ext_5 = False
im_24_ext_10 = False
im_24_ext_20 = True
im_24_ext_40 = False
im_24_ext_60 = False
im_100_div_70 = False
im_160_div_70 = False
im_160_div_100 = False
im_100 = False
im_160 = False

get_data = False
save_img = False
slider = False
save_video = False
hist = True

look = 5
if im_24:
    pref_sub = 'im24/snr24/24um_diff_'
    pref_bkgd = 'im24/bkgd24/24um_diff_'
    original = '../../Sky_Remove/Median_Removed/24um_medianRemoved.fits'
    vmin = 0; vmax = 3
    xdim = [190,230]; ydim = [130,180]
    folder = 'im24'
    save_name = folder+'/imgs/im24_'
    video_name = '24um_diffusion.avi'
    prefix = 'im24/imgs/im24_'
    start = 5
    stop = 5800 # Acutal number - > 5880/ Gonna try to keep it consistent. 
                # Things dont really move after 3000 anyway. 
                # If you wanted to make a longer video though youd do that here.
if im_24_ext_5: # extended by 5 pixels to the left
    pref_sub = 'im24_ext/snr24/24um_diff_'
    pref_bkgd = 'im24_ext/bkgd24/24um_diff_'
    original = '../../Sky_Remove/Median_Removed/24um_medianRemoved.fits'
    vmin = 0; vmax = 3
    xdim = [190,230]; ydim = [130,180]
    folder = 'im24_ext'
    save_name = folder+'/imgs/im24ext_'   
    video_name = '24um_extended_5_diffusion.avi'
    prefix = 'im24_ext/imgs/im24ext_'
    start = 5
    stop =  5800 # Acutal number - > 6660 
if im_24_ext_10: # extended by 5 pixels to the left
    pref_sub = 'im24_ext_10/snr24/24um_diff_'
    pref_bkgd = 'im24_ext_10/bkgd24/24um_diff_'
    original = '../../Sky_Remove/Median_Removed/24um_medianRemoved.fits'
    vmin = 0; vmax = 3
    xdim = [190,230]; ydim = [130,180]
    folder = 'im24_ext_10'
    save_name = folder+'/imgs/im24ext_10_'   
    video_name = '24um_extended_10_diffusion.avi'
    prefix = 'im24_ext_10/imgs/im24ext_10_'
    start = 5
    stop =  5800 # Acutal number - > 7690
if im_24_ext_20: # extended by 5 pixels to the left
    pref_sub = 'im24_ext_20/snr24/24um_diff_'
    pref_bkgd = 'im24_ext_20/bkgd24/24um_diff_'
    original = '../../Sky_Remove/Median_Removed/24um_medianRemoved.fits'
    vmin = 0; vmax = 3
    xdim = [190,230]; ydim = [130,180]
    folder = 'im24_ext_20'
    save_name = folder+'/imgs/im24ext_20_'   
    video_name = '24um_extended_20_diffusion.avi'
    prefix = 'im24_ext_20/imgs/im24ext_20_'
    start = 5
    stop =  1500 # Acutal number - > 7690   
if im_24_ext_40: # extended by 5 pixels to the left
    pref_sub = 'im24_ext_40/snr24/24um_diff_'
    pref_bkgd = 'im24_ext_40/bkgd24/24um_diff_'
    original = '../../Sky_Remove/Median_Removed/24um_medianRemoved.fits'
    vmin = 0; vmax = 3
    xdim = [190,230]; ydim = [130,180]
    folder = 'im24_ext_40'
    save_name = folder+'/imgs/im24ext_40_'   
    video_name = '24um_extended_40_diffusion.avi'
    prefix = 'im24_ext_40/imgs/im24ext_40_'
    start = 5
    stop =  6000 # Acutal number - > 7690  
if im_24_ext_60: # extended by 5 pixels to the left
    pref_sub = 'im24_ext_60/snr24/24um_diff_'   
    pref_bkgd = 'im24_ext_60/bkgd24/24um_diff_'
    original = '../../Sky_Remove/Median_Removed/24um_medianRemoved.fits'
    vmin = 0; vmax = 3
    xdim = [190,230]; ydim = [130,180]
    folder = 'im24_ext_60'
    save_name = folder+'/imgs/im24ext_60_'   
    video_name = '24um_extended_60_diffusion.avi'
    prefix = 'im24_ext_60/imgs/im24ext_60_'
    start = 5
    stop =  4195 # Acutal number - > 7690    
if im_70_ext_20:
    pref_sub = 'im70_ext_20/snr70/70um_diff_'
    pref_bkgd = 'im70_ext_20/bkgd70/70um_diff_'
    original = '../../Sky_Remove/Median_Removed/70um_medianRemoved.fits'
    vmin = 0; vmax = 30
    xdim = [210,275]; ydim = [200,260]
    folder = 'im70_ext_20'
    save_name = folder+'/imgs/im70ext_20_'
    video_name = '70um_extended_20_diffusion.avi'
    prefix = 'im70_ext_20/imgs/im70ext_20_'
    start = 5
    stop = 5990 # Acutal number - > 9995     
if im_70:
    pref_sub = 'im70/snr70/70um_diff_'
    pref_bkgd = 'im70/bkgd70/70um_diff_'
    original = '../../Sky_Remove/Median_Removed/70um_medianRemoved.fits'
    vmin = 0; vmax = 30
    xdim = [210,275]; ydim = [200,260]
    folder = 'im70'
    save_name = folder+'/imgs/im70_'
    video_name = '70um_diffusion.avi'
    prefix = 'im70/imgs/im70_'
    start = 5
    stop = 5800 # Acutal number - > 9995
if im_100_div_70:    
    pref_sub = 'ratio_100_to_70/snr100/100um_diff_'
    pref_bkgd = 'ratio_100_to_70/bkgd100/100um_diff_'
    original = '../../Sky_Remove/Median_Removed/100um_medianRemoved.fits'
    vmin = 0; vmax = 30
    xdim = [180,230]; ydim = [165,220]
    folder = 'ratio_100_to_70'
    save_name = folder+'/imgs/im100_'
    video_name = '100_divide_70_um_diffusion2.avi'
    prefix = 'ratio_100_to_70/imgs/im100_'
    start = 5
    stop = 6000
if im_160_div_70:    
    pref_sub = 'ratio_160_to_70/snr160/160um_diff_'
    pref_bkgd = 'ratio_160_to_70/bkgd160/160um_diff_'
    original = '../../Sky_Remove/Median_Removed/160um_medianRemoved.fits'
    vmin = 0; vmax = 2
    xdim = [110,145]; ydim = [105,135]
    folder = 'ratio_160_to_70'
    save_name = folder+'/imgs/im160_'
    video_name = '160_divide_70_um_diffusion.avi'
    prefix = 'ratio_160_to_100/imgs/im100_'
    start = 5
    stop =4515      
if im_160_div_100:    
    pref_sub = 'ratio_160_to_100/snr160/160um_diff_'
    pref_bkgd = 'ratio_160_to_100/bkgd160/160um_diff_'
    original = '../../Sky_Remove/Median_Removed/160um_medianRemoved.fits'
    vmin = 0; vmax = 2
    xdim = [110,145]; ydim = [105,135]
    folder = 'ratio_160_to_100'
    save_name = folder+'/imgs/im160_'
    video_name = '160_divide_100_um_diffusion.avi'
    prefix = 'ratio_160_to_100/imgs/im100_'
    start = 5
    stop =4515     
if im_100:
    pref_sub = 'im100/snr100/100um_diff_'
    pref_bkgd = 'im100/bkgd100/100um_diff_'
    original = '../../Sky_Remove/Median_Removed/100um_medianRemoved.fits'
    vmin = 0; vmax = 70
    xdim = [172,234]; ydim = [162,215]
    folder = 'im100'
    save_name = folder+'/imgs/im100_'
    video_name = '100um_diffusion.avi'
    prefix = 'im100/imgs/im100_'
    start = 5
    stop = 5245  
if im_160:
    pref_sub = 'im160/snr160/160um_diff_'
    pref_bkgd = 'im160/bkgd160/160um_diff_'
    original = '../../Sky_Remove/Median_Removed/160um_medianRemoved.fits'
    vmin = 0; vmax = 70
    xdim = [110,145]; ydim = [105,135]
    folder = 'im160'
    save_name = folder+'/imgs/im160_'
    video_name = '160um_diffusion.avi'
    prefix = 'im160/imgs/im160_'
    start = 5
    stop = 5800

app_sub = '_steps_snr.fits'
app_bkgd = '_steps_bkgd.fits'
colormap = 'gnuplot2'

steps_start = 5 # Starts at 5, need these for saving the imgs in batches
steps = 3000


import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from matplotlib.widgets import Slider
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy import units as u

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
    for i in range(steps_start,steps,5):
        with fits.open(pref_sub+str(i)+app_sub) as hdu_sub:
            sub_data.append(hdu_sub[0].data)
            sub_negs.append(count_negs(hdu_sub[0].data))
            neg_imgs.append(negative_img(hdu_sub[0].data))

        with fits.open(pref_bkgd+str(i)+app_bkgd) as hdu_bkgd:
            bkgd_data.append(hdu_bkgd[0].data)
            bkgd_negs.append(count_negs(hdu_bkgd[0].data))

        #hdu_sub = fits.open(pref_sub+str(i)+app_sub)
        #hdu_bkgd = fits.open(pref_bkgd+str(i)+app_bkgd)
        #sub_data.append(hdu_sub[0].data)
        #bkgd_data.append(hdu_bkgd[0].data)
        #sub_negs.append(count_negs(hdu_sub[0].data))
        #bkgd_negs.append(count_negs(hdu_bkgd[0].data))
        #hdu_sub.close()
    #hdu_bkgd.close()



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

    # Parameter Box
    r = np.arange(steps_start,steps,5)
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
        ax.clear()
        bx.clear()
        cx.clear()
        # Set image shapes
        ax.set_xlim(xdim)
        ax.set_ylim(ydim)
        bx.set_xlim(xdim)
        bx.set_ylim(ydim)
        cx.set_xlim(xdim)
        cx.set_ylim(ydim)

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
        r = np.arange(steps_start,steps,5)
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
        # Numbers here are from Karin's paper.
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
        #negs.append(count_negs(hdu[0].data))
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
    #def findPer2(): # Uses a slightly    
    
    plt.figure(1)
    
    plt.plot(x,sums)
    #ind_50 = findPer(sums,.5)
    #ind_90 = findPer(sums,.9)
    #ind_95 = findPer(sums,.95)
    #plt.scatter(x[ind_50],sums[ind_50],c='green', label='50%: '+str(round(sums[ind_50],0)))
    #plt.scatter(x[ind_90],sums[ind_90],c='yellow', label='90%: '+str(round(sums[ind_90],0)))
    #plt.scatter(x[ind_95],sums[ind_95],c='r', label='95%: '+str(round(sums[ind_95],0)))
    #plt.legend(bbox_to_anchor=(.7, .8), loc=2, borderaxespad=0.)
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
    #ind_50_ = findPer(pruned,.5)
    #ind_90_ = findPer(pruned,.9)
    #ind_75_ = findPer(pruned,.75)
    #plt.scatter(x[ind_50_],pruned[ind_50_],c='green', label='50%: '+str(round(pruned[ind_50_],0)))
    #plt.scatter(x[ind_75_],pruned[ind_75_],c='yellow', label='75%: '+str(round(pruned[ind_75_],0)))
    #plt.scatter(x[ind_90_],pruned[ind_90_],c='r', label='90%: '+str(round(pruned[ind_90_],0)))
    plt.xlabel('Steps')
    plt.ylabel('Sum of Values within a 22 arcsecond Radius')
    plt.title(folder)
    plt.legend(bbox_to_anchor=(.7, .8), loc=2, borderaxespad=0.)
    
    
    plt.figure(4)
    #lastSum = pruned[len(pruned)-1]
    lastSum = pruned[101]
    diffs = np.abs(pruned-lastSum)/lastSum * 100
    plt.plot(x,diffs)
    plt.xlabel('Steps')
    plt.ylabel('Percent Difference from Most Homogenized Within a 22 arcsecond Radius')
    plt.title(folder)
    #atol = 50
    #ind50 = findPer(diffs,.5) #np.where(np.isclose(diffs,50,atol=atol))[0][0]
    #ind90 = findPer(diffs,.9)#np.where(np.isclose(diffs,90,atol=atol))[0][0]
    #ind95 = findPer(diffs,.95)#np.where(np.isclose(diffs,95,atol=atol))[0][0]
    #plt.scatter(x[ind50],diffs[ind50],c='green', label=str(round(diffs[ind50],0))+'%')
    #plt.scatter(x[ind90],diffs[ind90],c='yellow', label=str(round(diffs[ind90],0))+'%')
    #plt.scatter(x[ind95],diffs[ind95],c='r', label=str(round(diffs[ind95],0))+'%')
    #plt.legend(bbox_to_anchor=(.7, .8), loc=2, borderaxespad=0.)
   # plt.xlim(1000,4000)
   # plt.ylim(0,100)
    
    plt.figure(5)
    origSum = np.sum(Prune(original))
    diffs_ = np.abs(pruned-origSum)/origSum * 100
    plt.plot(x,diffs_)
    plt.xlabel('Steps')
    plt.ylabel('Percent Difference from Original Within a 22 arcsecond Radius')
    #ind50_ = findPer(diffs_[::-1],.5)#np.where(np.isclose(diffs_,50,atol=atol))[0][0]
    #ind75_ = findPer(diffs_[::-1],.75)#np.where(np.isclose(diffs_,75,atol=atol))[0][0]
    #ind90_ = findPer(diffs_[::-1],.9)#np.where(np.isclose(diffs_,90,atol=atol))[0][0]
    #plt.scatter(x[ind50_],diffs_[ind50_],c='green', label=str(round(diffs_[ind50_],0))+'%')
    #plt.scatter(x[ind75_],diffs_[ind75_],c='yellow', label=str(round(diffs_[ind75_],0))+'%')
    #plt.scatter(x[ind90_],diffs_[ind90_],c='r', label=str(round(diffs_[ind90_],0))+'%')
    #plt.legend(bbox_to_anchor=(.7, .8), loc=2, borderaxespad=0.)
    plt.title(folder)
    
plt.show()
