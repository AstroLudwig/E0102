import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from matplotlib.widgets import Slider

steps = 6000
s = np.arange(5,steps,5)

def count_negs(data):
    n =  len(data[data<0])
    val = data[data<0]
    return n,val 

def neg_check(data):  
    e,f = np.shape(data)
    coord_x = [];coord_y = []; neg = []; 
    for i in range(e):
        for j in range(f):
            if np.isfinite(data[i,j]) and data[i,j] < 0.:
                coord_x.append(j)
                coord_y.append(i)
                neg.append(data[i,j])
    return neg

pref_sub_24 = 'Sky_Remove/Diff/All_steps/snr_24/24um_diff_'
pref_bkgd_24 = 'Sky_Remove/Diff/All_steps/bkgd_24/24um_diff_'
im_24 = 'Sky_Remove/Median_Removed/24um_medianRemoved.fits'

pref_sub_70 = 'Sky_Remove/Diff/All_steps/snr_70/70um_diff_'
pref_bkgd_70 = 'Sky_Remove/Diff/All_steps/bkgd_70/70um_diff_'
im70 = 'Sky_Remove/Median_Removed/70um_medianRemoved.fits'

app_sub = '_steps_snr.fits'
app_bkgd = '_steps_bkgd.fits'

sub_24 = []; bkgd_24 = []; sum_24 = []; sum_negs_24 = []; neg_val_24 = []
sub_70 = []; bkgd_70 = []; sum_70 = []; sum_negs_70 = []; neg_val_70 = []

for i in range(5,steps,5):
    sub_24_data = fits.open(pref_sub_24+str(i)+app_sub)[0].data
    bkgd_24_data = fits.open(pref_bkgd_24+str(i)+app_bkgd)[0].data
    sub_24.append(sub_24_data)
    bkgd_24.append(bkgd_24_data)
    n24,v24 = count_negs(sub_24_data)
    neg_val_24.append(np.nansum(v24))
    sum_negs_24.append(n24)
    sum_24.append(np.nansum(sub_24_data))
  
    
    sub_70_data = fits.open(pref_sub_70+str(i)+app_sub)[0].data
    bkgd_70_data = fits.open(pref_bkgd_70+str(i)+app_bkgd)[0].data
    sub_70.append(sub_70_data)
    bkgd_70.append(bkgd_70_data)
    n70,v70 = count_negs(sub_70_data)
    sum_negs_70.append(n70)
    sum_70.append(np.nansum(sub_70_data))
    neg_val_70.append(np.nansum(v70))

f, (ax,bx) = plt.subplots(1,2)
ax.plot(s,sum_negs_24)
ax.plot(s,sum_24)
bx.plot(s,sum_negs_70)
bx.plot(s,sum_70)
ax.set_title("Negatives in 24")
bx.set_title("Negatives in 70")
ax.set_ylabel("Number of Negatives in Subtraction")
ax.set_xlabel("Timesteps")


h, (ex,fx) = plt.subplots(1,2)
ex.plot(s,sum_24)
ex.plot(s,neg_val_24)
fx.plot(s,sum_70)
fx.plot(s,neg_val_70)

g, (cx,dx) = plt.subplots(1,2)

dum, vals_0_24 = count_negs(sub_24[0])
dum, vals_0_70 = count_negs(sub_70[0])
bins= 20
cx.hist(vals_0_24,bins=bins)
dx.hist(vals_0_70,bins=bins)
cx.set_title("Negatives in 24")
dx.set_title("Negatives in 70")
cx.set_xlabel("Negative Values")
cx.set_ylabel("Count")

def update(val):
    print("Updating...")
    num = int(slide.val)
    print(num)
    # Get Negs
    dum, vals_24 = count_negs(sub_24[num])
    dum, vals_70 = count_negs(sub_70[num])
    # Clear image
    cx.clear()
    dx.clear()


    # Fill images in
    cx.hist(vals_24,bins=bins)
    dx.hist(vals_70,bins=bins)
    
    cx.set_title("Negatives in 24")
    dx.set_title("Negatives in 70")
    cx.set_xlabel("Negative Values")
    cx.set_ylabel("Count")


    # Draw out
    g.canvas.draw_idle()


# Create Sliders
axcolor = 'aliceblue'
length = len(np.arange(1,steps,5)) - 1
ax_y = plt.axes([0.12, 0.01, 0.8, 0.03], facecolor=axcolor)
slide = Slider(ax_y, 'File: ', 0, length, valinit=0, valfmt='%0.0f')
slide.on_changed(update)
plt.show()