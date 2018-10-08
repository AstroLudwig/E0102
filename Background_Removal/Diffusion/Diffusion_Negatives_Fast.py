import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from matplotlib.widgets import Slider

steps_start = 5
steps = 9995
s = np.arange(steps_start,steps,5)

im24 = False
im24_5 = False
im24_10 = False
im24_20 = False
im70 = False
im70_20 = True

if im24:
    pref_sub = 'im24/snr24/24um_diff_'
    app_sub = '_steps_snr.fits'
    xdim = [190,230]; ydim = [130,180]
if im24_5:
    pref_sub = 'im24_ext/snr24/24um_diff_'
    app_sub = '_steps_snr.fits'
    xdim = [190,230]; ydim = [130,180] 
if im24_10:
    pref_sub = 'im24_ext_10/snr24/24um_diff_'
    app_sub = '_steps_snr.fits'
    xdim = [190,230]; ydim = [130,180] 
if im24_20:
    pref_sub = 'im24_ext_20/snr24/24um_diff_'
    app_sub = '_steps_snr.fits'
    xdim = [190,230]; ydim = [130,180]     
if im70:
    pref_sub = 'im70/snr70/70um_diff_'
    app_sub = '_steps_snr.fits'
    xdim = [210,275]; ydim = [200,260] 
if im70_20:
    pref_sub = 'im70_ext_20/snr70/70um_diff_'
    app_sub = '_steps_snr.fits'
    xdim = [210,275]; ydim = [200,260] 

def count_negs(data):
    x = np.arange(xdim[0],xdim[1])
    y = np.arange(ydim[0],ydim[1])
    x,y = np.meshgrid(x,y)
    newdata = data[y,x]

    a,b = np.where(np.isnan(newdata))
    newdata[a,b] = 0.
    
    e,count = np.shape(np.where(newdata< 0))
    c,d = np.where(newdata < 0)
    sum_n = np.sum(np.abs(newdata[c,d]))
   # print(count)
    g,h = np.where(newdata != 0)
    #print (newdata[g,h])
    return count, sum_n, np.mean(newdata)

pix_count = (xdim[1]-xdim[0])*(ydim[1]-ydim[0])


neg_sums = []; neg_counts = [];totals = []; mean_region = []
for i in range(steps_start,steps,5):
    with fits.open(pref_sub+str(i)+app_sub) as hdul:
        neg_count,neg_sum, avg = count_negs(hdul[0].data)
        mean_region.append(avg)
        neg_sums.append(neg_sum)
        neg_counts.append(neg_count)
        totals.append(np.nansum(hdul[0].data))

f, (ax,bx,cx) = plt.subplots(1,3)
ax.plot(s,neg_sums,c='b')
ax.plot(s,totals,c='Orange')
ax.grid
ax.set_title("Total Sum(Orange), Negative Absolute Sum (Blue)")
ax.set_xlabel("Time Steps")
bx.plot(s,neg_counts)
bx.set_title("Total Sum(Orange), Negative Count (Blue)")
bx.set_xlabel("Time Steps")
bx.plot(s,totals,c='Orange')
cx.plot(s,neg_counts)
cx.set_title(("Negative Count, {:2f} % max negative pixels in region").format(np.max(neg_counts)/pix_count*100))

"""
f, (ax,bx) = plt.subplots(1,2)
ax.plot(s,sum_negs_24)
bx.plot(s,sum_negs_70)
ax.set_title("Negatives in 24")
bx.set_title("Negatives in 70")
ax.set_ylabel("Number of Negatives in Subtraction")
ax.set_xlabel("Timesteps")

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
"""
plt.show()