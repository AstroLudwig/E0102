import matplotlib.pyplot as plt 
from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
from matplotlib import ticker
import matplotlib.cm as cm
from astropy.coordinates import SkyCoord
from astropy import units as u
import sys
sys.path.insert(0, '../Convolve_Regrid')
import Convolve
import regrid
############
# Switches #
############

prune = True
regrid_convolve = False
save_plot = False

c_dec = -72.03125; c_ra = 16.00875
d = 61 #kpc
arcs = 22
#############
# Colormaps #
#############

#colormap = 'gray'
#colormap = cm.greys_r
#colormap = 'bone'
colormap = 'gnuplot2'

#########
# Prune #
#########

# Take Final Supernova Files and Remove Everything Outside of a 22" Radius. 

def Prune(fname,sname):
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

    data[row[look],col[look]] = np.nan

    #fits.writeto(sname,data,hdr,overwrite=True)

    return data
if prune:
    Prune('24/24um_diff_3000_steps_snr.fits','24/24um_diff_3000_steps_prnd_snr.fits')
    Prune('70/70um_diff_3000_steps_snr.fits','70/70um_diff_3000_steps_prnd_snr.fits')
    Prune('100/100um_70modeled_ForceSlope_snr.fits','100/100um_70modeledFS_prnd_snr.fits')
    Prune('160/160um_70modeled_ForceSlope_snr.fits','160/160um_70modeledFS_prnd_snr.fits')
    
######################
# Load Data for Plot #
######################  

# Files
fo24 = '24/n76crop_smcsage_24_desaturated.fits'
fb24 = '24/24um_diff_3000_steps_bkgd.fits'
fs24 = '24/24um_diff_3000_steps_prnd_snr.fits'

fo70 = '70/e0102_pacs70_new.fits'
fb70 = '70/70um_diff_3000_steps_bkgd.fits'
fs70 = '70/70um_diff_3000_steps_prnd_snr.fits'

fo100 = '100/e0102_pacs100_new.fits'
fb100 = '100/100um_70modeled_ForceSlope_bkgd.fits'
fs100 = '100/100um_70modeledFS_prnd_snr.fits'


fo160 = '160/e0102_pacs160_new.fits'
fb160 = '160/160um_70modeled_ForceSlope_bkgd.fits'
fs160 = '160/160um_70modeledFS_prnd_snr.fits'

files = [fo24,fb24,fs24,fo70,fb70,fs70,fo100,fb100,fs100,fo160,fb160,fs160]
imgs = []; hdrs = []

for i in range(len(files)):
    imgs.append(fits.open(files[i])[0].data)
    hdrs.append(fits.open(files[i])[0].header)

##########################
# Regrid Convolve to 160 #
##########################  
def rc(kernel,img,regridto,sname):
    pref = '../Convolve_Regrid/Conv_Regrids/'
    Convolve.master_convolve(kernel,img,'../Convolve_Regrid/skip.fits')
    final = regrid.resample('../Convolve_Regrid/skip.fits',regridto,sname)
    return final

if regrid_convolve:
    k24_160 = '../Convolve_Regrid/Kernels/Kernel_LoRes_MIPS_24_to_PACS_160.fits'
    k70_160 = '../Convolve_Regrid/Kernels/Kernel_LoRes_PACS_70_to_PACS_160.fits'
    k100_160 = '../Convolve_Regrid/Kernels/Kernel_LoRes_PACS_100_to_PACS_160.fits'
    news24 = rc(k24_160,fs24,fo160,'final24_RC_160.fits')
    news70 = rc(k70_160,fs70,fo160,'final70_RC_160.fits')
    news100 = rc(k100_160,fs100,fo160,'final100_RC_160.fits')
##################
## Get Mask     ##
##################
# Masks
m24 = '../Xray/Mask_World_24_V2.txt' 
m70 = '../Xray/Mask_World_70.txt'   
m70 = '../Xray/Mask_World_100.txt'
m70 = '../Xray/Mask_World_160.txt'

def get_mask(mask_file,hdr):
    # Load mask from convolved and regridded xray file. 
    ra,dec = np.loadtxt(mask_file)
    # Convert mask from WCS to pixels
    w = WCS(hdr)
    x_m, y_m = w.all_world2pix(ra,dec,1)
    # Convert coordinates to integers. 
    xx = [int(i) for i in x_m]
    yy = [int(i) for i in y_m]
    
    return xx,yy

Masks_x = []
Masks_y = []


####################
## Set Plot Range ##
####################
def get_center(hdr):
    w = WCS(hdr)
    x,y = w.all_world2pix(c_ra,c_dec,1)
    x = int(round(x)); y = int(round(y))
    return x,y


# Set range based on center. 1/2 arcsecond up/down/left/right < - adjust
dec_range = [c_dec - 0.01, c_dec + 0.01]
ra_range = [c_ra - 0.05, c_ra + 0.05]

x_range = []
y_range = []
for i in range(len(hdrs)):
    # Note that when transferring WCS to pix, 
    # the xrange is out of order. the yrange is fine.
    w = WCS(hdrs[i])
    a,b = w.all_world2pix(ra_range,dec_range,1)
    x_range.append(np.asarray([int(round(a[1])),int(round(a[0]))]))
    y_range.append(np.asarray([int(round(b[0])),int(round(b[1]))]))

###############################
## Get Vmin/Vmax for Range   ##
###############################
n = 0 # to extend range 
def getV(data,xrange,yrange):

    x = np.arange(int(xrange[0]-n),int(xrange[1])+n)
    y = np.arange(int(yrange[0]-n),int(yrange[1])+n)
    x,y = np.meshgrid(x,y)
    x = x.flatten()
    y = y.flatten()
    grid = []
    for i in range(len(x)):
        grid.append(data[x[i],y[i]])
    vmin = np.nanmin(grid)
    vmax = np.nanmax(grid)
    return vmin,vmax

vmins = []
vmaxs = []
for i in range(len(imgs)):
    #if i % 2 == 0:
    va,vb = getV(imgs[i],x_range[i],y_range[i])
    print(va,vb)
    vmins.append(va)
    vmaxs.append(vb)
 #   elif i % 2 == 1: 
   #     vmins.append(np.nanmin(imgs[i]))
  #      vmaxs.append(np.nanmax(imgs[i]))
##################
## Plot Stuff   ##
##################
textsize = 10
def axoff(axis,n):
    ra = axis.coords['ra']
    dec = axis.coords['dec']
    if n == 0: # xaxis
        ra.set_ticks_visible(False)
        #ra.set_ticks_visible(True)
        ra.set_ticklabel_visible(False)
    if n == 1: # yaxis
        dec.set_ticks_visible(False)
        #dec.set_ticks_visible(True)
        dec.set_ticklabel_visible(False)

def axadjust(axis):
    ra = axis.coords['ra']
    dec = axis.coords['dec']
    ra.set_ticks(size=1)
    ra.set_ticklabel(size=textsize) 
    dec.set_ticks(size=1)
    dec.set_ticklabel(size=textsize)        


allvmin = 0       

# https://media.readthedocs.org/pdf/wcsaxes/latest/wcsaxes.pdf
#fig = plt.figure(figsize=(25.0, 20.0)) 
fig = plt.figure(figsize=(18.0, 21.0)) 
fig.set_tight_layout(True)

rows = 4; cols = 3
#################
## Originals   ##
#################
ex = 5; ex_x1 = 187 ; ex_x2 =234 
# 24 micron
ax = fig.add_subplot(rows,cols,1,projection=WCS(hdrs[0]))
cax = ax.imshow(imgs[0],origin='lower',vmin=allvmin,vmax=np.nanmax(imgs[2]),cmap=colormap)
ax.set_xlim(ex_x1,ex_x2)
ax.set_ylim(y_range[0])
cbar = fig.colorbar(cax)
cbar.ax.tick_params(labelsize= textsize)
axoff(ax,0)
axadjust(ax)
# 70 micron
bx = fig.add_subplot(rows,cols,4,projection=WCS(hdrs[3]))
cbx = bx.imshow(imgs[3],origin='lower',vmin=allvmin,vmax=vmaxs[3],cmap=colormap)
bx.set_xlim(x_range[3])
bx.set_ylim(y_range[3])
cbar = fig.colorbar(cbx)
cbar.ax.tick_params(labelsize= textsize)
axoff(bx,0)
#axoff(bx,1)
axadjust(bx)
#100 micron
cx = fig.add_subplot(rows,cols,7,projection=WCS(hdrs[6]))
ccx =cx.imshow(imgs[6],origin='lower',vmin=allvmin,vmax=vmaxs[6],cmap=colormap)
cx.set_xlim(x_range[6])
cx.set_ylim(y_range[6])
cbar = fig.colorbar(ccx)
cbar.ax.tick_params(labelsize= textsize)
axoff(cx,0)
axadjust(cx)
# 160 micron
dx = fig.add_subplot(rows,cols,10,projection = WCS(hdrs[9]))
cdx = dx.imshow(imgs[9],origin='lower',vmin=allvmin,vmax=vmaxs[9],cmap=colormap)
dx.set_xlim(x_range[9])
dx.set_ylim(y_range[9])
cbar = fig.colorbar(cdx)
cbar.ax.tick_params(labelsize= textsize)
cbar.locator = ticker.MaxNLocator(nbins=5)
cbar.update_ticks()
#axoff(dx,0)
#axoff(dx,1)
axadjust(dx)

############
## Bkgd   ##
############
# 24 micron
ex = fig.add_subplot(rows,cols,2,projection = WCS(hdrs[1]))
cex = ex.imshow(imgs[1],origin='lower',vmin=allvmin,vmax=vmaxs[1],cmap=colormap)
#ex.set_xlim(x_range[0][0]-ex,x_range[1][1])
ex.set_xlim(ex_x1,ex_x2)
ex.set_ylim(y_range[1])
cbar = fig.colorbar(cex)
cbar.ax.tick_params(labelsize= textsize)
axoff(ex,0)
axoff(ex,1)
axadjust(ex)
# 70 micron
fx = fig.add_subplot(rows,cols,5,projection = WCS(hdrs[4]))
cfx = fx.imshow(imgs[4],origin='lower',vmin=allvmin,vmax=vmaxs[4],cmap=colormap)
fx.set_xlim(x_range[4])
fx.set_ylim(y_range[4])
cbar = fig.colorbar(cfx)
cbar.ax.tick_params(labelsize= textsize)
cbar.locator = ticker.MaxNLocator(nbins=5)
cbar.update_ticks()
axoff(fx,0)
axoff(fx,1)
axadjust(fx)
# 100 micron
gx = fig.add_subplot(rows,cols,8,projection = WCS(hdrs[7]))
cgx = gx.imshow(imgs[7],origin='lower',vmin=allvmin,vmax=vmaxs[7],cmap=colormap)
gx.set_xlim(x_range[7])
gx.set_ylim(y_range[7])
cbar = fig.colorbar(cgx)
cbar.ax.tick_params(labelsize= textsize)
axoff(gx,0)
axoff(gx,1)
axadjust(gx)

# 160 microns
hx = fig.add_subplot(rows,cols,11,projection = WCS(hdrs[10]))
chx = hx.imshow(imgs[10],origin='lower',vmin=allvmin,vmax=vmaxs[10],cmap=colormap)
hx.set_xlim(x_range[10])
hx.set_ylim(y_range[10])
cbar = fig.colorbar(chx)
cbar.ax.tick_params(labelsize= textsize)

axoff(hx,1)
axadjust(hx)
#ra_hx = hx.coords['ra']
#ra_hx.set_major_formatter('hh:mm:ss')

####################
## Subtractions   ##
####################
# 24 microns
ix = fig.add_subplot(rows,cols,3,projection=WCS(hdrs[2]))
cix = ix.imshow(imgs[2],origin='lower',vmin=allvmin,vmax=np.nanmax(imgs[2]),cmap=colormap)
ix.set_xlim(ex_x1,ex_x2)
ix.set_ylim(y_range[2])
cbar = fig.colorbar(cix)
cbar.ax.tick_params(labelsize= textsize)
axoff(ix,0)
axoff(ix,1)
axadjust(ix)
# 70 microns
jx = fig.add_subplot(rows,cols,6,projection=WCS(hdrs[5]))
cjx = jx.imshow(imgs[5],origin='lower',vmin=allvmin,vmax=np.nanmax(imgs[5]),cmap=colormap)
jx.set_xlim(x_range[5])
jx.set_ylim(y_range[5])
cbar = fig.colorbar(cjx)
cbar.ax.tick_params(labelsize= textsize)
axoff(jx,0)
axoff(jx,1)
axadjust(jx)
# 100 microns
kx = fig.add_subplot(rows,cols,9,projection=WCS(hdrs[8]))
ckx =kx.imshow(imgs[8],origin='lower',vmin=allvmin,vmax=np.nanmax(imgs[8]),cmap=colormap)
kx.set_xlim(x_range[8])
kx.set_ylim(y_range[8])
cbar = fig.colorbar(ckx)
cbar.ax.tick_params(labelsize= textsize)
axoff(kx,0)
axoff(kx,1)
axadjust(kx)
# 160 microns
lx = fig.add_subplot(rows,cols,12,projection = WCS(hdrs[11]))
clx = lx.imshow(imgs[11],origin='lower',vmin=allvmin,vmax=np.nanmax(imgs[11]),cmap=colormap)
lx.set_xlim(x_range[11])
lx.set_ylim(y_range[11])
cbar = fig.colorbar(clx)
cbar.ax.tick_params(labelsize= textsize)
axoff(lx,1)
axadjust(lx)

fig.subplots_adjust(hspace=0.1)
fig.subplots_adjust(wspace=.01)

plt.savefig("plot.png")

plt.show()