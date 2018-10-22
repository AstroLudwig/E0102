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
sys.path.insert(0, '../Prune')
import prune as p
############
# Switches #
############
plot_1 = False
plot_2 = True

save_plot = False

c_dec = -72.03125; c_ra = 16.00875
d = 61 #kpc
arcs = 22
#############
# Colormaps #
#############

colormap = 'gnuplot2'
    
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
fs100 = '100/100um_70modeled_ForceSlope_prnd_snr.fits'


fo160 = '160/e0102_pacs160_new.fits'
fb160 = '160/160um_70modeled_ForceSlope_bkgd.fits'
fs160 = '160/160um_70modeled_ForceSlope_prnd_snr.fits'

# Convolved and Regrided And Pruned remnants
F24 = "24/Final_24_SNR_CR_Prnd.fits"
F70 = "70/Final_70_SNR_CR_Prnd.fits"
F100 = "100/Final_100_SNR_CR_Prnd.fits"
F160 = "160/160um_70modeled_ForceSlope_prnd_snr.fits"
Final = [F24,F70,F100,F160]


files = [fo24,fo70,fo100,fo160,fb24,fb70,fb100,fb160,fs24,fs70,fs100,fs160]

imgs = []; hdrs = []; wcs = [];


for i in range(12):
    imgs.append(fits.open(files[i])[0].data)
    hdrs.append(fits.open(files[i])[0].header)
    wcs.append(WCS(fits.open(files[i])[0].header))
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

x_range = []; x_r = []
y_range = []; y_r = []
for i in range(len(hdrs)):
    # Note that when transferring WCS to pix, 
    # the xrange is out of order. the yrange is fine.
    w = WCS(hdrs[i])
    a,b = w.all_world2pix(ra_range,dec_range,0)
    x_r.append([a[1],a[0]]); y_r.append([b[0],b[1]])
    x_range.append(np.asarray([int(round(a[1])),int(round(a[0]))]))
    y_range.append(np.asarray([int(round(b[0])),int(round(b[1]))]))

###############################
## Get Vmin/Vmax for Range   ##
###############################
n = 0
def getV(data,xrange,yrange):
    vmin = np.nanmin(data[yrange[0]+n:yrange[1]-n,xrange[0]+n:xrange[1]-n])
    vmax = np.nanmax(data[yrange[0]+n:yrange[1]-n,xrange[0]+n:xrange[1]-n])
    return vmin,vmax
# Get general vmins/vmaxs for a region around snr
vmins = []
vmaxs = []
for i in range(len(imgs)):
    newim = p.prune_extend(files[i],10)
    vmins.append(np.nanmin(newim))
    vmaxs.append(np.nanmax(newim))
vmins[8] = np.nanmin(imgs[8]); vmaxs[8] = np.nanmax(imgs[8])    
vmins[9] = np.nanmin(imgs[9]); vmaxs[9] = np.nanmax(imgs[9])
vmins[10] = np.nanmin(imgs[10]); vmaxs[10] = np.nanmax(imgs[10])
vmins[11] = np.nanmin(imgs[11]); vmaxs[11] = np.nanmax(imgs[11])
##########
## Plot ##
##########
if plot_1:
    # Initialize Plot
    fig = plt.figure(figsize=(15, 15)) 
    #Orginals
    O24 = plt.subplot2grid((4,3),(0,0),projection=wcs[0])
    O70 = plt.subplot2grid((4,3),(1,0),projection=wcs[0])
    O100 = plt.subplot2grid((4,3),(2,0),projection=wcs[0])
    O160 = plt.subplot2grid((4,3),(3,0),projection=wcs[0])
    #Background 
    B24 = plt.subplot2grid((4,3),(0,1),projection=wcs[0])
    B70 = plt.subplot2grid((4,3),(1,1),projection=wcs[0])
    B100 = plt.subplot2grid((4,3),(2,1),projection=wcs[0])
    B160 = plt.subplot2grid((4,3),(3,1),projection=wcs[0])
    #Subtracted Remnant
    S24 = plt.subplot2grid((4,3),(0,2),projection=wcs[0])
    S70 = plt.subplot2grid((4,3),(1,2),projection=wcs[0])
    S100 = plt.subplot2grid((4,3),(2,2),projection=wcs[0])
    S160 = plt.subplot2grid((4,3),(3,2),projection=wcs[0])
    # Store in Array
    plts = [O24,O70,O100,O160,B24,B70,B100,B160,S24,S70,S100,S160]
    # Plot loop for basic stuff
    #x_range =np.repeat()
    #vmins = np.repeat(0,12)
    #vmaxs = np.repeat(40,12)
    for i in range(12):

        P = plts[i].imshow(imgs[i],vmin=vmins[i],vmax=vmaxs[i])
        plts[i].set_xlim(x_r[i])
        plts[i].set_ylim(y_r[i])
        fig.colorbar(P, ax=plts[i],fraction=.03)
        
        if np.isin(i,[0,1,2,4,5,6,8,9,10]): 
            ra = plts[i].coords['ra']
            ra.set_ticklabel_visible(False)
        if np.isin(i,[4,5,6,7,8,9,10,11]):
            dec = plts[i].coords['dec']
            dec.set_ticklabel_visible(False)
    n = 4
    O24.set_xlim(x_r[0][0]-n,x_r[0][1]+n)
    B24.set_xlim(x_r[0][0]-n,x_r[0][1]+n)
    S24.set_xlim(x_r[0][0]-n,x_r[0][1]+n)
    #plt.subplots_adjust(left=0.077, right=0.969)
    plt.subplots_adjust(left=0.077, right=0.969, top=0.985,bottom=0.035,hspace=0,wspace=0.160)
    plt.savefig("BackgroundRemoval.png")


if plot_2:
    f, ((ax,bx),(cx,dx)) = plt.subplots(2,2)
    plts = [ax,bx,cx,dx]; tits = ["24 $\mu$m","70 $\mu$m","100 $\mu$m","160 $\mu$m"]
    for i in range(4):
        im = fits.open(Final[i])[0].data
        plts[i].imshow(im,vmin=0,vmax=9)
        plts[i].set_ylim(110,130)
        plts[i].set_xlim(115,140)
        plts[i].axis("off")
        plts[i].set_title(tits[i])

    plt.subplots_adjust(left=0.03, right=0.94, top=0.915,bottom=0.035,hspace=0.2,wspace=0.2)
    plt.savefig("ConvolvedFinalImgs.png")
plt.show()