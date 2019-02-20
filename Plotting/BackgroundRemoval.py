# -*- Copyright (c) 2019, Bethany Ann Ludwig, All rights reserved. -*-
"""
NAME:
	
PURPOSE:
	
"""
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt 
from astropy.wcs import WCS



xdim = [114,142]
ydim = [108,132]
xlim = [[114,142],[114,142],[114,142],[114,142]]
ylim = [[108,132],[108,132],[108,132],[108,132]]
title = ["24um","70um","100um","160um"]

f160 = "../Final_Files/160/160_SNR_Prune.fits"
def get_data(filename):
	return fits.open(filename)[0].data
def get_header(filename):
	return fits.open(filename)[0].header
def get_lim(filename):
	wcs = WCS(get_header(f160))
	ra,dec =  wcs.all_pix2world(xdim,ydim,1) 
	new_wcs = WCS(get_header(filename))
	x,y = new_wcs.all_world2pix(ra,dec,1)
	return [int(round(i)) for i in x],[int(round(i)) for i in y]
def adjust_lim(filename):
	x,y = get_lim(filename)
	width = x[1] - x[0]
	img_width = np.shape(get_data(filename))[1]
	ref_width = np.shape(get_data(f160))[1]
	standard_width = 28 * (ref_width/img_width)
	print(standard_width)
	if width > standard_width: 
		adjust = (width - standard_width ) / 2 
		print(adjust)
		x[0] = x[0]	+ adjust
		x[1] = x[1]	- adjust
	print(x)
	return x,y
def get_min(filename):
	data = get_data(filename)
	xdim,ydim = get_lim(filename)
	return np.nanmin(data[ydim[0]:ydim[1],xdim[0]:xdim[1]])
def get_max(filename):
	data = get_data(filename)
	xdim,ydim = get_lim(filename)
	return np.nanmax(data[ydim[0]:ydim[1],xdim[0]:xdim[1]])	
################
# Import Files #
################

# Originals 

original_files = ["../Original_Files/Infrared/n76crop_smcsage_24_desaturated.fits",
					"../Original_Files/Infrared/e0102_pacs70_new.fits",
						"../Original_Files/Infrared/e0102_pacs100_new.fits",
							"../Original_Files/Infrared/e0102_pacs160_new.fits"]
original_data = [get_data(x) for x in original_files ]
# Background 
background_files = ["../Final_Files/24/24um_diff_1375_steps_bkgd.fits",
						"../Final_Files/70/70um_diff_3620_steps_bkgd.fits",
							"../Final_Files/100/100um_70modeled_bkgd.fits",
								"../Final_Files/160/160um_70modeled_bkgd.fits"]
background_data = [get_data(x) for x in background_files ]								
# SNR  
subtraction_files = ["../Final_Files/24/24_SNR_Convolve_Regrid_Prune.fits",
						"../Final_Files/70/70_SNR_Convolve_Regrid_Prune.fits",
							"../Final_Files/100/100_SNR_Convolve_Regrid_Prune.fits",
								"../Final_Files/160/160_SNR_Prune.fits"]
subtraction_data = [get_data(x) for x in subtraction_files ]	

data = np.asarray(original_data + background_data + subtraction_data).reshape(3,4)
files = np.asarray(original_files + background_files + subtraction_files).reshape(3,4)
wcs = WCS(get_header(subtraction_files[3]))



counter = 1
fig = plt.figure(figsize=(11,8))
for i in range(4):
	for j in range(3):
		ax = fig.add_subplot(4,3,counter,projection=wcs)
		ax.text(.5,.5,str(counter))
		img = ax.imshow(get_data(files[j,i]),vmin=get_min(files[j,i]),vmax=get_max(files[j,i]))

		ax.set_xlim(get_lim(files[j,i])[0])
		print((get_lim(files[j,i])[0]))
		ax.set_ylim(get_lim(files[j,i])[1])
		print(get_lim(files[j,i])[1])





		if np.isin(counter,[1,2,3,4,5,6,7,8,9]):
			ax.tick_params(axis='x',which='both',bottom=False,top=False,labelbottom=False)

		if np.isin(counter,[2,3,5,6,8,9,11,12]):
			ax.tick_params(axis='y',which='both',bottom=False,top=False,labelbottom=False)	

		if np.isin(counter,[10,11,12]):
		 	ax.tick_params(axis='x',which='both',top=False,bottom=True)		

		if np.isin(counter,[1,4,7,10]):
		 	ax.tick_params(axis='y',which="both",left=True,right=False)
					

		counter += 1

plt.show()		


