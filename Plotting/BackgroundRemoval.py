# -*- Copyright (c) 2019, Bethany Ann Ludwig, All rights reserved. -*-
"""
NAME:
	
PURPOSE:
	
"""
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt 


def get_data(filename):
	return fits.open(filename)[0].data
def get_header(filename):
	return fits.open(filename)[0].header

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
# Diffusion stopped at 1375th step for 24 microns and at the 3620th step for 70 microns

# SNR  