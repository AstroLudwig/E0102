# -*- Copyright (c) 2018, Bethany Ann Ludwig, All rights reserved. -*-

#################################################
## ******************************************  ##
##  ~         Regridding Routine            ~  ##
## ******************************************  ##
#################################################

from astropy.io import fits
from reproject import reproject_interp
from reproject import reproject_exact 
from scipy.ndimage import zoom
import numpy as np
import Convolve 

###########################################################################
# Regrid is a script used to resample a particular file to another file's #
# dimensions.                                                             #
# The master function in this case is regrid. It takes two file locations #
# and a save name.                                                        #
# Known Issues:                                                           #
# ---> ValueError:                                                        #
# "Output WCS has celestial components but input WCS does not"            #
# Something went wrong when saving one of the headers for the             # 
# input or output files. Double check that the headers are corect.        # 
# ---> AttributeError:                                                    #
# "module '__main__' has no attribute '__spec__'"                         #
# This message is kind of cryptic. It ocassionally happens when using     #
# reproject_exact. I haven't been able to dermine the cause, but          #  
# restarting or reinstalling the package usually seems to fix it.         #
###########################################################################

# Regrid Function
def regrid_exact(input_string, output_string):
    
    projection = fits.open(output_string)[0].header
    data_arr = fits.open(input_string)[0].data
    data_hdr = fits.open(input_string)[0].header
    
    array,fp = reproject_exact((data_arr,data_hdr),projection,hdu_in=0)
    
    return array

# Master Function 
def resample(input_file,output_file,sname): 
    arr = regrid_exact(input_file,output_file)
    newfile = arr 

    fits.writeto(sname,newfile,fits.open(output_file)[0].header,overwrite=True)
    return newfile
