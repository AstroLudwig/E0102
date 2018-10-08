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

# I included a few files that can be used to test the functions. 
f1 = 'Background/LinSpace/Background/24micron_avgd_avg.fits' 
f2 = '../Original_Files/e0102_pacs100_comb.fits'
f3 = '../Original_Files/e0102_pacs160_comb.fits' 

# Quickest function but not as accurate as regrid_exact
def regrid(input_string, output_string):
    projection = fits.open(output_string)[0].header
    data = fits.open(input_string)[0]
    
    array,fp = reproject_interp(data,projection,order='bilinear')
        
    return array

# 24 microns doesn't want to convolve to 160 so there may be something going on
# with the kernel not regridding correctly.
# This function uses another method to regrid.
def regrid_kern(kern_string, img_string):
    i_x,i_y = Convolve.getscaleimg(img_string)
    k_x,k_y = Convolve.getscalekern(kern_string)
    k = fits.open(kern_string)
    array = zoom(k[0].data,i_x/k_x,order=1)    
    return array

# More accurate but takes longer
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

# When regridding it's necessary to 
# multiply by the array by the area change
# If the units are in JY. Including a function to handle this.  
# Get scale img retrieves the dimensions of a pixel 
# Resample combines getscaleimg and regrid_exact to return a correctly
# Resampled Array 
def getscaleimg(filename):
    fname = fits.open(filename)
    cdelt1 = np.absolute(fname[0].header['CDELT1']) #units: degrees/pixel
    cdelt2 = np.absolute(fname[0].header['CDELT2'])
    # may decide in the future that the scale is a ratio rather then equivalant
    scale_x = cdelt1 
    scale_y = cdelt2 
    if np.abs(scale_x - scale_y) > 0.000001: 
        print("Warning: cdelts are non matching.")
        print("Cdelt1 " + str(scale_x) + " Cdelt2 " +str(scale_y)+ " Difference "+ str(scale_x - scale_y))
    return scale_x, scale_y

# Master Function For Units of Jy 
def resample_JY(input_file,output_file,sname): 
    old_x, old_y =  getscaleimg(input_file)
    new_x, new_y = getscaleimg(output_file)
    
    arr = regrid(input_file,output_file)
    newfile = arr * np.power(new_x/old_x,2)

    fits.writeto(sname,newfile,fits.open(output_file)[0].header,overwrite=True)
    return newfile