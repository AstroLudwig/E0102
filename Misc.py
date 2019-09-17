import numpy as np 
from astropy.wcs import WCS
from astropy.io import fits

# Coordinate Transforms
def to_world(fname,pix_x,pix_y):
    return WCS(fname).all_pix2world(pix_x,pix_y,0)
def to_pix(fname,ra,dec):
    return WCS(fname).all_world2pix(ra,dec,0)

# Kernel Adjustments
def check_kernel(kernel_fname):
    if type(kernel_fname) is str:
        kernel = fits.open(kernel_fname)[0].data
    else:
        kernel = kernel_fname
    # Square
    if kernel.shape[0] != kernel.shape[1]:
        print("Additional Code Required. Square the Kernel.")
    # Max in Center
    ind = np.unravel_index(np.argmax(kernel),np.shape(kernel))    
    if ind[0] != ind[1] and np.abs(ind[0] / kernel.shape[0] - 0.5) > 0.002:
        print("Additional Code Required. Max value not centered.")
        print(np.abs(ind[0] / kernel.shape[0] - 0.5))
    # Normalize
    if np.nansum(kernel) != 1:
        kernel = kernel / np.nansum(kernel) 
    return kernel