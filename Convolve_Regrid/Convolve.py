import numpy as np
from astropy.io import fits
from astropy.convolution import convolve_fft
import scipy.ndimage
from reproject import reproject_interp

#################################################
## ******************************************  ##
##  ~         Convolution Routine           ~  ##
## ******************************************  ##
#################################################

test = True

"""
This is a convolution script which allows you to adjust the resolution of one 
image to match another. 

I wrote a master function at the end which combines all functions into one. 
You only need the names of the files you want to use and a name you want
the convolved file saved to.
"""

#######################
##  The Functions    ##
#######################

# The Kernel has a few properties that need to be checked out. 
# The sum of the kernel should always be one. (kerncheck)
# The max value of the kernel needs to be dead center. 
# So the array itself has to be odd. (square) 
# kernel variable refers to the actual array of the kernel 

#Ensures Kernel is equal to 1
def kerncheck(kernel):
    if np.nansum(kernel) != 1.:
        correctedker = kernel/np.nansum(kernel)
    else:
        correctedker = kernel
    return correctedker

#makes sure kernel is odd with the max value in the center
def square(ker):
    a,b = np.shape(ker)
    ind = np.unravel_index(np.argmax(ker),np.shape(ker))
    newKer = np.copy(ker)
    
    print('Initial kernal shape : '+str((a,b)))
    print('Initial max index : '+str((ind)))

    # The kernel matrix should be square. This appends
    # zeros to the appropriate place if it is not square.
    if a!=b:
        print('Squaring...') 
    while a != b:
        
        if a > b: 
            newKer = np.append(newKer,np.zeros((np.shape(newKer)[0],1)),axis=1)
        elif a < b: 
            newKer = np.append(newKer,[np.zeros(np.shape(newKer)[1])],axis=0)
        # need new a,b values that now we changed the array.
        a,b = np.shape(newKer)
        
    # Once the kernel is square it needs to be odd as well. 
    if (a%2) != 1: # If not odd, adds zeros to both sides.
        print('Oddifying...') 
        newKer = np.append(newKer,[np.zeros(np.shape(newKer)[1])],axis=0)
        newKer = np.append(newKer,np.zeros((np.shape(newKer)[0],1)),axis=1)
    
    # Kernel is odd and square. Now the max value needs to be in the center.     
    ind = np.unravel_index(np.argmax(newKer),np.shape(newKer))
    if ind[1] != np.floor(a/2): #need to get ind as an actual number     
        
        diff = int(np.abs(ind[1] - np.floor(a/2)))
        newKer = np.roll(newKer,diff, axis = 1)
    
    ind = np.unravel_index(np.argmax(newKer),np.shape(newKer))
    if ind[0] != np.floor(b/2):
        diff = int(np.abs(ind[0] - np.floor(b/2)))
        newKer = np.roll(newKer,diff, axis = 0)
    
    # need final dimension values. 
    c,d = np.shape(newKer)
    e,f = np.unravel_index(np.argmax(newKer),np.shape(newKer))
    # series of printed checks to make sure everything makes sense when you convolve. 
    print('Kernel shape after squaring: '+str((c,d)))
    if c % 2 != 1 and d % 2 != 1:
        print('WARNING: Kernel dimensions are not odd.')
    print('Kernel index of max value after squaring: '+ str((e,f)))
    if e != np.floor(c/2) and f != np.floor(d/2): # Note that you have to subtract 1 because the 0th index is counted!!
        print('WARNING: Max value may not be centered.')
    print('Kernel value at middle index: '+str(newKer[e,f]))
    print('Kernel max values ' + str(np.nanmax(newKer)))
    
    return newKer

# Get the cdelts out of the header, check that they are close to equal.
# This gives the pixel shape for the science images.  
def getscaleimg(fname):
    cdelt1 = np.absolute(fname[0].header['CDELT1']) #units: degrees/pixel
    cdelt2 = np.absolute(fname[0].header['CDELT2'])
   # print(("Cdelt 1: {} Cdelt 2: {}").format(cdelt1,cdelt2))
    # may decide in the future that the scale is a ratio rather then equivalant
    scale_x = cdelt1 
    scale_y = cdelt2 
    if np.abs(scale_x - scale_y) > 0.000001: 
        print("Warning: cdelts are non matching.")
        print("Cdelt1 " + str(scale_x) + " Cdelt2 " +str(scale_y)+ " Difference "+ str(scale_x - scale_y))
    return scale_x, scale_y


# Get the cds out of the header, check that they are close to equal. 
# There may be an infinitesimal difference in the numbers. 
# Just got to make sure it isn't more then that and that two of the values are zero.
# It could imply a rotation matrix otherwise which could be a problem for this method.  
def getscalekern(kernname):
    cd1 = np.absolute(kernname[0].header['CD1_1']) #units: degrees/pixel
    cd2 = np.absolute(kernname[0].header['CD2_2'])
    cd1_2 = np.absolute(kernname[0].header['CD1_2']) # should be zero
    cd2_1 = np.absolute(kernname[0].header['CD2_1']) # should be zero
    scale_x = cd1    
    scale_y = cd2  
    if np.abs(scale_x - scale_y) > 0.000001:    
        print("Warning: CDs are non matching.")
        print("Cd1 " + str(cd1) + " Cdelt2 " +str(cd2)+ " Difference "+ str(cd1 - cd2))
    if cd1_2 & cd2_1 != 0:
        print("Warning: CD1_2 and CD2_1 are nonzero. Rotational Matrix may need to be accounted for.")
    return scale_x, scale_y

# Matches pixel scale of kernel to the image needed.
# Need the name of the imported file, before data/hdr seperation
def resample(kern, kernelfile,imgfile,test): 
    scalekern_x, scalekern_y = getscalekern(kernelfile) 
    scaleimg_x, scaleimg_y =  getscaleimg(imgfile)
    
    # The scaling ratio is chosen to be img / kern.     
    scaling = scaleimg_x / scalekern_x         
    print('Resampling with scaling factor: '+str(scaling))
    
    # Resampling:
    # A spline of degree 1 is piecewise linear!
    # 1/scaling is important here. Not well listed in the docs.
    newfile = scipy.ndimage.zoom(kern, 1/scaling, order=1) 
    
    print(("Scale for Kernel: {} Scale for img: {} New Kernel Scale: {}").format(scalekern_x,scaleimg_x,scalekern_x*(np.shape(kern)[0]/np.shape(newfile)[0]) ))
    # Left an option here where you can test to see if resampling makes sense
    # Outputs a small sample range before and after resampling.
    if test: 
        print('Before:') # Print out a test view before resampling
        ind_x, ind_y = np.unravel_index(np.argmax(kern),np.shape(kern))
        print(kern[(ind_x-2):(ind_x+3), (ind_x-2):(ind_x+3)])  
        print('After:') # Print out a test view after resampling
        ind_ax, ind_ay = np.unravel_index(np.argmax(newfile),np.shape(newfile))
        print(newfile[(ind_ax-2):(ind_ax+3), (ind_ax-2):(ind_ax+3)])  
        
    
    return newfile

# Saves file
def savefile(name,file,header):
    #fits.writeto(name, file, clobber=True)
    fits.writeto(name, file, header, overwrite=True)
    print("File saved.")

# Updates header with new scaling
def newhead(scalestr,namestr,file,data,header):
    file[0].header['CDELT1'] = scalestr    
    file[0].header['CDELT2'] = scalestr
    fits.writeto(namestr+'.fits', data, header, clobber = True)
    print("Header updated.")


#######################
##  Master Function  ##
#######################

# All you need for this function is the string of the file names
# and the name of the file you want saved. 
# Everything else is done for you.     
    
def master_convolve(kernelfilename, sciencefilename, savename):
    # Step 0: Admin 
    kernel_file = fits.open(kernelfilename)
    kernel = kernel_file[0].data
    science_file = fits.open(sciencefilename) 
    science = science_file[0].data

    # Step 1: Regrid the kernel to match pixel scales if necessary.
    a,b = getscalekern(kernel_file)
    c,d = getscaleimg(science_file)
    if a - c == 0 and b - d == 0:  
        print ('Pixel Scales Match.')  
        print (' ')
    else:
        kernel = resample(kernel,kernel_file,science_file,'False')
        print (' ')
        
    # Step 2 : Check that the kernel sums to 1. Ok to do this more then once. 
    kernel = kerncheck(kernel)
    # Step 3: Check that the kernel is odd with max value in the center. 
    kernel = square(kernel)
    # Step 4 : Double check kernel sum. 
    kernel = kerncheck(kernel)
    print("Kernel sum: "+str(np.nansum(kernel)))
    if np.isclose(np.nansum(kernel),1.) == False:
        print('WARNING: Kernel sum was not equal to 1.')
    # Step 5: Convolve.
    convolved_file = convolve_fft(science,kernel,allow_huge=True)
    # Step 6: Save convolved file.
    savefile(savename,convolved_file, science_file[0].header)
    # Step 7: Return convolved file as variable. 
    return convolved_file
 
