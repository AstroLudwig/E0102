from mpl_toolkits.mplot3d import Axes3D
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt 
from matplotlib.widgets import Slider,RadioButtons

##############
# Functions #
##############
def getMap(ChiMap,interval):
	shell = np.copy(ChiMap)
	i,j =  np.where(shell>np.min(shell)+interval)
	shell[i,j] = np.nan
	return shell
def getRange(pixel):
	# Set Parameter Ranges
	T = np.arange(2,70,1) # Kelvin
	M = 10**np.arange(-4,0.1,.1)

	leftMost=0; rightMost=0;
	topMost=0; bottomMost=0;

	firstFinitePixel = True
	for i in range(len(T)):
		for j in range(len(M)):
			if np.isfinite(pixel[i][j]):
				if firstFinitePixel:
					leftMost = i
					rightMost = i
					topMost = j
					bottomMost = j
					firstFinitePixel = False
				if i < leftMost:
					leftMost = i
				if i > rightMost:
					rightMost = i
				if j < bottomMost:
					bottomMost = j
				if j > bottomMost:
					bottomOst = j

	tempRange = T[rightMost] - T[leftMost]
	massRange = M[topMost] - M[bottomMost]
	return tempRange, massRange
	
#########
# Main #
#########

# Load Chi Map for AMC Full
chiMap = np.load("Sols/Parameter_Chi_Map.npy")

# Initialize Arrays for Various Confidence Levels
chimap = []
chimap99 = []; chimap95 = []; chimap90 = []; chimap68 = []
top = []; bottom = []; left = []; right = [] 
# Keep track of what pixel we're on. 
Px = []; Py = []
for i in range(244):
	for j in range(254):
		# Only look at valid pixels.
		if np.sum(chiMap[i,j]) > 0.:	
			# Initialize Maps
			chimap99.append(getMap(chiMap[i,j],6.63))
			chimap95.append(getMap(chiMap[i,j],4))
			chimap90.append(getMap(chiMap[i,j],2.71))
			chimap68.append(getMap(chiMap[i,j],1)) 
			chimap.append(chiMap[i,j])
			Px.append(j);Py.append(i);
np.savetxt("Sols/achi.txt",chimap[0])
T = np.arange(2,70,1) # Kelvin
M = 10**np.arange(-4,0.1,.1)
T,M = np.meshgrid(T,M)

zmin = np.nanmin(chimap[0])

T = np.arange(2,70,1) # Kelvin
M = 10**np.arange(-4,0.1,.1)

f, (ax,bx) = plt.subplots(1,2)
extent = (M[0],M[len(M)-1],T[0],T[len(T)-1])
Try = map(chimap68,(M,T))
ax.imshow(chimap68[0])

ax.contour(chimap[0], levels=(zmin+1,zmin+2.71,zmin+4,zmin+6.63))
ax.set_ylim(20,50)
ax.set_xlim(0,20)


plt.show()