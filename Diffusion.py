import numpy as np 

# Diffusion
def ApplyBoundaryConditions(original_data, new_data, index_set, timestep) :
    for indexPair in index_set :
        # Converts between index in original image and index in smaller template img. 
        # Diffusion occurs in smaller template image and then gets replaced in the 
        i = int(indexPair[1]) - minY
        j = int(indexPair[0]) - minX
        new_data[i,j,timestep] = original_data[int(indexPair[1]), int(indexPair[0])]
def StepForward(data, timestep) :
    xLim, yLim, junk = np.shape(data)
    # Kappa should be small to avoid converging too quicly, 
    # but is arbitrary as we are using the 
    # homoginzed solution as a final stopping place. 
    alpha = 0.1
    for i in range(1, xLim-1) :
        for j in range(1, yLim-1) :
            # Discretizing heat equation
            # Backward Finite Difference in Time
            # Second Order Central Difference in Space
            data[i,j,timestep] = data[i,j,timestep-1] + alpha * (
                    (data[i+1,j,timestep-1] - 2*data[i,j,timestep-1] +
                    data[i-1,j,timestep-1]) + (data[i,j+1,timestep-1] -
                        2*data[i,j,timestep-1] + data[i,j-1,timestep-1]
                            ))        
def InsideSet(index_array, i, j) :
    sameJ = np.where(np.isclose(index_array[:, 1], j))[0]
    sameI = np.where(np.isclose(index_array[:, 0], i))[0]
    if len(sameJ) == 0 or len(sameI) == 0 :
        return False
    possibleIs = index_array[sameJ][:,0]
    possibleJs = index_array[sameI][:,1]
    minI = min(possibleIs)
    minJ = min(possibleJs)
    maxI = max(possibleIs)
    maxJ = max(possibleJs)
    if minI <= i and i <= maxI and minJ <= j and j <= maxJ :
        return True
    return False

def diffusion_inpainting(x,y,data):
    x = x.astype(int)
    y = y.astype(int)
    
    global minY, minX
        
    minX = int(min(x))
    maxX = int(max(x)+1)
    minY = int(min(y))
    maxY = int(max(y)+1)
        
    xyZip = list(zip(x,y))
    xyArray = np.asarray(xyZip)

    saveEveryXFrames = 5
    maxNumberOfSteps = 10000

    saveCounter = 0
    intensity_values = np.zeros([maxY-minY, maxX-minX, maxNumberOfSteps])
    # Fill in the mask coordinates with their original intensity values. 
    ApplyBoundaryConditions(data, intensity_values, xyZip, 0);
    
    count = maxNumberOfSteps / 10
    
    insideIndices = []
    shell = np.copy(data)
    for i in range(shell.shape[0]):
        for j in range(shell.shape[1]):
            if InsideSet(xyArray, i, j) :
                insideIndices.append((i, j))
                    
    subs, shells = [],[]
    for t in range(1,maxNumberOfSteps):
        # Do Diffusion Step
        StepForward(intensity_values, t);
        # Re-plug in the mask intensitity values.
        ApplyBoundaryConditions(data, intensity_values, xyZip, t);
        
        #if t % count == 0 :
            #print("Current iteration = " + str(t))
        
        for idx in range(len(insideIndices)):
            i = insideIndices[idx][0]
            j = insideIndices[idx][1]
            shell[j,i] = intensity_values[j-minY, i-minX, t]

        # Remove the background (shell) from the original image.
        shells.append(np.copy(shell))
        subs.append(np.copy(data) - np.copy(shell))
        
    
    return subs, shells