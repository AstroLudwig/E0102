import numpy as np 

def extend_mask(extend,edge_x,edge_y,mask_x,mask_y):
    mid_x = (np.max(edge_x)+np.min(edge_x))/2; mid_y = (np.max(edge_y)+np.min(edge_y))/2
    mask_x = np.asarray(mask_x);     mask_y = np.asarray(mask_y)
    push_x = np.copy(mask_x);        push_y = np.copy(mask_y)

    
    for i in range(len(mask_x)): 
        if mask_x[i] < mid_x:
            push_x[i] = mask_x[i] - extend
        if mask_x[i] > mid_x:
            push_x[i] = mask_x[i] + extend
        if mask_y[i] > mid_y:
            push_y[i] = mask_y[i] + extend
        if mask_y[i] < mid_y:
            push_y[i] = mask_y[i] - extend

    # ----> Top and Bottom 
    fpx = mid_x - 1 - extend
    lpx = mid_x + 1 + extend 
    fpy = np.max(mask_y) + extend
    fppy = np.min(mask_y) - extend
    fyrr = np.repeat(fpy,lpx-fpx)
    fypr = np.repeat(fppy,lpx-fpx)           
    fxrr = np.arange(fpx,lpx)
    push_x = np.append(push_x,fxrr)
    push_y = np.append(push_y,fyrr)
    push_x = np.append(push_x,fxrr)
    push_y = np.append(push_y,fypr)
    # -----> Left and Right
    kpx = np.max(mask_x) + extend 
    kppx = np.min(mask_x) - extend
    kpy = mid_y - 1 - extend
    kppy = mid_y + 1 + extend 
    kxrr = np.repeat(kpx,kppy - kpy)
    kxpr = np.repeat(kppx,kppy - kpy)
    kyrr = np.arange(kpy,kppy)
    push_x = np.append(push_x,kxrr)
    push_y = np.append(push_y,kyrr)
    push_x = np.append(push_x,kxpr)
    push_y = np.append(push_y,kyrr)

    return push_x,push_y

# Reduce Mask 
def mask_reduction(x,y):
    # Create mask arrays. 
    maskx_x, maskx_y = [],[]
    for i in range(np.max(x)): 
        # remove duplicates in X          
        look = np.where(x == i)[0]
        dim = np.shape(look)
        if dim[0] > 0: 
             maskx_y.append(y[look[0]])
             maskx_y.append(y[look[len(look)-1]])
             maskx_x.append(x[look[0]])
             maskx_x.append(x[look[len(look)-1]])

    masky_x, masky_y = [],[]
    for i in range(np.max(y)):  
        # remove duplicates in Y
        look = np.where(y == i)[0]
        dim = np.shape(look)
        if dim[0] > 0: 
             masky_y.append(y[look[0]])
             masky_y.append(y[look[len(look)-1]])
             masky_x.append(x[look[0]])
             masky_x.append(x[look[len(look)-1]])

    mask_x = maskx_x + masky_x         
    mask_y = maskx_y + masky_y

    return np.asarray(list(zip(mask_x,mask_y)))