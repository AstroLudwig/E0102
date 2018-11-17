#!/usr/bin/env python
"""
PURPOSE: Some handy extra tools for matplotlib colormaps
AUTHOR: Dylan Gregersen
DATE: Thu Dec 18 23:41:30 2014
"""
# ########################################################################### #

# import modules 
import numpy as np 
import pylab as plt

# ########################################################################### #

class ScalarColorMappable (plt.cm.ScalarMappable):
    """ Wrapper around matplotlib.cm.ScalarMappable 
    Set the range of the data or inferred from a set of values
            
    """       
    def __init__ (self,values=None,cmap=plt.cm.jet,vmin=None,vmax=None,clip=False): #@UndefinedVariable called plt.cm.jet
        """ 
    
        Parameters
        ----------
        cmap : `matplotlib.colormap`
        values : np.ndarray, (N,)
            values you want to map, if both vmin and vmax are given this is not necessary 
        vmin : float, None
            limits of values to map
        vmax : float,None
        
        
        Examples
        --------
        >>> values = np.linspace(0.2,5,10)
        >>> get_color = ScalarColorMappable(values,cmap=plt.cm.jet)
        >>> ax = plt.figure().add_subplot(111)
        >>> x = np.linspace(-10,10,10000)
        >>> for v in values:
        >>>    y = 2.0*x**v        
        >>>    c = get_color(v)
        >>>    ax.plot(x,y,color=c)
        >>> plt.colorbar(get_color)
        >>> plt.show()
        
        
        Notes
        -----
        __1__ the values parameter is only necessary if vmin or vmax is None
    
        """
        # get vmin and vmax
        if vmin is None:
            vmin = np.nanmin(values)
        if vmax is None:
            vmax = np.nanmax(values)
        if vmax <= vmin:
            raise ValueError("vmin must be < vmax, (vmin,vmax)={}".format((vmin,vmax)))
        # get the normalization object
        norm = plt.matplotlib.colors.Normalize(vmin,vmax,clip=clip)        
        # initialize the matplotlib scalar
        plt.cm.ScalarMappable.__init__(self,norm=norm,cmap=cmap)        
        self.set_array(np.linspace(vmin,vmax,cmap.N))        
        self.__call__ = self.to_rgba

class MonoColormap (plt.matplotlib.colors.Colormap,object):
    """
    This creates a matplotlib colormap that returns a single color
    
    Parameters
    ----------
    color : string
        The color to be used, recognizes rgb, hex, and common string names (see
        `matplotlib.colors.cnames.keys` for a list
    alpha : float, 0 to 1
        Sets the transparency of the color if not set by color
    name : None or string
        Name of the colormap, if None then named 'single_color'
    
    """
    
    def __init__ (self,color,alpha=1.0,name=None):
        # get the color and store it        
        self._color_converter = plt.matplotlib.colors.ColorConverter()        
        self._color = self._color_converter.to_rgb(color) # color tuple         
        # get some name for the colormap
        if name is None:
            name = "single_color"
        # initalize colormap
        plt.matplotlib.colors.Colormap.__init__(self,name)
        self.set_alpha(alpha)    
        self.monochrome = True     
             
    def _init (self):
        """ Initialize the lookup table with all the same color"""
        if len(self._color) != 3:
            raise ValueError("This should be 3 long, don't know why it's not")
        element = list(self._color)
        if len(element) == 3:
            element.append(self.get_alpha())
        self._lut = np.asarray([element for _ in xrange(self.N+3)])
        self._isinit = True
        self._set_extremes()    

    def get_alpha (self):
        return self._alpha
    
    def set_alpha (self,alpha):
        if alpha < 0: 
            alpha = 0.0
        elif alpha >= 1:
            alpha = 1.0
        self._alpha = alpha
    
    def get_color (self):
        """ Returns the value in self.color, the single color for this map """
        return self.color 
    
    def set_color (self,color): 
        """ Set the color of this monocromatic map
        
        See `matplotlib.colors.ColorConverter` and `to_rgb` for more info
        
        """        
        self._color = self._color_converter.to_rgb(color)
        self._init()
        
    @property
    def color (self):   
        return self._color
        
    @color.setter   
    def color (self, color):          
        self.set_color(color)
    
    @color.getter
    def color (self):
        return self._color

def truncate_cmap (cmap,n_min=0,n_max=256):
    """ Generate a truncated colormap 
    Colormaps have a certain number of colors (usually 255) which is 
    available via cmap.N
    This is a simple function which returns a new color map from a 
    subset of the input colormap. For example, if you only want the jet 
    colormap up to green colors you may use 
    tcmap = truncate_cmap(plt.cm.jet,n_max=150)
    This function is especially useful when you want to remove the white
    bounds of a colormap 
    Parameters 
    cmap : plt.matplotlib.colors.Colormap 
    n_min : int 
    n_max : int 
    Return 
    truncated_cmap : plt.matplotlib.colors.Colormap
        
    """
    color_index = np.arange(n_min,n_max).astype(int)
    colors = cmap(color_index)
    name = "truncated_{}".format(cmap.name)
    return plt.matplotlib.colors.ListedColormap(colors,name=name)

# ########################################################################### #

# I have the following added to my PYTHONSTARTUP file so that I can use 
# them interactively when the need arises. Particularly `show_mpl_cmaps()`

def mpl_grayify_cmap(cmap):
    """Return a grayscale version of the colormap
    Source: https://jakevdp.github.io/blog/2014/10/16/how-bad-is-your-colormap/
    Thank's Jake!
    """
    cmap = plt.cm.get_cmap(cmap)
    colors = cmap(np.arange(cmap.N))
    
    # convert RGBA to perceived greyscale luminance
    # cf. http://alienryderflex.com/hsp.html
    RGB_weight = [0.299, 0.587, 0.114]
    luminance = np.sqrt(np.dot(colors[:, :3] ** 2, RGB_weight))
    colors[:, :3] = luminance[:, np.newaxis]
    
    return cmap.from_list(cmap.name + "_grayscale", colors, cmap.N)

def mpl_get_colormaps ():
    """ Return a dictionary of lists of colormap names
    The dictionary keys refer to types of colormaps (sequential,diverging,etc)
    And the lists are matplotlib colormaps names you can get to via plt.cm.???
    Returns
    categorized_cmaps : dict of lists of strings
    Inspired from: http://matplotlib.org/examples/color/colormaps_reference.html
    """
    from collections import OrderedDict
    cmaps = OrderedDict(\
        [('Sequential',
            ['Blues', 'BuGn', 'BuPu',
            'GnBu', 'Greens', 'Greys', 'Oranges', 'OrRd',
            'PuBu', 'PuBuGn', 'PuRd', 'Purples', 'RdPu',
            'Reds', 'YlGn', 'YlGnBu', 'YlOrBr', 'YlOrRd']
            ),
        ('Sequential (2)', 
            ['cool', 'afmhot', 'autumn', 'bone', 'copper',
            'gist_heat', 'gray', 'hot', 'pink',
            'spring', 'summer', 'winter']
            ),
        ('Diverging',      
            ['BrBG', 'coolwarm', 'bwr',  'PiYG', 'PRGn', 'PuOr',
            'RdBu', 'RdGy', 'RdYlBu', 'RdYlGn', 'Spectral',
            'seismic']
            ),
        ('Qualitative',    
            ['Accent', 'Dark2', 'Paired', 'Pastel1',
            'Pastel2', 'Set1', 'Set2', 'Set3']
            ),
        ('Miscellaneous',  
            ['gist_earth', 'terrain', 'ocean', 'gist_stern',
            'brg', 'CMRmap', 'cubehelix',
            'gnuplot', 'gnuplot2', 'gist_ncar',
            'nipy_spectral', 'jet', 'rainbow',
            'gist_rainbow', 'hsv', 'flag', 'prism']
            )]
        )
    return cmaps 

def show_mpl_cmaps (with_gray=True):
    """ Create and show a plot of the matplotlib colormaps 
    Parameters
    with_gray : bool 
        If True then the grayified version of the cmap is also plotted
    Inspired from: http://matplotlib.org/examples/color/colormaps_reference.html
    """
    import numpy as np 
    import pylab as plt

    gradient = np.linspace(0, 1, 256)
    gradient = np.vstack((gradient, gradient))

    cmaps = mpl_get_colormaps()
    naxes = sum(map(len,cmaps.values()))
    
    if with_gray:
        ncols = 2
    else:
        ncols = 1
    fig,axes = plt.subplots(naxes,ncols,figsize=(5,9),dpi=100)
    fig.subplots_adjust(top=0.99, bottom=0.01, left=0.2, right=0.99)
    if not with_gray:
        axes = axes.reshape((-1,1))

    i = 0
    for category in cmaps:
        for j,name in enumerate(cmaps[category]):
            axes[i][0].imshow(gradient, aspect='auto', cmap=plt.get_cmap(name))
            axes[i][0].set_axis_off()
            pos = list(axes[i][0].get_position().bounds)
            x_text = pos[0] - 0.01
            y_text = pos[1] + pos[3]/2.
            fig.text(x_text, y_text, name, va='center', ha='right', fontsize=9)
            if j == 0:
                x_text = pos[0] + pos[2]+0.01
                y_text = pos[1] + pos[3]/2.
                fig.text(x_text,y_text,category,va='top',ha='left',rotation="vertical",fontsize=10)

            if with_gray:
                axes[i][1].imshow(gradient, aspect='auto', cmap=mpl_grayify_cmap(plt.get_cmap(name)))        
                axes[i][1].set_axis_off()
                pos = list(axes[i][0].get_position().bounds)
            i += 1            
plt.show() 