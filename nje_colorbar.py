###############################################################################
# Author: Chris Slocum
# Name: colorbar.py
# Version: 1.0

# Description:
# Functions to creatae and modify Python colormaps
###############################################################################

def cmap_discretize(cmap, N):

  import numpy as np
  import matplotlib.colors as colors
  import matplotlib.pyplot as plt

  """
	Return a discrete colormap from the continuous colormap cmap.
	Requires numpy and matplotlib.colors module

	cmap: colormap instance, e.g., jet
	N : number of colors returned

	Example:
	discrete_jet = cmap_discretize('jet', 5)
  """

  if type(cmap) == str:
	cmap = plt.cm.get_cmap(cmap)
  colors_i = np.concatenate((np.linspace(0, 1., N), (0.,0.,0.,0.)))
  colors_rgba = cmap(colors_i)
  indices = np.linspace(0, 1., N+1)
  cdict = {}
  for ki,key in enumerate(('red','green','blue')):
	cdict[key] = [ (indices[i], colors_rgba[i-1,ki], \
			colors_rgba[i,ki]) for i in xrange(N+1) ]
  # Return colormap object.
  return colors.LinearSegmentedColormap(cmap.name + "_%d"%N, cdict, 1024)

#NAME
#    Custom Colormaps for Matplotlib
#PURPOSE
#    This program shows how to implement make_cmap which is a function that
#    generates a colorbar.  If you want to look at different color schemes,
#    check out https://kuler.adobe.com/create.
#PROGRAMMER(S)
#    Chris Slocum
#REVISION HISTORY
#    20130411 -- Initial version created
#    20140313 -- Small changes made and code posted online
#    20140320 -- Added the ability to set the position of each color
#
#    make_cmap takes a list of tuples which contain RGB values. The RGB
#    values may either be in 8-bit [0 to 255] (in which bit must be set to
#    True when called) or arithmetic [0 to 1] (default). make_cmap returns
#    a cmap with equally spaced colors.
#    Arrange your tuples so that the first color is the lowest value for the
#    colorbar and the last is the highest.
#    position contains values from 0 to 1 to dictate the location of each color.

#  EXAMPLE
#  Create a list of RGB tuples
#colors = [(0,0,0), (0,0,0),(102,0,204),(0,0,255), (51,153,255), \
#	(0,255,0), (255,255,255), (255,255,0), (255,130,0), \
#	(255,0,0), (255,0,0), (0,0,0), (255,255,255)]  # This example uses the 8-bit RGB
#  Call the function make_cmap which returns your colormap
#my_cmap = make_cmap(colors, bit=True)

def make_cmap(colors, position=None, bit=False, max=256):

    import matplotlib as mpl
    import numpy as np
    
    bit_rgb = np.linspace(0,1,max)
    if position == None:
        position = np.linspace(0,1,len(colors))
    else:
        if len(position) != len(colors):
            sys.exit("position length must be the same as colors")
        elif position[0] != 0 or position[-1] != 1:
            sys.exit("position must start with 0 and end with 1")
    if bit:
        for i in range(len(colors)):
            colors[i] = (bit_rgb[colors[i][0]],
                         bit_rgb[colors[i][1]],
                         bit_rgb[colors[i][2]])
    cdict = {'red':[], 'green':[], 'blue':[]}
    for pos, color in zip(position, colors):
        cdict['red'].append((pos, color[0], color[0]))
        cdict['green'].append((pos, color[1], color[1]))
        cdict['blue'].append((pos, color[2], color[2]))

    cmap = mpl.colors.LinearSegmentedColormap('my_colormap',cdict,256)
    return cmap





