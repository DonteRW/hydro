'''
###############################################################################
Author: Nicholas Elmer (UAH/SPoRT)
Date: July 2016
Name: nje_WRF-Hydro_QCtool_grid.py
Version: 1.0

Description:
Quality control CHANNELGRID.nc and LAKEGRID.nc files from 
    WRF-Hydro ArcGIS preprocessing tool for grid-based routing.
###############################################################################
'''

import numpy as np
from mpl_toolkits.basemap import Basemap as bmap
import os.path as os
from netCDF4 import Dataset as netcdf

from nje_read_netcdf import getdata
from nje_colorbar import make_cmap, cmap_discretize 

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.figure as mplfig
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg

import sys
if sys.version_info[0] < 3:
    from Tkinter import *
else:
    from tkinter import *
    
def netcdf_write(srcfileobj, dstfile, newVarGrid):
    print dstfile
    # open for writing
    dstfileobj = netcdf(dstfile, 'w', format=srcfileobj.file_format)
    # copy dimensions
    for dname, the_dim in srcfileobj.dimensions.iteritems():
        dstfileobj.createDimension(dname, \
                len(the_dim) if not the_dim.isunlimited() else None)
    # copy attributes
    dstfileobj.setncatts({k:srcfileobj.getncattr(k) \
            for k in srcfileobj.ncattrs()})
    # copy variables
    for v_name, varin in srcfileobj.variables.iteritems():
        outVar = dstfileobj.createVariable(v_name, varin.datatype, varin.dimensions)
        outVar.setncatts({k:varin.getncattr(k) for k in varin.ncattrs()})
        if v_name == 'CHANNELGRID' or v_name == 'LAKEGRID':
            outVar[:] = newVarGrid
        else:
            outVar[:] = varin[:]
    print dstfileobj
    dstfileobj.close()
    return

def getbounds(ind, xelem, yelem):
    # create mesh
    xshape = np.shape(xelem)[0]
    yshape = np.shape(yelem)[0]
    xgrid, ygrid = np.meshgrid(xelem, yelem)
    # get bounds and create 2 pixel buffer
    minlatind, maxlatind = np.min(ind[1])-2, np.max(ind[1])+2
    minlonind, maxlonind = np.min(ind[0])-2, np.max(ind[0])+2
    # make sure indices are within array bounds
    if minlonind <= 0: minlonind = 0
    if minlatind <= 0: minlatind = 0
    if maxlonind >= yshape: maxlonind = yshape - 1
    if maxlatind >= xshape: maxlatind = xshape - 1
    # retrieve indices of gridbox
    try:
      Xsubset = xgrid[minlonind:maxlonind+1, minlatind:maxlatind+1]
      Ysubset = ygrid[minlonind:maxlonind+1, minlatind:maxlatind+1]
    except IndexError:
      Xsubset = xgrid[minlonind:, minlatind:]
      Ysubset = ygrid[minlonind:, minlatind:]

    subshape = np.shape(Xsubset)
    bound = [Xsubset[0,0], Xsubset[-1,-1],  \
                Ysubset[0,0], Ysubset[-1,-1] ]
    center = [Xsubset[subshape[0]/2,subshape[1]/2], \
                 Ysubset[subshape[0]/2,subshape[1]/2]]
    ind = [minlonind, maxlonind, minlatind, maxlatind]
    return (bound, center, ind)
        
def haversine(lon, lat, lon_2Darray, lat_2Darray):
    """
    Calculate the great circle distance between two points 
    on the earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians
    lon_2Darray = np.deg2rad(lon_2Darray)
    lat_2Darray = np.deg2rad(lat_2Darray)
    latlon = np.radians([lon, lat])
    # haversine formula 
    dlon = lon_2Darray - latlon[0] 
    dlat = lat_2Darray - latlon[1] 
    a = np.sin(dlat/2)**2 + np.cos(latlon[1]) * np.cos(lat_2Darray) * np.sin(dlon/2)**2
    arcdist = 2 * np.arcsin(np.sqrt(a)) 
    minind = np.where(np.abs(arcdist - np.min(arcdist) == 0))
    return minind

class Tool(object):

  def _save(self):
    print 'Saving...'
    try:
      # specify output filenames
      chan_dst = os.join(self.filepath, 'CHANNELGRID_QC_currentLake%i.nc' %(self.lakenum[self.thisLakeInd-1]))
      lake_dst = os.join(self.filepath,    'LAKEGRID_QC_currentLake%i.nc' %(self.lakenum[self.thisLakeInd-1]))
    except IndexError:
      # specify output filenames
      chan_dst = os.join(self.filepath, 'CHANNELGRID_QC_currentLake0.nc')
      lake_dst = os.join(self.filepath,    'LAKEGRID_QC_currentLake0.nc')
      
    try:
      # save modified channelgrid and lakegrid files
      Cwrite = netcdf_write(self.ncfileCHANr, chan_dst, self.channelgrid)
      Lwrite = netcdf_write(self.ncfileLAKEr, lake_dst, self.lakecopy)
    except TypeError as typerr:
        print typerr
        self._quit(save=False)
    return

  def _quit(self, save=False):
    if save: self._save()
    try:
      if self.ncfileCHANr._isopen: self.ncfileCHANr.close()
      if self.ncfileLAKEr._isopen: self.ncfileLAKEr.close()
    except (RuntimeError, AttributeError) as runerr:
      print runerr
    finally:
      self.fig.canvas.mpl_disconnect(self.cid) 
      self.fig.canvas.mpl_disconnect(self.kid)
      self.master.quit()     # stops mainloop
      self.master.destroy()  # this is necessary on Windows to prevent
                             # Fatal Python Error: PyEval_RestoreThread: NULL tstate
    print 'Quitting...'
    return
    
  def plot_subset(self, whichLake=1, direction='', zoom=''):
    # whichLake: -1 for previous, 0 for current, 1 for next
    # go back one
    if whichLake == -1:
        if self.thisLakeInd > 1:
            self.thisLakeInd -= 2
        else:
            self.thisLakeInd = len(self.lakenum) - 1
        newLake = True
    # stay at current
    elif whichLake == 0:
        self.thisLakeInd -= 1
        newLake = False
    # go to next lake
    else:
        newLake = True
        if self.thisLakeInd >= len(self.lakenum):
            self.thisLakeInd = 0
    # get lake indices from grid
    lakeind = np.where(self.lakegrid == self.lakenum[self.thisLakeInd])
    # get lake bounds
    lakbound, lakcenter, self.lakinds_full = \
            getbounds(lakeind, self.xdim, self.ydim)

    # determine which lakinds to use
    if newLake:
        lakinds_which = self.lakinds_full
    else:
        lakinds_which = self.lakinds_current

    # zoom if needed
    if zoom == 'in':
        # zoom in 10% each side
        xzoom=int((lakinds_which[1] - lakinds_which[0])*0.1)
        yzoom=int((lakinds_which[3] - lakinds_which[2])*0.1)
        if xzoom == 0 or yzoom == 0:
            self.lakinds_zoom = lakinds_which
            print '...maximum zoom reached...'
        else:
            self.lakinds_zoom = [lakinds_which[0]+xzoom, lakinds_which[1]-xzoom, \
                             lakinds_which[2]+yzoom, lakinds_which[3]-yzoom]
        # TODO need to modify lakbound too for mapping
    elif zoom == 'out':
        # zoom out 10% each side
        xzoom=int((lakinds_which[1] - lakinds_which[0])*0.1)
        yzoom=int((lakinds_which[3] - lakinds_which[2])*0.1)
        if xzoom == 0 or yzoom == 0:
            self.lakinds_zoom = [lakinds_which[0]-2, lakinds_which[1]+2, \
                             lakinds_which[2]-2, lakinds_which[3]+2]
        elif xzoom >= self.shape[0]/100 or yzoom >= self.shape[1]/100:
            xzoom, yzoom = 0, 0
            print '...minimum zoom reached...'
        else:
            self.lakinds_zoom = [lakinds_which[0]-xzoom, lakinds_which[1]+xzoom, \
                             lakinds_which[2]-yzoom, lakinds_which[3]+yzoom]
    else:
        self.lakinds_zoom = lakinds_which
    
    xind_diff = self.lakinds_zoom[3] - self.lakinds_zoom[2]
    yind_diff = self.lakinds_zoom[1] - self.lakinds_zoom[0]
    
    # pan left/right/up/down by quarter frame if requested
    if   direction == 'up':
        self.lakinds = [self.lakinds_zoom[0]-yind_diff/4, self.lakinds_zoom[1]-yind_diff/4, \
                        self.lakinds_zoom[2], self.lakinds_zoom[3]]
    elif direction == 'down':
        self.lakinds = [self.lakinds_zoom[0]+yind_diff/4, self.lakinds_zoom[1]+yind_diff/4, \
                        self.lakinds_zoom[2], self.lakinds_zoom[3]]
    elif direction == 'left':
        self.lakinds = [self.lakinds_zoom[0], self.lakinds_zoom[1], \
                        self.lakinds_zoom[2]-xind_diff/4, self.lakinds_zoom[3]-xind_diff/4]
    elif direction == 'right':
        self.lakinds = [self.lakinds_zoom[0], self.lakinds_zoom[1], \
                        self.lakinds_zoom[2]+xind_diff/4, self.lakinds_zoom[3]+xind_diff/4]
    else:
        self.lakinds = self.lakinds_zoom
    
    # make sure new lakinds are within array bounds
    if self.lakinds[0] <= 0:
        tempdiff = 0 - self.lakinds[0]
        self.lakinds[0] = 0
        self.lakinds[1] = yind_diff
    if self.lakinds[2] <= 0:
        tempdiff = 0 - self.lakinds[2]
        self.lakinds[2] = 0
        self.lakinds[3] = xind_diff
    if self.lakinds[1] >= self.shape[0]:
        tempdiff = self.lakinds[1] - self.shape[0]
        self.lakinds[0] = int(self.shape[0])-(yind_diff+1)
        self.lakinds[1] = self.shape[0]-1
    if self.lakinds[3] >= self.shape[1]:
        tempdiff = self.lakinds[3] - self.shape[1]
        self.lakinds[2] = int(self.shape[1])-(xind_diff+1)
        self.lakinds[3] = self.shape[1]-1

    # save current lakinds
    self.lakinds_current = self.lakinds
    # get lakbound and lakcenter in lat/lon coordinates too
    cLon, cLat         = self.map1( lakcenter[0], lakcenter[1], inverse=True)
    edgeLon, edgeLat   = self.map1(lakbound[0:2], lakbound[2:], inverse=True)
    
    # plot lake location on full map
    self.lakeloc = self.map1.scatter([lakcenter[0]], [lakcenter[1]], \
                marker='*', color='r', edgecolor='r', s=200.0, \
                zorder=10, ax=self.ax1)
    self.ax2.set_title('Lake %i' %(self.lakenum[self.thisLakeInd]))
    submerge = self.merge[self.lakinds[0]:self.lakinds[1]+1, \
            self.lakinds[2]:self.lakinds[3]+1]
    
    # draw inset
    if self.map_subset: # not yet supported
        self.map2.llcrnrx = lakbound[0]
        self.map2.llcrnry = lakbound[3]
        self.map2.urcrnrx = lakbound[1]
        self.map2.urcrnry = lakbound[2]
        print np.shape(self.xdim), np.shape(self.ydim), np.shape(self.channelgrid)
        goodLakeInd = np.where( \
                (self.lakecenter[:,0] >= edgeLon[0]) & \
                (self.lakecenter[:,0] <= edgeLon[1]) & \
                (self.lakecenter[:,1] >= edgeLat[1]) & \
                (self.lakecenter[:,1] <= edgeLat[0]) )
        print goodLakeInd
        for feature in goodLakeInd[0]:
            goodLake = self.polygons[feature]
            xfeat, yfeat = self.map2(*self.map1(goodLake[:,0], goodLake[:,1], inverse=True))
            self.map2.plot(xfeat, yfeat, ax=self.ax2, color='b', zorder=20)

        self.zoomgrid = self.map2.imshow(submerge, origin='upper', \
            vmin=self.subplot_vmin, vmax=self.subplot_vmax, \
            interpolation='none', cmap=self.cmap, ax=self.ax2, zorder=1)
    else:
        self.zoomgrid = self.ax2.imshow(submerge, origin='upper', \
                    vmin=self.subplot_vmin, vmax=self.subplot_vmax, \
                    interpolation='none', cmap=self.cmap, zorder=1)
                    
    self.canvas.draw()
    self.thisLakeInd += 1
    return

  def key(self, event):
    if event.key == 'n' or event.key == ' ': # go to next lake
        if self.thisLakeInd >= 0 and self.thisLakeInd < len(self.lakenum):
            print 'Advancing to next lake...Lake %i' \
                                %(self.lakenum[self.thisLakeInd])
        else:
            print 'Wrapping around to first lake...Lake %i' %(self.lakenum[0])
        self.ax2.clear()
        self.lakeloc.remove()
        self.lakinds_current = self.lakinds_reset
        self.plot_subset(whichLake=1)
    elif event.key == 'q' or event.key == 'escape': # quit
        self._quit(save=True)
    elif event.key == 's' or event.key == 'w': # save
        self._save()
    elif event.key == 'u' or event.key == 'up': # pan up
        self.ax2.clear()
        self.lakeloc.remove()
        self.plot_subset(whichLake=0, direction='up')
    elif event.key == 'd' or event.key == 'down': # pan down
        self.ax2.clear()
        self.lakeloc.remove()
        self.plot_subset(whichLake=0, direction='down')
    elif event.key == 'l' or event.key == 'left': # pan left
        self.ax2.clear()
        self.lakeloc.remove()
        self.plot_subset(whichLake=0, direction='left')
    elif event.key == 'r' or event.key == 'right': # pan right
        self.ax2.clear()
        self.lakeloc.remove()
        self.plot_subset(whichLake=0, direction='right')
    elif event.key == 'z' or event.key == 'i': #zoom in
        self.ax2.clear()
        self.lakeloc.remove()
        self.plot_subset(whichLake=0, zoom='in')
    elif event.key == 'f': #full extent
        self.ax2.clear()
        self.lakeloc.remove()
        self.lakinds_current = self.lakinds_full
        self.plot_subset(whichLake=0)
    elif event.key == 'y' or event.key == 'o': #zoom out
        self.ax2.clear()
        self.lakeloc.remove()
        self.plot_subset(whichLake=0, zoom='out')
    elif event.key == 'g' or event.key == 'j': #go/jump to lake
        # ask for user number input
        gotoLake = raw_input('Enter Lake Number: ')
        gotoLake = int(gotoLake)
        if gotoLake >= 0 and gotoLake <= self.lakenum[-1]:
            self.thisLakeInd = np.argmin(np.abs(np.array(self.lakenum)-gotoLake))
            print 'Moving to Lake %i...' %(self.lakenum[self.thisLakeInd])
        elif gotoLake > self.lakenum[-1]:
            self.thisLakeInd = len(self.lakenum)-1
            print 'Moving to Lake %i...' %(self.lakenum[-1])
        else:
            self.thisLakeInd = 0
            print 'Moving to Lake 1...'          
        self.ax2.clear()
        self.lakeloc.remove()
        self.plot_subset(whichLake=1)
    elif event.key == 'backspace' or event.key == 'b': # go back one lake
        print self.thisLakeInd
        if self.thisLakeInd > 1 and self.thisLakeInd <= len(self.lakenum):
            print 'Returning to prev lake...Lake %i' \
                                %(self.lakenum[self.thisLakeInd - 2])
        else:
            print 'Wrapping around to last lake...Lake %i' %(self.lakenum[-1])
        self.ax2.clear()
        self.lakeloc.remove()
        self.lakinds_current = self.lakinds_reset
        self.plot_subset(whichLake=-1)
    else: print 'Pressing the %s key does nothing...' %repr(event.key)
    return

  def click(self, event):
    templakinds = self.lakinds
    if event.button == 1:
      try:
        if self.map_subset: # not currently supported
            cell_lon, cell_lat = self.map2(event.xdata, event.ydata, inverse=True)
            nearest = haversine(cell_lon, cell_lat, longitude, latitude)
            xcell, ycell = nearest[0], nearest[1]
        else:
            # event.xdata matches with ycell for some reason, and vice-versa
            xLcell, yLcell = int(np.round(event.xdata)), int(np.round(event.ydata))
            xcell, ycell = int(self.lakinds[0]+yLcell), int(self.lakinds[2]+xLcell)
        # orig back to new lake
        if (self.merge[xcell, ycell] == 1.0) & \
                    (self.lakenum[self.thisLakeInd - 1] in self.channelgrid):
            self.merge[xcell, ycell]=6.0     # green
            self.lakecopy[xcell, ycell] = self.lakenum[self.thisLakeInd - 1]
        # orig chan to new back
        elif (self.merge[xcell, ycell] == 2.0):
            self.merge[xcell, ycell]=4.0     # pink
            self.channelgrid[xcell, ycell] = self.fillval
        # orig lake to new back
        elif (self.merge[xcell, ycell] == 3.0) & \
             (self.lakegrid[xcell, ycell] == self.lakenum[self.thisLakeInd - 1]):
            self.merge[xcell, ycell]=5.0     # white
            self.lakecopy[xcell, ycell] = self.fillval
        elif (self.merge[xcell, ycell] == 3.0) & \
             (self.lakegrid[xcell, ycell] != self.lakenum[self.thisLakeInd - 1]):
            print '  Not current lake...'
        # from orig chan / new back to new lake
        elif (self.merge[xcell, ycell] == 4.0):
            self.merge[xcell, ycell]=7.0     # purple
            self.lakecopy[xcell, ycell] = self.lakenum[self.thisLakeInd - 1]
        # from orig lake / new back to orig lake
        elif (self.merge[xcell, ycell] == 5.0):
            if (not self.lakenum[self.thisLakeInd - 1] in self.channelgrid):
                ghostLakeInd = np.where(self.lakegrid==self.lakenum[self.thisLakeInd-1])
                self.merge[ghostLakeInd]=3.0     # aqua
                self.lakecopy[ghostLakeInd] = self.lakenum[self.thisLakeInd-1]
                ghostNodeInd = np.where(self.chanorig==self.lakenum[self.thisLakeInd-1])
                self.merge[ghostNodeInd]=0.0     # yellow
                self.channelgrid[ghostNodeInd] = self.lakenum[self.thisLakeInd-1]
                self.lakecopy[ghostNodeInd] = self.lakenum[self.thisLakeInd-1]
            else:
                self.merge[xcell, ycell]=3.0     # aqua
                self.lakecopy[xcell, ycell] = self.lakenum[self.thisLakeInd - 1]
        # from orig back / new lake to new back or new chan
        elif (self.merge[xcell, ycell] == 6.0):
            if (self.strorder[xcell, ycell] != -9999):
                self.merge[xcell, ycell]=8.0     # orange
                self.lakecopy[xcell, ycell] = self.fillval
                self.channelgrid[xcell,ycell] = 0.0
            else:
                self.merge[xcell, ycell]=1.0     # gray
                self.lakecopy[xcell, ycell] = self.fillval
        # from orig chan / new back / new lake to orig chan
        elif (self.merge[xcell, ycell] == 7.0):
            self.merge[xcell, ycell]=2.0     # red
            self.channelgrid[xcell, ycell] = 0.0
            self.lakecopy[xcell, ycell] = self.fillval
        # delete lake
        elif (self.merge[xcell, ycell] == 0.0) & \
                (self.lakecopy[xcell,ycell] == self.lakenum[self.thisLakeInd-1]):
            newLakeInd = np.where(self.lakecopy==self.lakenum[self.thisLakeInd-1])
            self.merge[newLakeInd]=1.0
            self.lakecopy[newLakeInd] = self.fillval 
            self.channelgrid[newLakeInd] = self.chanorig[newLakeInd] #self.fillval
            ghostLakeInd = np.where(self.lakegrid==self.lakenum[self.thisLakeInd-1])
            self.merge[ghostLakeInd]=5.0     # aqua
            self.lakecopy[ghostLakeInd] = self.fillval         
            self.merge[xcell, ycell] =-1.0     # black
            self.channelgrid[xcell, ycell] = self.fillval
            self.lakecopy[xcell, ycell] = self.fillval
        # restore lake
        elif (self.merge[xcell, ycell] == -1.0):
            ghostLakeInd = np.where(self.lakegrid==self.lakenum[self.thisLakeInd-1])
            self.merge[ghostLakeInd]=3.0     # aqua
            self.lakecopy[ghostLakeInd] = self.lakenum[self.thisLakeInd-1]
            self.channelgrid[ghostLakeInd] = self.chanorig[ghostLakeInd]
            self.merge[xcell, ycell]=0.0     # yellow
            self.channelgrid[xcell, ycell] = self.lakenum[self.thisLakeInd-1]
            self.lakecopy[xcell, ycell] = self.lakenum[self.thisLakeInd-1]
        # create new node
        elif (self.merge[xcell, ycell] >= -1.0) & \
                    (not self.lakenum[self.thisLakeInd - 1] in self.channelgrid):
            print 'Cannot create new nodes...'
        # from new chan to orig back
        elif (self.merge[xcell, ycell] == 8.0):
                self.merge[xcell, ycell]=1.0     # gray
                self.lakecopy[xcell, ycell] = self.fillval
                self.channelgrid[xcell, ycell] = self.fillval
        else:
            pass
        self.ax2.clear()
        self.lakeloc.remove()
        self.plot_subset(whichLake=0)
      except:
        print 'Error. Try Again... '
        self.lakinds = templakinds
    else:
        print 'Left click only!'
    return
        
  def __init__(self, master, gridpath, lakefile, map_subset=False, \
                  open_saved_grids=False, currentLake=0):
                      
    # define attributes
    self.master = master
    self.filepath = gridpath
    self.lake_shapefile = os.splitext(lakefile)[0]
    self.map_subset = map_subset
    if self.map_subset:
        print 'Mapping subset not yet supported.'
        self.map_subset = False
    
    # set up interactive figure
    self.fig = mplfig.Figure(figsize=(20,15))
    self.ax2 = self.fig.add_subplot(121)
    self.ax1 = self.fig.add_subplot(122)  
    self.canvas = FigureCanvasTkAgg(self.fig, master=self.master)
    
    # set up event connections
    self.cid = self.canvas.mpl_connect('button_press_event', self.click)  
    self.kid = self.canvas.mpl_connect('key_press_event', self.key)
    quit_but = Button(master=self.master, text='Quit', command=self._quit)
    quit_but.pack(side=BOTTOM)
    self.canvas.get_tk_widget().pack()
    
    # create color map for plotting
    colors = [(0,0,0), (255,255,0), (180,180,180), (255,0,0), (0,255,255), \
        (255,200,200), (255,255,255), (0,255,0), (200,200,255), (255,120,0)]
    pltcmap = make_cmap(colors, bit=True)
    self.cmap = pltcmap
    self.subplot_vmin =-1.
    self.subplot_vmax = 8.
    
    # read in channelgrid.nc, lakegrid.nc, lakes.shp
    if open_saved_grids:
        chan_src = os.join(self.filepath, \
                'CHANNELGRID_QC_currentLake%i.nc' %(currentLake))
        lake_src = os.join(self.filepath, \
                'LAKEGRID_QC_currentLake%i.nc' %(currentLake))
    else:
        chan_src = os.join(self.filepath, 'CHANNELGRID.nc')
        lake_src = os.join(self.filepath, 'LAKEGRID.nc')
    strord_src = os.join(self.filepath, 'str_order.nc')

    self.ncfileCHANr = netcdf(chan_src,'r')
    self.channelgrid = self.ncfileCHANr.variables['CHANNELGRID'][:]
    self.chanrange   = [np.min(self.channelgrid), np.max(self.channelgrid)]
    
    self.ncfileLAKEr = netcdf(lake_src,'r')
    self.lakegrid    = self.ncfileLAKEr.variables['LAKEGRID'][:]
    self.lakerange   = [np.min(self.lakegrid), np.max(self.lakegrid)]
    
    self.ncfileORDER = netcdf(strord_src,'r')
    self.strorder    = self.ncfileORDER.variables['STREAMORDER'][:]
    self.ncfileORDER.close()
    
    self.latitude    = getdata(os.join(self.filepath,'latitude.nc'), ['LATITUDE'])['LATITUDE']['data']
    self.latrange    = [np.min(self.latitude), np.max(self.latitude)]
    self.longitude   = getdata(os.join(self.filepath,'longitude.nc'), ['LONGITUDE'])['LONGITUDE']['data']
    self.lonrange    = [np.min(self.longitude), np.max(self.longitude)]
    self.fillval = np.min(self.channelgrid)
    
    # make copy of lakegrid and channelgrid for plotting purposes
    self.lakecopy = np.copy(self.lakegrid)
    self.chanorig = np.copy(self.channelgrid)
    
    # get shape and extent of full domain
    self.shape = np.shape(self.channelgrid)
    self.lakinds_current = [0, self.shape[1], 0, self.shape[0]]
    self.lakinds_reset = [0, self.shape[1], 0, self.shape[0]]
    
    # get projection information
    self.mapproj = self.ncfileCHANr.variables['lambert_conformal_conic']
    self.lon_0 = self.mapproj.longitude_of_central_meridian
    self.lat_0 = self.mapproj.latitude_of_projection_origin
    self.lat_1 = self.mapproj.standard_parallel[0]
    self.lat_2 = self.mapproj.standard_parallel[1]
    xdim = self.ncfileCHANr.variables['x'][:]
    ydim = self.ncfileCHANr.variables['y'][:]
    # shift xdim and ydim to match basemap convention
    self.xdim = xdim - np.min(xdim)
    self.ydim = ydim - np.min(ydim)
    xdim = ydim = None
    
    # merge channelgrid and lake grid for displaying subset
    self.merge = np.ones_like(self.channelgrid, dtype=float)
    self.lakind  = np.where(self.lakegrid    > 0)
    self.nodeind = np.where(self.channelgrid > 0)
    self.chanind = np.where(self.channelgrid == 0)
    # background equals 1.0
    self.merge[self.chanind]  = 2.0	#streams - red
    self.merge[self.lakind]   = 3.0	#lakes   - aqua
    self.merge[self.nodeind]  = 0.0   #nodes   - yellow
    
    # create full domain map
    parallels = np.arange(  20.,  55., 1.)
    meridians = np.arange(-130., -60., 1.)
    self.map1 = bmap(width=np.max(self.xdim)-np.min(self.xdim), \
                height=np.max(self.ydim)-np.min(self.ydim), \
                projection='lcc', resolution='l', \
                lat_1=self.lat_1, lat_2=self.lat_2, \
                lat_0=self.lat_0, lon_0=self.lon_0, \
                )
    if self.map_subset: # not yet supported
        self.map2 = bmap(width=np.max(self.xdim)-np.min(self.xdim), \
                height=np.max(self.ydim)-np.min(self.ydim), \
                projection='lcc', resolution='l', \
                lat_1=self.lat_1, lat_2=self.lat_2, \
                lat_0=self.lat_0, lon_0=self.lon_0, \
                )

    # read in lake polygons from shapefile
    self.lakepoly = self.map1.readshapefile(self.lake_shapefile, \
                    'lakes', color='b', linewidth=1.0, ax=self.ax1)
    self.map1.drawparallels(parallels, labels=[1,0,0,0], ax=self.ax1)
    self.map1.drawmeridians(meridians, labels=[0,0,0,1], ax=self.ax1)

    # calculate center of lake polygons within domain
    polygons = []
    tempcenter = []
    for fshape, lcind in zip(self.map1.lakes, range(len(self.map1.lakes))):
        farray = np.array(fshape)
        featlon, featlat = self.map1(farray[:,0], farray[:,1], inverse=True)
        lcenter = [np.mean(featlon), np.mean(featlat)]
        if ( (lcenter[0] >= self.lonrange[0]-0.05) & \
             (lcenter[0] <= self.lonrange[1]+0.05) & \
             (lcenter[1] >= self.latrange[0]-0.05) & \
             (lcenter[1] <= self.latrange[1]+0.05) ):
          tempcenter.append(lcenter)
          polygons.append(farray)
    self.polygons = polygons
    self.lakecenter = np.array(tempcenter, dtype=float)
    
    # find first lake and go to it on inset map
    self.lakenum = sorted(set(self.lakegrid[self.lakind]))
    if currentLake > 0 and currentLake <= self.lakenum[-1]:
        self.thisLakeInd = np.argmin(np.abs(np.array(self.lakenum)-currentLake))
        print 'Beginning at Lake %i...' %(self.lakenum[self.thisLakeInd])
    elif currentLake > self.lakenum[-1]:
        self.thisLakeInd = len(self.lakenum)-1
        print 'Beginning at Lake %i...' %(self.lakenum[-1])
    else:
        self.thisLakeInd = 0
        print 'Beginning at Lake 1...'
    self.plot_subset(whichLake=1)
    # now ready for mouse and keyboard events

'''
###############################################################################
Author: Nicholas Elmer (UAH/SPoRT)
Date: July 2016
Version: 1.0

Description:
Creates Tkinter window and opens widget for interactive QC of channelgrid and 
  lakegrid from WRF-Hydro ArcGIS Preprocessing Tool.  Use WRF-Hydro_reorder_lakes.py
  to update LAKEPARM.TBL to account for changes and reorder lakes consecutively.

Input:
gridpath,           string,    full path and filename to WRF-Hydro ArcGIS 
                                Preprocessing Tool output folder (unzipped)
lakefile,           string,    full path and filename to lake polygon shape file. 
                                It is suggested that the lake polygons are clipped 
                                to the domain extent prior to use to limit overhead.
map_subset,         boolean,   overlay lake polygons on subset image to assist
                                with quality control (not yet supported).
                                Default is False.
open_saved_grids,   boolean,   open a grid already saved using this tool.
                                Default is False.
currentLake,        integer,   if open_saved_grids = False, the first lake is
                                plotted in the widget.  If open_saved_grids = True,
                                the value must match the number contained in the
                                saved grid's filename.  The lake number specified 
                                by that value will also be plotted first in the
                                widget.
###############################################################################
'''
def main(gridpath, lakefile, map_subset=False, open_saved_grids=False, \
            currentLake=0):
    root = Tk()
    root.wm_title("WRF-Hydro Grid Correction Tool")
    tool = Tool(root, gridpath, lakefile, map_subset=map_subset, \
            open_saved_grids=open_saved_grids, currentLake=currentLake)
    root.mainloop()


if __name__ == '__main__':
    # specify whether to start from scratch or use a saved grid
    open_saved_grids = False 
    currentLake = 1  # must match filename
    
    # specify domain and data path
    hydrodata = '.'
    lakename = 'lakes.shp'
    
    # define path to files
    gridpath = os.join(hydrodata)
    lakepath = os.join(hydrodata)
    lakefile = os.join(lakepath, lakename)
    
    # run the widget
    main(gridpath, lakefile, open_saved_grids=open_saved_grids, currentLake=currentLake)
