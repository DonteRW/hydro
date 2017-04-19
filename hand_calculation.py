# -*- coding: utf-8 -*-
"""
Created on Mon Sep 12 12:51:05 2016

@author: nelmer
"""

import arcpy, os, traceback, sys, numpy
from arcpy import env
from arcpy.sa import *

#execfile("C:\Users\\nelmer\Documents\PythonScripts\hand_calculation_orig.py")
norder=1 #stream order - should be 1 unless running modified HAND method
dir="\Users\nelmer\Documents\HAND\"
fdir=os.path.join(dir, 'flowdirection')
strfname = 'stream_elev' #Raster containing elevation of channel network.
                         #Can be calculated from channel network using Raster Calculator:
                         #Ex. Con("channelgrid" == 1, "DEM")
                         #Where channelgrid is your channel network raster, which defines
                         # a channel pixel as 1 and a non-channel pixel as NoData.

for iorder in numpy.arange(norder)+1:
  print 'streamorder %i' %iorder
  try:
    demG=os.path.join(dir, strfname)

    d=arcpy.Describe(demG)
    ##SR=d.spatialReference
    def showPyMessage():
        arcpy.AddMessage(str(time.ctime()) + " - " + message)

    dirArray = arcpy.RasterToNumPyArray(fdir,"","","",-9999)
    dem = arcpy.RasterToNumPyArray(demG,"","","",-9999)
    nRows,nCols=dirArray.shape
    blankArray=arcpy.RasterToNumPyArray(demG,"","","",-9999)
    cellsTotal=nCols*nRows
    d=arcpy.Describe(fdir)
    origin=d.extent.lowerLeft
    cSize=arcpy.Raster(fdir).meanCellHeight
##  directions to find cell neighbour
    fDirs=(1,2,4,8,16,32,64,128)
    dCol=(1,  1,  0, -1, -1,-1, 0,1)
    dRow=(0, -1, -1, -1,  0,  1, 1,1)
##  flipped 
    dRow=(0,  1,  1,  1,  0, -1, -1,-1)
##  main loop
    arcpy.SetProgressor("step", "", 0, nRows)
    for nRow in range (nRows):
            for nCol in range (nCols):
                S=dirArray[nRow,nCol]
                if S in (-1,-9999): continue
                nR,nC=nRow,nCol
                cells=[(nR,nC)]
                z=dem[nR,nC]
                while True:
                    direction=dirArray[nR,nC]
                    i=fDirs.index(direction)
                    dX=dCol[i];nC+=dX
##                    if nC in (0,nCols): break
                    if nC not in range(nCols): break
                    dY=dRow[i];nR+=dY
##                    if nR in (0,nRows): break
                    if nR not in range(nRows): break
                    cells.append((nR,nC))
                    S=dirArray[nR,nC]
                    z=dem[nR,nC]
                    if S<-9998 or z>-9998:
                        break
                    if S==-1:
                        z=blankArray[nR,nC]
                        break
                for nR,nC in cells:
                    blankArray[nR,nC]=z
                    dirArray[nR,nC]=-1
            arcpy.SetProgressorPosition()
    myRaster = arcpy.NumPyArrayToRaster(blankArray,origin,cSize,cSize)
    oneGrid=Con(myRaster<>-9999,myRaster)
    oneGrid.save(os.join(dir,'mapped_DEM_HAND_%i' %iorder))
    del dem,blankArray
  except:
    message = "\n*** PYTHON ERRORS *** "; showPyMessage()
    message = "Python Traceback Info: " + traceback.format_tb(sys.exc_info()[2])[0]; showPyMessage()
    message = "Python Error Info: " +  str(sys.exc_type)+ ": " + str(sys.exc_value) + "\n"; showPyMessage()