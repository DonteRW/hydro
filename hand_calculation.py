# -*- coding: utf-8 -*-
"""
Created on Mon Aug 15 12:23:30 2016

@author: nelmer
"""

import arcpy, os, traceback, sys, numpy
from arcpy import env
from arcpy.sa import *

def showPyMessage():
    arcpy.AddMessage(str(time.ctime()) + " - " + message)

#execfile("C:\Users\\nelmer\Documents\PythonScripts\hand_calculation.py")
norder=1
fdir="\Users\\nelmer\Documents\\HAND\\flowdirection"
dirArray = arcpy.RasterToNumPyArray(fdir,"","","",-9999)
nRows,nCols=dirArray.shape
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
for iorder in [1]: #numpy.arange(1,norder,1):
    demG="\Users\\nelmer\Documents\\HAND\\stream_elev%i" %(iorder)
    d=arcpy.Describe(demG)
    ##SR=d.spatialReference
    try:
        dem = arcpy.RasterToNumPyArray(demG,"","","",-9999)
        blankArray=arcpy.RasterToNumPyArray(demG,"","","",-9999)
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
        oneGrid.save("\Users\\nelmer\Documents\\HAND_NWM\\dist_hand%i" %(iorder))
    except:
        print 'Error some place'
        message = "\n*** PYTHON ERRORS *** "; showPyMessage()
        message = "Python Traceback Info: " + traceback.format_tb(sys.exc_info()[2])[0]; showPyMessage()
        message = "Python Error Info: " +  str(sys.exc_type)+ ": " + str(sys.exc_value) + "\n"; showPyMessage()