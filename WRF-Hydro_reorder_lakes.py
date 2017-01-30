'''
###############################################################################
Author: Nicholas Elmer (UAH/SPoRT)
Date: July 2016
Name: WRF-Hydro_reorder_lakes.py
Version: 1.0

Description:
Reorder lakes consecutively in CHANNELGRID.nc, LAKEGRID.nc and LAKEPARM.TBL
    output from WRF-Hydro ArcGIS Preprocessing Tool.
###############################################################################
'''

import numpy as np
import os.path as os
from netCDF4 import Dataset as netcdf

'''
###############################################################################
Author: Nicholas Elmer (UAH/SPoRT)
Date: July 2016
Version: 1.0

Description:
Reorder lakes consecutively in CHANNELGRID.nc, LAKEGRID.nc and LAKEPARM.TBL
    output from WRF-Hydro ArcGIS Preprocessing Tool.

Input:
filepath,           string,     full path to WRF-Hydro ArcGIS 
                                Preprocessing Tool output folder (unzipped)
inputgrid_suffix,   string,     suffix in input grid filename. The suffix for 
                                the CHANNELGRID and LAKEGRID files must match.
                                Ex. if grid filename is CHANNELGRID_currentLake12.nc,
                                inputgrid_suffix = '_currentLake12'.
outputgrid_suffix,  string,     same as inputgrid_suffix, but for the reordered
                                grid output filename.
newLakefile,        string,     filename of reordered LAKEPARM.TBL file.  
                                If newLakefile = 'LAKEPARM.TBL', the input version 
                                will be overwritten.
###############################################################################
'''

def main(filepath, inputgrid_suffix, outputgrid_suffix, newLakefile):

    lakeparmTBL = os.join(filepath, 'LAKEPARM.TBL')
    newlakeparmTBL = os.join(filepath, newLakefile)
    
    for gridtype in ['CHANNEL', 'LAKE']:
        dstfile = os.join(filepath,    '%sGRID%s.nc' %(gridtype,outputgrid_suffix))
        srcfile = os.join(filepath,    '%sGRID%s.nc' %(gridtype,inputgrid_suffix))
        srcfileobj = netcdf(srcfile, 'r')
        dstfileobj = netcdf(dstfile, 'w', format=srcfileobj.file_format)
        
        oldgrid = srcfileobj.variables['%sGRID' %gridtype][:]
        newgrid = np.copy(oldgrid)
        
        unique = np.unique(oldgrid)
        unique = unique[unique>0]
        nlakes = len(unique)
        if gridtype == 'CHANNEL':
            print 'There are %i lakes to process...' %nlakes
            print 'The lakes were numbered nonconsecutively from %i to %i.' \
                %(np.min(unique), np.max(unique))
        elif gridtype == 'LAKE':
            # check to make sure LAKEPARM.TBL is correct
            with open(lakeparmTBL, 'r') as oldf:
                TBLrows = []
                for line in oldf:
                    TBLrows.append(line)
            oldf.close()
            print len(TBLrows)
            newf = open(newlakeparmTBL, 'w')
            newf.write(TBLrows[0])
            rows_written = 0
            for r in np.arange(len(TBLrows)-1)+1:
                row = TBLrows[r]
                ilake = int(row.split('\t', 1)[0])
                if ilake in unique:
                    newf.write(TBLrows[r])
                    rows_written += 1
                else:
                    continue
            if rows_written != nlakes:
                raise ValueError('Error some place!')
            newf.close()
        else:
            pass
        
        for lakenum, i in zip(unique, np.arange(len(unique))+1):
            if lakenum != i:
                newgrid[oldgrid==lakenum] = i
            else:
                continue
        
        if gridtype == 'CHANNEL':
            print 'The lakes are now numbered consecutively from 1 to %i.' %nlakes
    
        print 'Writing new %sgrid file...' %gridtype.lower()
        
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
            if v_name == '%sGRID' %gridtype:
                outVar[:] = newgrid
            else:
                outVar[:] = varin[:]
        print dstfileobj
        dstfileobj.close()
        srcfileobj.close()
        
    print 'Done!'

if __name__ == '__main__':
    #specify filepath and filenames
    filepathA = '\Users\\nelmer\Documents\\MyResearch\\WRF-Hydro\\Data\\GIStool'
    filepathB = 'NorthAlabama\\WRF-Hydro_routing-grid_100m'
    filepath = os.join(filepathA, filepathB)
    newLakefile = 'LAKEPARMnew.TBL'
    ingrid = '_QC'
    outgrid = '_reorder'
    main(filepath, ingrid, outgrid, newLakefile)
