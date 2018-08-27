#!/disk41/jjung_linux/util/python/anaconda3/bin/python

# step0 return (i,j,k,hr) of model results from air plane path
# required input files are :
# [1] CAMx 3D met file
# [2] CAMx landuse file which has 'TOPO_M' 
# [3] Aircraft measurement icartt data file
#
# Original script is from zliu.
# Change calculation of finding layer index that flight path passed from pressure to height (2017-5-18, jjung)

####################################################################
#module declaration
import math
import numpy as np
import sys
from jd3 import julian_date , caldat
from pyproj import Proj
from PseudoNetCDF import PNC
from PseudoNetCDF.camxfiles.Memmaps import uamiv
from PseudoNetCDF.icarttfiles.ffi1001 import ffi1001

####################################################################
def get_ijk ( lcc, dxy, nx, ny, height_mod , lat_plane , lon_plane , height_plane , lemis ) :
    #if not (math.isnan(height_plane) or math.isnan(lon_plane) or math.isnan(lat_plane)):
    if not (math.isnan(lon_plane) or math.isnan(lat_plane)):
      lcpx, lcpy = lcc(lon_plane,lat_plane)
      ii = int(lcpx/dxy) # minus one grid cell as python starts from 0
      jj = int(lcpy/dxy) # minus one grid cell as python starts from 0
      if (ii < 0) or (ii >= nx) or (jj < 0) or (jj >= ny):
        print ("i or j index is out of domain")
        print ("i, j, nx, ny = {}, {}, {}, {}".format(ii+1, jj+1, nx, ny))
        print ("This point will be ignored")
        ii = -999 ; jj = -999 ; kk = -999
      #print("ii,jj,nx,ny={},{},{},{}".format(ii,jj,nx,ny))
      else: # the plane is within domain horizontally
        top_hgt1d = height_mod[:,jj,ii]
        nz = height_mod.shape[0]
        if not lemis: 
          if (height_plane > top_hgt1d[nz-1]) or math.isnan(height_plane) :
            print ("Airplane is located above the model top or missing")
            print ("Airplane altitude = {}".format(height_plane))
            print ("Model top (sum of topo and layer height) = {}".format(top_hgt1d[nz-1]))
            ii = -999 ; jj = -999 ; kk = -999
          else:
            kk = np.searchsorted(top_hgt1d,[height_plane,],side='right')[0]
        else: # altitude will not be checked for emissions evaluation. It is always assigned to the first layer 
          kk = 0
    else :
      ii = -999 ; jj = -999 ; kk = -999
    #print  ii , jj , kk
    ijk = [ii, jj, kk]
    return ijk

#######################################################################
def get_plane_xyht (file_plane,itzon,yyyy,mm,dd) :
    data_plane = ffi1001(file_plane)
    lats = data_plane.variables['LATITUDE']
    lons360 = data_plane.variables['LONGITUDE']
    lons = [((x + 180) % 360) - 180 for x in lons360 ]
    height = data_plane.variables['ALTP']
    height = [ x for x in height ] #array
    fday = data_plane.variables['UTC'][:]/86400. - itzon/24.
    jday = [ x + julian_date (yyyy,mm,dd, 0 , 0 , 0 ) for x in fday ]
    #print ('jday[0]={}'.format(jday[0]))
    return lats , lons , height , jday

###################################################################
def jtime2yyyymmddhh ( jtime ) :

    ctime = caldat ( jtime )
    yyyy = str(ctime[0]).zfill(4) ; mm = str(ctime[1]).zfill(2) ; dd = str(ctime[2]).zfill(2)
    yyyymmdd = yyyy + mm + dd 
    hh = str(ctime[3]).zfill(2)
    return [yyyymmdd, hh]

####################################################################
def diag_plane() :

    lemis = False
    if str(sys.argv[1]).lower() == 'true' :
      lemis = True
    outfile = str(sys.argv[2])
    met3d_a = str(sys.argv[3])
    met3d_z = str(sys.argv[4])
    file_met2d = str(sys.argv[5])
    file_plane = str(sys.argv[6])
    itzon = int(sys.argv[7])
    date8c = str(sys.argv[8])
    file_met3d = met3d_a + date8c + met3d_z

#   Read the aircraft flight track coordinates ( lat , lon , height , julian time )
    yyyy = int(date8c[:4])
    mm = int(date8c[4:6])
    dd = int(date8c[6:8])

    lats , lons , height , jtimes = get_plane_xyht ( file_plane , itzon, yyyy, mm, dd)

    npts = len (lats)
    print ('no. of data points = {}'.format(npts))

    ijkdhs = np.zeros((npts,6))

    # Get the model grid cell that includes the current point on track at the hour
    data_met3d = uamiv( file_met3d )
    data_met2d = uamiv( file_met2d )

    # Extract domain defintion from data_met3d
    nx = len(data_met3d.dimensions['COL'])
    ny = len(data_met3d.dimensions['ROW'])
    nl = len(data_met3d.dimensions['LAY'])
    dxy = float(getattr(data_met3d,'XCELL')) # Simply assume XCELL = YCELL
    lon0 = str(getattr(data_met3d,'XCENT'))
    lat0 = str(getattr(data_met3d,'YCENT'))
    lat1 = str(getattr(data_met3d,'P_ALP'))
    lat2 = str(getattr(data_met3d,'P_BET'))
    x0 = str(-getattr(data_met3d,'XORIG'))
    y0 = str(-getattr(data_met3d,'YORIG'))
    lcc = Proj('+proj=lcc +a=6370000, +b=6370000, +lon_0='+lon0+' +lat_0='+lat0+' +lat_1='+lat1+' +lat_2='+lat2+' +x_0='+x0+' +y_0='+y0)

    # Read height variables
    height0 = np.asarray ( data_met3d.variables['ZGRID_M'] ) #meter
    topo    = np.asarray ( data_met2d.variables['TOPO_M'] ) #meter
    heights_mod  = np.zeros ( data_met3d.variables['ZGRID_M'][0,:,:,:].shape,data_met3d.variables['ZGRID_M'].dtype ) #meter

#   Loop through all the data points along flight tracks   
    for i in range ( npts ) :
        lat_plane = lats[i] ; lon_plane = lons[i] ; height_plane = height[i] ; jtime = jtimes[i]     
        #print ('jtime= {}'.format(jtime))
        yyyymmddhh = jtime2yyyymmddhh ( jtime )
        yyyymmdd = yyyymmddhh [ 0 ]
        hour = int(yyyymmddhh[ 1 ])
        
        for il in range(nl):
          heights_mod [il,:,:] = height0 [hour,il,:,:] + topo [0,0,:,:]

        ijk = get_ijk ( lcc, dxy, nx, ny, heights_mod , lat_plane , lon_plane , height_plane , lemis ) 
        ii = ijk[0] ; jj = ijk[1] ; kk = ijk[2]             
        #if np.min (ijk) >0 :
          # grid index reduced by one intentionally as python starts from 0. (2015-5-18, jjung) 
          #print (i , yyyymmdd , ii , jj , kk , hour , lat_plane , lon_plane , height_plane ,heights_mod[kk,jj,ii])
        ijkdhs[i,:] = ii,jj,kk,yyyymmdd,hour,int(date8c)

    np.savetxt(outfile,ijkdhs,fmt="%d")

####################################################################
if __name__ == '__main__':
    diag_plane()
