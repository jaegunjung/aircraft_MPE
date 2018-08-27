#!/disk41/jjung_linux/util/python/anaconda3/bin/python

# extract model outputs along flight tracks by picking the nearest model box to each flight track point for comparison.
# required input files are :
# [1] CAMx 3D output files
# [2] Aircraft measurement icartt data file
#
# Original script is from zliu.
# Changed according to changes in step0 (2017-6-3, jjung)

####################################################################
#module declaration
import math
import numpy as np
import sys
import datetime
import csv
from jd3 import julian_date , caldat
from pyproj import Proj
from PseudoNetCDF import PNC
from PseudoNetCDF.camxfiles.Memmaps import uamiv
from PseudoNetCDF.icarttfiles.ffi1001 import ffi1001
from itertools import groupby
from numpy import newaxis

#######################################################################
def get_mod_var_model ( infile , var , ijkdh , npts, sfactor) :
    i = ijkdh[:,0] ; j = ijkdh[:,1] ; k = ijkdh[:,2] ; h = ijkdh[:,4]
    data_var4d = uamiv ( infile ).variables[var]
    data_varpts = np.zeros(npts)-999
    for ipt in range (npts):
      if np.min(ijkdh[ipt,:]) >= 0:
        data_varpts[ipt] = data_var4d [int(h[ipt]),int(k[ipt]),int(j[ipt]),int(i[ipt])] * sfactor
    del data_var4d
    return data_varpts
########################################################################
def get_plane_var ( var , file_plane , ijkdh , npts , chkneg) :
    data_plane = ffi1001(file_plane)
    varall = data_plane.variables[var]
    varall = [ x for x in varall ] #array
    data_varpts = np.zeros(npts)-999
    for ipt in range (npts):
      if np.min(ijkdh[ipt,:]) >= 0:
        if chkneg :
          if varall[ipt] >= 0 :
            data_varpts[ipt] = varall[ipt]
        else:
          data_varpts[ipt] = varall[ipt]
    del data_plane
    return data_varpts
########################################################################
def jtime2yyyymmddhh ( jtime ) :

    ctime = caldat ( jtime )
    yyyy = str(ctime[0]).zfill(4) ; mm = str(ctime[1]).zfill(2) ; dd = str(ctime[2]).zfill(2)
    yyyymmdd = yyyy + mm + dd 
    hh = str(ctime[3]).zfill(2)
    return [yyyymmdd, hh]
########################################################################
def diag_plane() :
    outfile = str(sys.argv[1])
    infile_a = str(sys.argv[2])
    infile_z = str(sys.argv[3])
    file_plane_a = str(sys.argv[4])
    file_plane_z = str(sys.argv[5])
    flight_path = str(sys.argv[6])
    mapping = str(sys.argv[7])
    sfactor = float(sys.argv[8]) # To convert model unit from ppm to ppb
    lyyyyjjj = False
    if str(sys.argv[9]).lower() == 'true' :
      lyyyyjjj = True # CAMx file name has yyyyjjj
    obs_spc=[]; model_spc=[]
    with open(mapping, 'r') as csvfile:
      lines = csv.reader(csvfile, delimiter=',', quotechar='|')
      for line in lines:
        obs_spc.append(line[0])
        model_spc.append(line[1].split(":")) 
    nspc = len(model_spc)
    if len(obs_spc) != nspc:
      print ("no. of obs species is not same as no. of model species")
      print ("no. of obs species = {}".format(len(obs_spc)))
      print ("no. of model species = {}".format(nspc))
    obs_alt = 'ALTP' # m
    obs_lon = 'LONGITUDE' # deg
    obs_lat = 'LATITUDE' # deg
    chkneg = True # Drop negative from observation data. Set to false if it is a variable such as longitude.

#   Read the text file that includes ijkdh (icell, jcell, kcell, date, and hour)
    ijkdhs = np.genfromtxt(flight_path)
    npt = len(ijkdhs)
    # No. of days and data points for Model
    day = sorted(list(set(ijkdhs[:,3])))
    nday = len(day)
    npts = [len(list(group)) for key, group in groupby(ijkdhs[:,3])]
    print("day={}".format(day))

    # No. of days and data points for Obs
    oday = sorted(list(set(ijkdhs[:,5])))
    noday = len(oday)
    nopts = [len(list(group)) for key, group in groupby(ijkdhs[:,5])]

    obs_val = np.zeros([nspc,npt])-999; model_val  = np.zeros([nspc,npt])-999
    alt_val = np.zeros(0)
    lon_val = np.zeros(0)
    lat_val = np.zeros(0)

    # Loop for observation data
    icnt = 0
    for iody in range ( noday ) :
        indbeg = sum(nopts[0:iody])
        indend = sum(nopts[0:iody+1])
        ijkdh = ijkdhs [ indbeg:indend , : ]
        date8c = str(int(oday[iody]))
        print ("date8c={}".format(date8c))
        print("nopts[iody]={}".format(nopts[iody]))
        file_plane = file_plane_a + date8c + file_plane_z
        obs_val0 = np.zeros([nspc,nopts[iody]])-999
        for ispc in range (nspc):
          obs_val0[ispc,:] = get_plane_var(obs_spc[ispc],file_plane,ijkdh,nopts[iody],chkneg)
        for ipt in range (nopts[iody]):
          for ispc in range (nspc):
            obs_val[ispc,icnt] = obs_val0[ispc,ipt]
          icnt += 1
        alt_val= np.append(alt_val,get_plane_var(obs_alt,file_plane,ijkdh,nopts[iody],chkneg))
        lon_val= np.append(lon_val,get_plane_var(obs_lon,file_plane,ijkdh,nopts[iody],False))
        lat_val= np.append(lat_val,get_plane_var(obs_lat,file_plane,ijkdh,nopts[iody],chkneg))

    # Loop for model data
    icnt = 0
    for idy in range ( nday ) :
        indbeg = sum(npts[0:idy])
        indend = sum(npts[0:idy+1])
        ijkdh = ijkdhs [ indbeg:indend , : ]
        print("npts[idy]={}".format(npts[idy]))
        date8c = str(int(day[idy]))
        print ("date8c={}".format(date8c))
        infile = infile_a + date8c + infile_z
        if lyyyyjjj:
          yyyy = int(date8c[0:4])
          mm = int(date8c[4:6])
          dd = int(date8c[6:8])
          jdate7c = str(datetime.date(yyyy,mm,dd).strftime("%Y%j"))
          print ("jdate7c={}".format(jdate7c))
          infile = infile_a + jdate7c + infile_z
        #file_plane = file_plane_a + date8c + file_plane_z
        model_val0 = np.zeros([nspc,npts[idy]])-999
        for ispc in range (nspc):
          for isubspc in range(len(model_spc[ispc])-1):
            if isubspc == 0: # if the first subspc, construct model_val0[ispc,:] 
              model_val0[ispc,:] = get_mod_var_model(infile,model_spc[ispc][1],ijkdh,npts[idy],sfactor) #model_spci[ispc][0] is used for representative variable name
            else: # from the second subspc, add each elements
              temp_array = get_mod_var_model(infile,model_spc[ispc][isubspc+1],ijkdh,npts[idy],sfactor)
              for ipt in range (npts[idy]):
                model_val0[ispc,ipt] += temp_array[ipt]
          #obs_val0[ispc,:] = get_plane_var(obs_spc[ispc],file_plane,ijkdh,npts[idy],chkneg)
        for ipt in range (npts[idy]):
          for ispc in range (nspc):
            model_val[ispc,icnt] = model_val0[ispc,ipt]
          icnt += 1
        #alt_val= np.append(alt_val,get_plane_var(obs_alt,file_plane,ijkdh,npts[idy],chkneg))
        #lon_val= np.append(lon_val,get_plane_var(obs_lon,file_plane,ijkdh,npts[idy],False))
        #lat_val= np.append(lat_val,get_plane_var(obs_lat,file_plane,ijkdh,npts[idy],chkneg))

    np.savez ( outfile+'.npz' , day=ijkdhs[:,3] , alts=alt_val , model=model_val , obs=obs_val )
    data_array = np.column_stack((ijkdhs[:,3] , lon_val, lat_val, alt_val))
    modobs_spc = ["" for x in range(nspc)]
    fmt_str = '%d,%1.5f,%1.5f,%1.3f' + ',%1.3f,%1.3f'*nspc
    for ispc in range (nspc) : # combine model and obs species names for text file header
      modobs_spc[ispc] = ','.join([model_spc[ispc][0],obs_spc[ispc]])
      #print("model_val[ispc,:]{}".format(model_val[ispc,:]))
      #print("obs_val[ispc,:]{}".format(obs_val[ispc,:]))
      data_array = np.concatenate([data_array,model_val[ispc,:,newaxis]],axis=1)
      data_array = np.concatenate([data_array,obs_val[ispc,:,newaxis]],axis=1)
    np.savetxt ( outfile+'.txt' , data_array, fmt=fmt_str, header="Day,Lon,Lat,Altitude[m]," + ','.join(modobs_spc[ispc] for ispc in range(nspc)))
####################################################################
if __name__ == '__main__':
    diag_plane()
