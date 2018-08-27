#!/usr/bin/env python

# Draw scatter plots and vertical profiles
# Original script is from zliu.
# Altered to reflect changes in .npz (7/3/2017 jjung)

import math
import csv
import numpy as np
import os
import sys
from scipy import stats
from matplotlib import pyplot as plt

def get_vertprof(data,alt,altbins) :
    alt  = alt[data>0]
    data = data[data>0]
    ind_alts = np.digitize(alt, altbins)
    data_altbin_nspls = [data[ind_alts == i].shape[0] for i in range(1, len(altbins))]
    data_altbin_means = [np.mean(data[ind_alts == i]) for i in range(1, len(altbins))]
    data_altbin_medians = [np.median(data[ind_alts == i]) for i in range(1, len(altbins))]
    data_altbin_25 = [np.percentile(data[ind_alts == i],25) for i in range(1, len(altbins))]
    data_altbin_75 = [np.percentile(data[ind_alts == i],75) for i in range(1, len(altbins))]

    return [data_altbin_means,data_altbin_25,data_altbin_medians,data_altbin_75]

def get_mpe_linregress ( obsdata , moddata ) :
    obsdata = obsdata[obsdata != -999 ]
    moddata = moddata[moddata != -999 ]
    obsdata = np.asarray(obsdata) ; moddata = np.asarray(moddata)
    slope, intercept, r_value, p_value, std_err = stats.linregress(obsdata,moddata)
    return slope, intercept, r_value, p_value, std_err

def get_mpe_stats ( obsdata , moddata ) :
    obsdata = obsdata[obsdata != -999 ]
    moddata = moddata[moddata != -999 ]
    obsdata = np.asarray(obsdata) ; moddata = np.asarray(moddata)
    nmb = np.sum(moddata-obsdata)/(np.sum(obsdata))*100.
    nme = (np.sum(np.abs(moddata-obsdata)))/(np.sum(obsdata))*100.

    #r2  = (np.corrcoef(moddata,obsdata)[0])*(np.corrcoef(moddata,obsdata)[0]) 
    return nme,nmb#,r2    

def plot_vertprof_scatters ( data , flight , var, var_ord, runname, ovar, unit, scale) :
   day = data['day']
   if flight != 'ALL_FLIGHTS' : 
      condition = (day - int(flight)) == 0
   else :
      condition = (day - day) == 0

   model=data['model']
   obs=data['obs']
   alts=data['alts']
   print ("obs={}".format(obs))
   print ("data['obs']={}".format(data['obs']))
   print ("data['obs'][2][condition]={}".format(data['obs'][var_ord][condition]))
   print ("len(data['day'])={}".format(len(data['day'])))
   print ("var={}".format(var))
   yy = data['alts'][condition]
   npt = len(data['day'])
   xx0 = data['obs'][var_ord][condition]*scale
   xx1 = data['model'][var_ord][condition]
   for ipt in range (npt):
     if min(xx0[ipt],xx1[ipt]) < 0: # If either model or obs is missing, set both missing to plot
       xx0[ipt] = -999.; xx1[ipt] = -999.
   altbins = [0,500,1000,1500,2000,2500,3000,3500,4000,4500,5000,5500,6000,6500,7000,7500]
   altaxis = [0.5*(altbins[i]+altbins[i+1]) for i in range (len(altbins)-1)]
   varlim =  [0,np.max([xx0,xx1])]
   campname = '$\mathregular{SAS}$'
   # Specify run name
   if runname == 'megan2.1_org_01':
     run4title = 'CAMx with MEGAN v2.1'
   elif runname == 'megan3_02':
     run4title = 'CAMx with MEGAN v3'
   else:
     print ("your runname is {}".format(runname))
     print ("This runname is not appropriate.")
     exit()
   if var =='O3': 
      varlabel = '$\mathregular{O_3 [' + unit + ']}$'
      figtitle = '$\mathregular{O_3}$ during '+ campname +':: Aircraft measurement vs. ' + run4title + ' ' + flight
   else:
      varlabel = "$\mathregular{" + var + "[" + unit + "]}$"
      figtitle = "$\mathregular{" + var + "}$ during "+ campname +":: Aircraft measurement vs. " + run4title + ' ' + flight

   slope_model1, intcpt_model1, r_val_model1, p_val_model1, std_err_model1 = get_mpe_linregress ( xx0 , xx1 )
   r2_model1 = r_val_model1 * r_val_model1 
   
   nme_model1, nmb_model1 = get_mpe_stats(xx0,xx1)  
   if np.max(xx1)>0 :\
       prof_model1 = get_vertprof(xx1,data['alts'][condition],altbins)
   if np.max(xx0)>0 :\
       prof_p3b  = get_vertprof(xx0,data['alts'][condition],altbins)


   fig = plt.figure(figsize=(10, 5))
   fig.subplots_adjust(left=0.1, right=0.95, wspace=0.25,bottom=0.15, top=0.9)
   fig.suptitle(figtitle)


   # Vertical Profiles
   ax = fig.add_subplot(121)

   ax.plot(prof_model1[1],altaxis,color='red',linewidth=2,linestyle=':')
   #ax.plot(prof_p3b[1],altaxis,color='black',linewidth=2,linestyle=':')
   ax.fill_betweenx(altaxis,prof_p3b[1],prof_p3b[3],facecolor='grey',alpha=0.5)
   ax.plot(prof_model1[2],altaxis,color='red',linewidth=3,linestyle='-')
   ax.plot(prof_p3b[2],altaxis,color='black',linewidth=3,linestyle='-')
   ax.plot(prof_model1[3],altaxis,color='red',linewidth=2,linestyle='--')
   #ax.plot(prof_p3b[3],altaxis,color='black',linewidth=2,linestyle='--')

   ax.set_xlim(0,np.max([prof_model1,prof_p3b])*1.1)

   ax.set_ylim(0,np.max(altbins))
   ax.set_xlabel(varlabel)
   ax.set_ylabel('Altitude [m]')
   ax.annotate('SAS',xy=(0.7,0.10),color='black',xycoords='axes fraction')
   ax.annotate(run4title,xy=(0.7,0.05),color='red',xycoords='axes fraction')
   ax.annotate('25%(dots), 50%(solid), and 75%(dashed) percentiles',\
                xy=(0.1,0.95),xycoords='axes fraction',fontsize=8) 

   # Scatter plots
   ax = fig.add_subplot(122)
   ax.plot(xx0, xx1, marker='.',markersize=2,alpha=0.5,linestyle='none',color='red')
   ax.set_xlim(varlim[0], varlim[1])
   ax.set_ylim(varlim[0], varlim[1])
   ax.set_ylabel('Modeled '+varlabel)
   ax.set_xlabel('Observed '+varlabel)
   ax.annotate(run4title + ': NME=%.1f%% NMB=%.1f%% $\mathregular{R^2}$=%.3f' %(nme_model1,nmb_model1,r2_model1),\
                xy=(0.06,0.85),color='red',xycoords='axes fraction',fontsize=8)


   ax.plot(varlim,varlim,linestyle='--',color='black')
   ax.plot(varlim,varlim,linestyle='--',color='black')
   ax.annotate('y=x',xy =(0.85,0.95),xycoords='axes fraction',fontsize=8) 
   ax.plot(varlim,[i * 2 for i in varlim],linestyle='--',color='black')
   ax.annotate('y=2x',xy =(0.5,0.95),xycoords='axes fraction',fontsize=8)
   ax.plot(varlim,[i * .5 for i in varlim],linestyle='--',color='black')
   ax.annotate('y=0.5x',xy =(0.82,0.35),xycoords='axes fraction',fontsize=8)
   ax.plot(varlim,[i*slope_model1+intcpt_model1 for i in varlim ],linestyle='-',color='red',linewidth = 2)

   if not os.path.exists('./png/'+runname): os.makedirs('./png/'+runname)
   plt.savefig('./png/'+runname+'/'+flight+'_'+ovar+'_vert_scatt.png', dpi=300)
   plt.close()

def plot():
   mapping = str(sys.argv[1])
   runname = str(sys.argv[2])
   unit = str(sys.argv[3])
   ovars = []; vars = []
   scales = []
   with open(mapping, 'r') as csvfile:
     lines = csv.reader(csvfile, delimiter=',', quotechar='|')
     for line in lines:
       ovars.append(line[0]) # Create a list of observation species 
       scales.append(float(line[2])) # Create a list of observation species units
       vars.append(line[1].split(":")[0]) # Create a list of model representative species
   data =  np.load ('./2_extract/out_data/camx/vars_' + runname + '.vs.plane_mrg60.npz')
   #flights = ['ALL_FLIGHTS','20130603','20130605',...]
   flights = ['ALL_FLIGHTS']

   for m in range (len(vars)):
     var = vars[m]
     var_ord = m
     scale = scales[m]
     ovar = ovars[m]
     print (var)
     print (scale)
     for n in range ( len(flights) ) :
       flight = flights[n]
       print (flight)
       plot_vertprof_scatters( data , flight , var, var_ord, runname, ovar, unit, scale)

if __name__ == '__main__':
    plot()
