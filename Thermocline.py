#! /usr/bin/env ipython                                                                              
''' This Script will plot the thermocline depth with latitude'''
from scipy.io import netcdf
#from numba import autojit
from scipy.interpolate import interp1d
from scipy import interpolate
import numpy as np
import matplotlib.pyplot as plt
import os
#import glob
import sys
##~~~~~~~~~~Set up and sanity checks~~~~~~~~~~~## 
''' We have 8 runs:                                                                                  
    3,10,30,100,300,1000,3000,and Closed '''
if len(sys.argv) < 1:
   print '''\n    OOPSIE                                                                             
            usage:  give me a year to calculate over \n'''
   sys.exit(1)
Year = sys.argv[1]
tau=['3','10','30','100','300','1000','3000','Closed']
#Check what files exist                                                                              
''' Let me know what files are missing '''
check=0
runs=[]
for i in range(len(tau)):
    flist=str(tau[i])+'daynokpp/'+Year+'all.nc'
    if not os.path.exists(flist):
       print '\n WARNING: '+flist+' does not exist! (skipping this tau...)'
       check+=0
    else:
       check+=1
       runs.append(i)
if check<1:
   print "\n No files found - no point. I'm stopping this..."
   sys.exit(1)
##~~~~~~~~~~~~~~~~Get the data~~~~~~~~~~~~~~~~~~~~~~###                                              
Runs=np.array(runs)
#Grab a grid file to get model X,Y,Z arrays  
file2=netcdf.netcdf_file(str(tau[Runs[0]])+'daynokpp/grid.nc','r')
Z=file2.variables['Z']
Z=Z[:]
Yc=file2.variables['Y']
Yc=Yc[:]*1
Zp=file2.variables['Zp1']
Zp=Zp[:]
yc=Yc/1000
Tdepthall=np.zeros((len(Runs),len(Yc)))
Z2=interp1d(Z,Z,axis=0)
Znew=np.linspace(int(Z[0]),int(Z[-1]),120)
Zexp=Z2(Znew)
Zp2=interp1d(Zp,Zp,axis=0)
Zpnew=np.linspace(int(Zp[0]),int(Zp[-1]),121)
Zpexp=Zp2(Zpnew)
dT=np.zeros((len(Zexp),len(Yc)))
dz=Zpexp[0:len(Zpexp)-1]-Zpexp[1:len(Zpexp)]
Tdepth=np.zeros(len(Yc))
for i in range(len(Runs)):
    fname=str(tau[Runs[i]])+'daynokpp/'+Year+'all.nc'
    file2read = netcdf.NetCDFFile(fname,'r')
    Temp=file2read.variables['THETA']
    Temp=Temp[:]
    Tav=np.mean(Temp,axis=0)
    Tavlat=np.mean(Tav,axis=2)
    T2=interp1d(Z,Tavlat,axis=0)
    Tnew=np.linspace(int(Z[0]),int(Z[-1]),120)
    Texp=T2(Tnew)
    DT=Texp[0:len(Zexp)-1,:]-Texp[1:len(Zexp),:]
    for j in range(len(Yc)):
        #for k in range(len(Zexp)-1):
        #    dT[k,j]=DT[k,j]/dz[k]
        dT[0:-1,j]=np.apply_along_axis(np.divide,0,DT[:,j],dz[0:-1]) 
        if np.sum(dT[:,j])==0:
            Tdepth[j]=np.nan
            b=0
        else:
            #b=Zexp[dT[:,j]<=0.0001][0]
            b=Zexp[dT[:,j]==np.nanmax(dT[50::,j])][0] 
            Tdepth[j]=b
    plt.plot(yc,Tdepth,linewidth=2)
    plt.title(str(tau[Runs[i]])+"daynokpp. Average depth of thermocline for years "+Year+")")
    plt.xlabel(r'Distance (km)')
    plt.ylabel('Depth (m)') 
    x=( os.path.expanduser('~')+"/Figures/"+str(tau[Runs[i]])+"daynokpp")
    if not os.path.exists(x):
          os.makedirs(x)
    plt.savefig(x+'/Thermocline'+Year+'.png')
    plt.clf()
    Tdepthall[i,:]=Tdepth
plt.plot(yc,np.transpose(Tdepthall))
plt.title(" Average depth of thermocline for years "+Year+")")
plt.xlabel(r'Distance (km)')
plt.ylabel('Depth (m)')
runnames=[]
for i in range(len(Runs)):
    runnames.append(tau[Runs[i]])
plt.legend(runnames,loc='upper center', shadow=True)
x=( os.path.expanduser('~')+"/Figures")
if not os.path.exists(x):
         os.makedirs(x)
plt.savefig(x+'/Thermocline'+Year+'.png')
plt.clf()
