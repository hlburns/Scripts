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
    flist=str(tau[i])+'daynokpp/'+Year+'Surf.nc'
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
X=file2.variables['X']
X=X[:]
Yc=file2.variables['Y']
Yc=Yc[:]*1
yc=Yc/1000
hmlall=np.zeros((len(Runs),len(Yc)))
for i in range(len(Runs)):
    fname=str(tau[Runs[i]])+'daynokpp/'+Year+'Surf.nc'
    file2read = netcdf.NetCDFFile(fname,'r')
    hml=file2read.variables['MXLDEPTH']
    hml=hml[:] 
    hml=np.squeeze(hml)
    hmlav=np.mean(hml,1) 
    plt.plot(yc,hml,linewidth=2)
    plt.title(str(tau[Runs[i]])+"daynokpp. Average mixlayer depth for years "+Year+")")
    plt.xlabel(r'Distance (km)')
    plt.ylabel('Depth (m)') 
    x=( os.path.expanduser('~')+"/Figures/"+str(tau[Runs[i]])+"daynokpp")
    if not os.path.exists(x):
          os.makedirs(x)
    plt.savefig(x+'/Hml'+Year+'.png')
    plt.clf()
    hmlall[i,:]=hmlav
plt.plot(yc,-np.transpose(hmlall))
plt.title(" Average depth of mixedlayer for years "+Year)
plt.xlabel(r'Distance (km)')
plt.ylabel('Depth (m)')
runnames=[]
for i in range(len(Runs)):
    runnames.append(tau[Runs[i]])
plt.legend(runnames,loc='lower center', shadow=True)
x=( os.path.expanduser('~')+"/Figures")
if not os.path.exists(x):
         os.makedirs(x)
plt.savefig(x+'/Hml'+Year+'.png')
plt.clf()
