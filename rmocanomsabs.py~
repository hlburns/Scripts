#! /usr/bin/env ipython 
################################################################                                     
##                                                            ##                                     
##      Anomalies for RMOC from the 3day relaxtion time       ##                                     
##                        Helen Burns                         ##   
################################################################

''' Please give a run name and the this will then calculate 
    the RMOC in temperature layers and draw a picture for 
    all files in that folder and save the figures in ~/Figures/Folder.                               

    For spinup timeseries please use the Spin up diagnostics notebook.''' 
###################################################################   
###################################################################   
#--Import modules--# 
from scipy.io import netcdf
import numpy as np
import matplotlib.pyplot as plt
import os
import csv
import sys
import glob
from Calculations import remap
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
       if i==0:
           print '\n ERROR: '+flist+' does not exist! Are you sure you want to continue?\n y or n:'
           answer = raw_input() 
           if answer == 'n':
              print '\n Stopping...'
              sys.exit(1)  
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
Zp1=file2.variables['Zp1']
Zp=Zp1[:]*1
Z=file2.variables['Z']
Z=Z[:]
Y=file2.variables['Yp1']
Y=Y[:]*1
y=Y/1000
X=file2.variables['X']
X=X[:]*1
Yc=file2.variables['Y']
Yc=Yc[:]*1
yc=Yc/1000
#Plotting things                                                                                  
Q_levs = (np.arange(-25,25)+0.5)
Psi_levs = Q_levs / 10
Q_ticks = np.arange(-8,9,2.)
Psi_ticks = Q_ticks / 10
# Find Base RMOC
fname3day=str(tau[Runs[0]])+'daynokpp/'+Year+'all.nc'
Rho = np.genfromtxt(str(tau[Runs[0]])+'daynokpp/Temp', delimiter = ',')
print '\n Calculating base RMOC...'
Psimap1,Psied,Zexp=remap(fname3day,Rho,X,Yc,Z,Y,Zp)
print '\n OK '+str(tau[Runs[0]])+'day set as base see figure for reference. \n Im now calculating the annomalies, this should take a while ....'
#Find annomalies!
for i in range(len(Runs)):
    fname=str(tau[Runs[i]])+'daynokpp/'+Year+'all.nc' 
    Psimap,Psied,Zexp=remap(fname,Rho,X,Yc,Z,Y,Zp)
    Psimapanom=Psimap-Psimap1
    cf=plt.contourf(yc,Zexp,Psimapanom,Psi_levs,cmap=plt.cm.seismic)
    cbar = plt.colorbar(cf, ticks=Psi_ticks, shrink=0.8)
    plt.clim(-1,1)
    plt.title(str(tau[Runs[i]])+"day Remapped RMOC annomally")
    plt.xlabel('Distance (km)')
    plt.ylabel('Depth (m)')
    cbar.ax.set_ylabel('$\psi_{res}$ anomally (Sv)')
    x=( os.path.expanduser('~')+"/Figures/")
    plt.savefig(x+str(tau[Runs[i]])+'daynokpp/'+str(Year)+'RMOCanom'+str(tau[Runs[i]])+".png")
    plt.clf() 
    print '\n'+str(tau[Runs[i]])+'day done.'
