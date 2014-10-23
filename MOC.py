#! /usr/bin/env ipython 
################################################################
##                                                            ##
##        MOC General purpose script for MOC calculation      ##      
##                        Helen Burns                         ##
################################################################


''' Please give a run name in order to file away the output in ~/Figures/Runname      
    This script will then calculate the MOC and draw a picture for all files in 
    that folder and save the figures in ~/Figures/Folder.  
                               
    For spinup timeseries please use the Spin up diagnostics notebook.'''                  
      
###################################################################                                  
###################################################################                                  
#--Import modules--#
from scipy.io import netcdf
import numpy as np
from pylab import clf, plot, show, floor, ceil, imshow
from matplotlib.pyplot import xlabel, ylabel, legend, savefig, colorbar, title, clim, pcolor, cm, contourf
import os
import csv
import sys
import glob
#--Take terminal inputs--#
if len(sys.argv) < 2:
   print '''\n    OOPSIE  
            usage: Go to the folder with netcdf V output 
            (named *all.nc with a grid folder located ../ or
            within the folder (for working directories) \n'''
   sys.exit(1)
OP = sys.argv[1]
#--Set folder structure--#
x=os.getcwd()
lists=glob.glob(x+'/'+str(OP)+'/*all.nc')
#--Main For loop--#
#For every .nc file in the folder
#Read in netcdf variables
#Decide resolution
#Find grid info
#Calculate the streamfunction
#Draw a picture and save it
if not lists:
   print '''\n    OOPSIE: 
            no *all.nc folders here!!      
            usage: Go to the folder with netcdf V output                                            
            (named *all.nc with a grid folder located ../ or                                                                                
            within the folder (for working directories) \n'''
   sys.exit(1)
file2=netcdf.netcdf_file(x+'/'+str(OP)+'/grid.nc','r')
Zp1=file2.variables['Zp1']
Zp=Zp1[:]*1
file2=netcdf.netcdf_file(x+'/'+str(OP)+'/grid.nc','r')
Z=file2.variables['Z']
Z=Z[:]
Y=file2.variables['Yp1']
Y=Y[:]
X=file2.variables['X']
X=X[:]
if len(X)>400 :
    Q_levs = (np.arange(-0,150,10.))
    Psi_levs = Q_levs / 10
    Q_ticks = np.arange(0,150,20.)
    Psi_ticks = Q_ticks / 10
else:
    Q_levs = (np.arange(-0,30)+0.5)
    Psi_levs = Q_levs / 10
    Q_ticks = np.arange(0,30,3.)
    Psi_ticks = Q_ticks / 10
for file in lists:
       file2read = netcdf.NetCDFFile(file,'r')
       V=file2read.variables['VVEL']
       V=V[:]
       #Find timestep and resolution
       time=file2read.variables['T']
       ti=time[:]
       dx=Y[1]-Y[0]
       Vtave=np.mean(V,axis = 0)
       #Vtave[Vtave==0]=np.nan
       Vzone=np.nansum(Vtave,axis = 2)*dx
       dz=Zp[0:len(Zp)-1]-Zp[1:len(Zp)]
       # Got rid of for loop here (much quicker!!)                                             
       psi2=np.apply_along_axis(np.multiply,0,Vzone,dz)
       psi=np.cumsum(-psi2[::-1,:],axis=0)
       npad = ((0,1), (0,0))
       # Pad with zeros at bottom                                                                    
       psi = np.pad(psi[::-1,:], pad_width=npad, mode='constant', constant_values=0)
       y =Y/1000
       Psi=psi/10**6 #Convert to Sv
       start=int(np.divide(ti[0],(86400*360)))
       end=int(np.divide(ti[-1],(86400*360)))
       cf=contourf(y,Zp,Psi,Psi_levs,cmap=cm.seismic) #Use b2r colourmap
       clim(-15,15) # Put 0 to white
       cbar = colorbar(cf, ticks=Psi_ticks, shrink=0.8)
       title("MOC years "+str(start)+"-"+str(end))
       xlabel('Distance (km)')
       ylabel('Depth (m)')
       cbar.ax.set_ylabel('$\overline{\psi} \,\, (sv)$')
       x=( os.path.expanduser('~')+"/Figures/"+OP)
       if not os.path.exists(x):
          os.makedirs(x)
       y=x+"/MOC"+str(start)+"-"+str(end)+".png"
       savefig(y)
       clf()
       
