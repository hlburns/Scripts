#! /usr/bin/env ipython 
################################################################                                     
##                                                            ##                                     
##      Eddy General purpose script for MOC calculation       ##                                     
##                        Helen Burns                         ##   
################################################################


''' Please give a run name and the this will then calculate 
    the Eddy MOC and draw a picture for 
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
from pylab import *
#--Take terminal inputs--# 
if len(sys.argv) < 2:
   print '''\n    OOPSIE                                                                             
            usage: *all.nc *Psi.nc grid.nc Temp must be here'''
   sys.exit(1)
OP = sys.argv[1]
#--Set folder structure--# 
x=os.getcwd()
lists=glob.glob(x+'/*Psi.nc')
#--Main For loop--#                                                           
#For every .nc file in the folder                                             
#Read in netcdf variables                                                     
#Decide resolution and time step                                              
#Find grid info                                                              
#Calculate the streamfunction                                                 
#Draw a picture and save it
if not lists:
   print '''\n    OOPSIE:                                                                            
            Wrong folder you ninny!*all.nc *Psi.nc grid.nc Temp must be here'''
   sys.exit(1)
def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx] 
for file in lists:
    file2read = netcdf.NetCDFFile(file,'r')
    lvrho=file2read.variables["LaVH1TH"]
    lvrho=lvrho[:]
    Y=file2read.variables['Yp1']
    Y=Y[:]
    time=file2read.variables['T']
    ti=time[:]
    #Stop memory issues!!
    if len(ti)>300:
       continue
    dx=Y[1]-Y[0] # Find Resolution
    VT=np.sum(lvrho*dx,axis=3) #integrate Vdx along x
    VTfdz=np.cumsum(VT[:,::-1,:],axis=1) #sum up the water column
    psi=np.mean(VTfdz[:,::-1,:],axis=0)/10**6 #Time average and put into Sv and put back in right order
    y=Y/1000
    Rho = np.genfromtxt('Temp', delimiter = ',') 
    nolayers=len(psi[:,1])
    Rho=Rho[0:nolayers]#The layers package bins a layer so adjust for that
    start=int(np.divide(ti[0],(86400*360)))#Find run start and stop times
    end=int(np.divide(ti[-1],(86400*360)))
    #Now find the corresponding all file! (turn to list and replace Psi with all!)
    a=list(file)
    a[-6]='a'
    a[-5]='l'
    a[-4]='l'
    A=''.join(a)
    file2 = netcdf.NetCDFFile(A,'r')
    Temp=file2.variables['THETA']
    Temp=Temp[:]
    Tav=np.mean(Temp,axis=0)
    V=file2.variables['VVEL']
    V=V[:]*1
    file2 = netcdf.NetCDFFile('grid.nc','r')
    lm=file2.variables['HFacC']
    Yc=file2.variables['Y']
    Yc=Yc[:]
    Z=file2.variables['Z']
    Z=Z[:]
    X=file2.variables['X']
    X=X[:]
    Z=Z[:]
    lm=lm[:]
    lmc=np.array(lm)
    lmc[lmc<1]=np.nan
    Tavlat=np.mean(Tav,axis=2)
    Yc=Yc/1000
    lmav=np.mean(lmc,axis=2)
    Vtave=np.mean(V,axis = 0)
    Zp1=file2.variables['Zp1']
    Zp=Zp1[:]
    Vtave[Vtave==0]=np.nan
    Vzone=np.nansum(Vtave*dx,axis = 2)
    dz=Zp[0:len(Zp)-1]-Zp[1:len(Zp)]
    # No more super slow forloop!
    psi2=np.apply_along_axis(np.multiply,0,Vzone,dz)
    psi3=np.cumsum(-psi2[::-1,:],axis=0)
    npad = ((0,1), (0,0))
    psi4 = np.pad(psi3, pad_width=npad, mode='constant', constant_values=0)
    y =Y/1000
    Psi=psi/10**6 #Convert to Sv 
    Psi2=psi4/10**6
    Psimap=np.zeros(np.shape(Psi2))
    for i in range(len(Yc)):
        for k in range(len(Z)):
            #D=Densav[k,i]
            D=Tavlat[k,i]
            if np.isnan(D):
               Psimap[k,i]=np.nan
            else:
                P=psi[:,i]
                I=find_nearest(Rho, D)
                b=nonzero(Rho==I)[0][0]
                Psimap[k,i]=P[b]
    #Psimapped=np.multiply(Psimap,lmav) HFacS is req for topo runs
    Psimap[Psimap==0]=np.nan
    Psied=Psimap-Psi2
    x=( os.path.expanduser('~')+"/Figures/"+OP)
    if not os.path.exists(x):
          os.makedirs(x)
    contourf(y,Zp,Psied,50,cmap=cm.seismic) #Use b2r colourmap                                       
    clim(-2.5,2.5) # Put 0 to white                           
    cbar = colorbar()
    title("Eddy MOC year "+str(start)+" "+OP)
    xlabel('Distance (km)')
    ylabel('Density')
    cbar.ax.set_ylabel('Psi (sv)')
    plt.savefig(x+"/PsiEddy"+str(start)+"-"+str(end)+".png")
    plt.clf()
    plot1=plt.contourf(y,Zp,Psimap,150,cmap=cm.seismic)
    cbar = plt.colorbar()
    clim(-1,1)
    plt.title("Remapped Stream Function"+str(start)+OP)
    plt.xlabel('Distance (km)')
    plt.ylabel('Depth')
    cbar.ax.set_ylabel('$\psi_{res}$ (Sv)')
    plt.savefig(x+"/Psiremapped"+str(start)+"-"+str(end)+".png")
    plt.clf()    

