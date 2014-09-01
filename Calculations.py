from scipy.io import netcdf
from numba import autojit
import numpy as np
import matplotlib.pyplot as plt
import os
import glob
import sys
from scipy.interpolate import interp1d
from scipy import interpolate
## MOC ##                                                                                            
def MOC(file,Y,Z,Zp):
       '''
       Calculates MOC from the filename (to read in
       changing V) and Yp1, Z and Zp from the grid
       file. In the future you should add in to read
       Y,Z,Zp from the file if they're not already 
       there. 
       Outputs Psi in Sv and Initial time of file
       in years.
       '''
       file2read = netcdf.NetCDFFile(file,'r')
       V=file2read.variables['VVEL']
       V=V[:]
       #Find timestep and resolution                                                                 
       time=file2read.variables['T']
       ti=time[:]
       dx=Y[1]-Y[0]
       Vtave=np.mean(V,axis = 0)
       Vtave[Vtave==0]=np.nan
       Vzone=np.nansum(Vtave,axis = 2)*dx
       dz=Zp[0:len(Zp)-1]-Zp[1:len(Zp)]
       # Got rid of for loop here (much quicker!!)                                                   
       psi2=np.apply_along_axis(np.multiply,0,Vzone,dz)
       psi=np.cumsum(-psi2[::-1,:],axis=0)
       npad = ((0,1), (0,0))
       # Pad with zeros at bottom                                                                    
       psi = np.pad(psi, pad_width=npad, mode='constant', constant_values=0)
       Psi=(psi/10**6) #Convert to Sv                                                  
       T=ti[0]
       start=int(np.divide(ti[0],(86400*360)))
       end=int(np.divide(ti[-1],(86400*360)))
       return Psi,T

## RMOC ##                                                                                           
def RMOC(file,Y):
    ''' 
    Calculates the RMOC in temperature space
    from the filename and Y from the grid file.
    In future add the option to read in Y from 
    file if not present.
    Outputs Psi in Sv and initial time of file. 
    '''
    file2read = netcdf.NetCDFFile(file,'r')
    lvrho=file2read.variables["LaVH1TH"]
    lvrho=lvrho[:]
    Y=file2read.variables['Yp1']
    Y=Y[:]
    time=file2read.variables['T']
    ti=time[:]
    dx=Y[1]-Y[0] # Find Resolution                                                                   
    VT=np.nansum(lvrho*dx,axis=3) #integrate Vdx along x                                             
    VTfdz=np.cumsum(VT[:,::-1,:],axis=1) #sum up the water column                                    
    psi=np.nanmean(VTfdz[:,::-1,:],axis=0)/10**6 #Time average and put into Sv                       
    #Take the minimum as the strength of ABBW Cells                                                  
    start=int(np.divide(ti[0],(86400*360)))#Find run start and stop times                            
    end=int(np.divide(ti[-1],(86400*360)))
    T=ti[0]
    return psi,T

def find_nearest(array,value):
    ''' What it says on the tin really.
    '''
    idx = (np.abs(array-value)).argmin()
    return array[idx]

def regrid(Variable):
    ''' 
    Regrid from Yp1 to Y
    Stickthis to a numba function for
    for loops!
    '''
    Vc=(Variable[:,0:-1]+Variable[:,1::])/2
    return Vc

def remap(file,Rho,X,Yc,Z,Y,Zp):
    '''
    Remap the RMOCT into Z=space 
    Read in the *all.nc and the X, Y,Z
    from the grid file. Calculates RMOCT
    and gets the temp feild
    Then interpolates both the temperature
    and Psi field to give a smooth result.
    Calclates the MOC as well and iterpolates
    it to the same fine grid so we can take
    the timeaveraged part from the resdiual 
    to give the eddy overturning as well with 
    out the high temporal frequency data - nice! 
    '''
    numba_regrid = autojit()(regrid)
    numba_regrid.func_name = "numba_regrid"
    #Make the file names to read
    a=list(file)
    a[-6]='P'
    a[-5]='s'
    a[-4]='i'
    A=''.join(a)
    if not os.path.isfile(A):
       print "Warning: file "+A+" not found, so i'm skipping this!"
       return
    file2read = netcdf.NetCDFFile(A,'r')
    lvrho=file2read.variables["LaVH1TH"]
    lvrho=lvrho[:]
    dx=Y[1]-Y[0] # Find Resolution                                                                   
    VT=np.sum(lvrho*dx,axis=3) #integrate Vdx along x                                                
    VTfdz=np.cumsum(VT[:,::-1,:],axis=1) #sum up the water column                                    
    psi=np.mean(VTfdz[:,::-1,:],axis=0)/10**6 #Time average and put into Sv and put back in right ord\er                                                                                                  
    y=Y/1000
    nolayers=len(psi[:,1])
    Rho=Rho[0:nolayers]#The layers package bins a layer so adjust for that                           
    file2 = netcdf.NetCDFFile(file,'r')
    Temp=file2.variables['THETA']
    Temp=Temp[:]
    Tav=np.mean(Temp,axis=0)
    V=file2.variables['VVEL']
    V=V[:]*1
    Tavlat=np.mean(Tav,axis=2)
    Yc=Yc/1000
    Vtave=np.mean(V,axis = 0)
    Vtave[Vtave==0]=np.nan
    Vzone=np.nansum(Vtave*dx,axis = 2)
    dz=Zp[0:len(Zp)-1]-Zp[1:len(Zp)]
    # No more super slow forloop!                                                                    
    psi2=np.apply_along_axis(np.multiply,0,Vzone,dz)
    psi3=np.cumsum(-psi2[::-1,:],axis=0)
    npad = ((0,1), (0,0))
    psi4 = np.pad(psi3, pad_width=npad, mode='constant', constant_values=0)
    y =Y/1000
    Psi=psi  #Convert to Sv                                                                          
    Psi2=psi4/10**6
    Psi=numba_regrid(Psi)
    Psi2=numba_regrid(Psi2)
    #Remap Part                                                                                      
    #Expand temperature co-ordinates (30 lvls to 168 lvls)  
    Z2=interp1d(Z,Z,axis=0)
    Znew=np.linspace(int(Z[0]),int(Z[-1]),168)
    Zexp=Z2(Znew)
    T2=interp1d(Z,Tavlat,axis=0)
    Tnew=np.linspace(int(Z[0]),int(Z[-1]),168)
    Texp=T2(Tnew)
    R2=interp1d(Rho,Rho,axis=0)
    Rnew=np.linspace(Rho[0],Rho[-1],168)
    Rexp=R2(Rnew)
    P2=interp1d(Rho,psi,axis=0)
    Pnew=np.linspace(Rho[0],Rho[-1],168)
    Pexp=P2(Pnew)
    Psimap=np.zeros(np.shape(Texp))
    for i in range(len(Yc)):
        for k in range(len(Zexp)):
            D=Texp[k,i]
            if np.isnan(D):
               Psimap[k,i]=np.nan
            else:
                P=Pexp[:,i]
                I=find_nearest(Rexp, D)
                b=np.nonzero(Rexp==I)[0][0]
                Psimap[k,i]=P[b]
    #Psimapped=np.multiply(Psimap,lmav) HFacS is req for topo runs                                   
    Psimap[Psimap==0]=np.nan
    #Now put the MOC into 168 lvls to make the eddy plot                                             
    P3=interp1d(Zp,Psi2,axis=0)
    P3new=np.linspace(Z[0],Zp[-1],168)
    P3exp=P3(P3new)
    Psied=Psimap-P3exp
    return Psimap,Psied,Zexp
       
# Regridding
def regridy(Variable):
    Vc=(Variable[:,:,0:-1,:]+Variable[:,:,1::,:])/2
    return Vc
def regridx(Variable):
    Vc=(Variable[:,:,:,0:-1]+Variable[:,:,:,1::])/2
    return Vc
numba_regridy = autojit()(regridy)
numba_regridy.func_name = "numba_regridy"
numba_regridx = autojit()(regridx)
numba_regridx.func_name = "numba_regridx"
## EKE ##                                                      
def EKEf(U,V,dx):
       '''
        Calculate EKE                                                                               
       '''
       numba_regridy = autojit()(regridy)
       numba_regridy.func_name = "numba_regridy"
       numba_regridx = autojit()(regridx)
       numba_regridx.func_name = "numba_regridx"
       Vc=numba_regridy(V)
       Uc=numba_regridx(U)
       k=0.5*(Uc**2+Vc**2)
       EKE=np.nanmean(k)
       return EKE
