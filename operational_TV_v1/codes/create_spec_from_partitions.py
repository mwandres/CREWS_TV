# -*- coding: utf-8 -*-
"""
Created on Fri Oct 18 17:44:42 2019
Code to create directional wave spectrum from partitions
Minimum 1 partitions and 1 time steps required at this point.
See the test script below for how to use this code.
@author: moritzw
"""
import numpy as np
import numpy.matlib as npm
###############################################################################
######################### Define some functions ###############################
###############################################################################

def normpdf_python(x, mu, sigma):
   return( 1/(sigma*np.sqrt(2*np.pi))*np.exp(-1*(x-mu)**2/ (2*sigma**2) ))

def check_if_nan_and_replace_with_zero(X1,X2,X3,X4):
    if np.isnan(X1):
        X1 = 0.0
        X2 = 0.0
        X3 = 0.0
        X4 = 0.0
    elif X1 < -1:
        X1 = 0.0
        X2 = 0.0
        X3 = 0.0
        X4 = 0.0
    return(X1,X2,X3,X4)

###############################################################################
###################### Create directional spectrum ############################
###############################################################################
      
def directional_distribution( dirs, si, th ):
    # INPUT
    # dirs: vector of directions defining output discretization (degrees)
    # si: directional spread (degrees)
    # th: mean direction (degrees)
    #
    # OUTPUT
    # D: distribution of energy for each direction (adimensional, sum(D*dth(radians))=1)
    #
    
    # NO wrapped because partitions have very small spread
    dirs = dirs.reshape(1,-1)
    concatenated_dirs = np.squeeze([abs(th-dirs).T, 360-abs(th-dirs).T]).T

    ind = np.argmin(concatenated_dirs.sum(1))
    theta2=concatenated_dirs[:,ind]
    dth = 2*np.pi/np.size(dirs)
    D = normpdf_python(0, theta2, si)
    D = D/(np.sum(D))/dth
    D = np.squeeze(D.T)
    return(D)

###############################################################################
####################### Create frequency spectrum #############################
###############################################################################
    
def frequency_spectrum( freqs, hs, tp, gamma, sa, sb):
    # INPUT
    # freqs: vector of frequencies defining output discretization (s^-1)
    # hs: significant wave height (m)
    # tp: peak period (s)
    #
    # OPTIONAL INPUT
    # gamma: peak enhancement parameter (adimensional) 
    # sa: shape parameter for frequencies < peak frequency (adimensional)
    # sb: shape parameter for frequencies > peak frequency (adimensional)
    #    
    # OUTPUT
    # Sf: distribution of energy for each frequency (m^2 * s)
    #    
    
    fp = np.squeeze(1/tp)
    s = np.zeros(np.size(freqs))
    s[freqs<=fp] = sa
    s[freqs>fp] = sb
    
    Agamma = 1-0.287*np.log(gamma)
    
    SfPM =  5/16 * np.square(hs) / (tp**4 * freqs**5) * np.exp(-1.25*(tp * freqs)**(-4))
 
    peakFactor = gamma**np.exp( -0.5 * ((tp*freqs-1)/s)**2 )
    Sf = Agamma * SfPM * peakFactor
    
    return(Sf)
    

###############################################################################
################# Create directional frequency spectrum #######################
###############################################################################
    
def freqdir_spectrum(freqs, dirs, phs, ptp, pgamma, psi, pth):
    # INPUT
    # freqs: vector of frequencies defining output discretization (s^-1)
    # dirs: vector of directions defining output discretization (degrees)
    # partitions: struct with fields hsi, tpi, sii, thi being i the part number
    # phs: matrix of significant wave heights (m), rows=times, cols=partitions
    # ptp: matrix of peak periods (s), rows=times, cols=partitions
    # pgamma: vector of gammas (adimensional)
    # psi: vector of directional spreads (degrees)
    # pth: vector of mean directions (degrees)
    #
    # OUTPUT
    # Sthf: matrix of direction(rows) - frequency(cols) spectrum (m^2 * s)
    # freqsMatrix: matrix of frequencies (s^-1)
    # dirsMatrix: matrix of directions (degrees)
    freqs = freqs[:]
    dirs = dirs[:]
    # Sthf = np.zeros((np.size(dirs), np.size(freqs), phs.shape[0]))
   
    [freqsMatrix, dirsMatrix] = np.meshgrid(freqs, dirs)
    non0= npm.where(phs!=0)
    phs=phs[non0]
    ptp=ptp[non0]
    psi=psi[non0]
    pth=pth[non0]
    sort_index = np.argsort(-phs)
    phs=phs[sort_index]
    ptp=ptp[sort_index]
    psi=psi[sort_index]
    pth=pth[sort_index]
    
    Sthf = np.zeros((np.size(dirs), np.size(freqs)))
    
    # for itime in range(0,phs.shape[0]):
    for ipart in range(0,len(phs)):
        hsi = phs[ipart]
        tpi = ptp[ipart]
        sii = psi[ipart]
        thi = pth[ipart]
        
        hsi,thi,tpi,sii = check_if_nan_and_replace_with_zero(hsi,thi,tpi,sii)
        if hsi == 0:
            Sfi = np.zeros(np.size(freqs))
        else:
            Sfi = frequency_spectrum(freqs, hsi, tpi, gamma = 3.3, sa=0.07, sb=0.09)
      
        Di = directional_distribution(dirs, sii, thi)
        
        Sthfi = np.matlib.repmat(Sfi[:].T,np.size(dirs),1) * np.matlib.repmat(Di,np.size(freqs),1).T
        Sthfi[np.isnan(Sthfi)] = 0
        Sthfi = Sthfi * np.pi/180
        Sthf= Sthf+ Sthfi
        
        # 4*np.sqrt(sum(freqs*0.1*sum(Sthf*10)))    
        
    return(Sthf, freqsMatrix, dirsMatrix)


###############################################################################
############################# Test script #####################################
###############################################################################
### Uncomment the below to test the functions.    
###  
#from matplotlib import pyplot as plt
#def plot_dir_freq_spec(freqs,theta,A):
#    per = 1/freqs
#    theta_rad = theta*np.pi/180
#    [sth,sr]=np.meshgrid(theta_rad,per)
#    r = np.arange(0.0001,1.,0.0005)
#    fig, ax = plt.subplots(1,1, subplot_kw=dict(projection='polar'))
#    ax.set_theta_zero_location("N")
#    ax.set_theta_direction(-1)
#    im = ax.contourf(sth,sr,A,r,cmap = 'jet',vmin = 0,vmax = 0.5) 
#    ax.set_rmax(30)
#    ax.set_rticks([5, 10, 15, 20, 25, 30])
#    v = np.linspace(0, 0.5, 6, endpoint=True)
#    cb = fig.colorbar(im,ticks = v)
#    cb.set_clim(0,1)
#    cb.set_label('Wave Energy'+ r'$[m^2 Hz^{-1} deg^{-1}]$')
#    return()
#
#freqs = np.arange(0.025,1.025,.025)
#dirs = np.arange(0,361,10)  
#pgamma = 3.3
##
#phs = np.array([[2.5,1,1]])
#ptp = np.array([[12,7,17]])
#psi = np.array([[20,60,15]])
#pth = np.array([[220,100,5]])
###phs = np.array([[2.2,1,3],[2.2,1,3]])
###ptp = np.array([[14,7,10],[14,7,10]])
###psi = np.array([[20,40,50],[20,40,50]])
###pth = np.array([[220,100,5],[220,100,100]])
##phs = np.array([[2.2]])
##ptp = np.array([[14]])
##psi = np.array([[20]])
##pth = np.array([[220]])
####
####
#Sthf, freqsMatrix, dirsMatrix = freqdir_spectrum(freqs, dirs, phs, ptp, pgamma, psi, pth)
#plot_dir_freq_spec(freqs,dirs,np.squeeze(Sthf).T)


