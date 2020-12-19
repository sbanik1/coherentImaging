# -*- coding: utf-8 -*-
"""
This code simulates the blurring caused in Coherent Absorption Imaging
It convolves the light field with the ASF = sqrt(PSF)
Created on Thu July 24 12:30:00 2020
@author: Swarnav Banik
sbanik1@umd.edu
"""
from __future__ import division
import startUp
import GaussianBeam as GB
import OpticalElements as OE
import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt
import os as os

# %% Add path and initialize ##################################################
startUp
if (os.path.isdir(startUp.saveDir) == False):
    os.mkdir(startUp.saveDir)

# %% Input Parameters #########################################################
# Probe and Imaging Properties =================================================
lamda = 589                                     # Wavelength [nm]
NA = 0.1                                       # numerical aperture
w0 = 100/3                                      # Waist of incoming collimated probe beam [um]
P = 1                                           # Total power in probe beam [mW]
Rad_airy = 1.22*lamda*1e-3/(2*NA)               # radius of the airy disk [um]
sigma = 3*(589*10**-9)**2/(2*np.pi)*10**12      # Absorption cross section [um^2]
# Atomic Distribution Properties ==============================================
Rad_ring = 50                                   # radius of ring [um]
Thk_ring = 20                                   # thichkness of ring [um]
N = 300e+3                                      # Number of atoms
#Array Parameters
Nx = 2**8                                       # No. of points along X
Ny = 2**8                                       # No. of points along Y

# %% Common Functions ##############################################
    
def plotInt(num,X,Y,Int,ROI_rad,alpha_X,alpha_Y,**kwargs):
    # Plot the probe beam and how its altered by the atomic distribution n2D
    # Inputs:
    #   num: figure number
    #   X,Y: X and Y as 2D arrays of plane of the beam [um]
    #   Probe: Intensity profile of the probe beam w/o atoms as 2D array 
    #   Atoms: Intensity profile of the probe beam with atoms as 2D array 
    #   n2D: atomic distribution
    #   alpha_X and alpha_Y: Dimensionless constants that help define an ROI
    if (Int.shape != X.shape or X.shape != Y.shape): 
        raise Exception('4Fimaging::plotProbe::X, Y, Probe, Atoms and n2D should have same dimensions.')
    title = 1    
    for key, value in kwargs.items(): 
     if key == 'title':
         title = value  
    [Ny,Nx] = X.shape
    #Plot Data
    fig = plt.figure(num) 
    fig.clf()
    ax1=fig.add_subplot(111)
    c = ax1.pcolor(X[0,:], Y[:,0],Int, cmap='Blues', vmin=np.min(Int), vmax=np.max(Int))
    if title == 1:
        ax1.set_title('Intensity: \n Total Power = {0} $\mu$W'.format(round(np.sum(Int)*1e3,2)))
    fig.colorbar(c, ax=ax1, label = 'Beam Intensity')
    ax1.axis([-alpha_X*ROI_rad, alpha_X*ROI_rad, -alpha_Y*ROI_rad, alpha_Y*ROI_rad,])
    ax1.set_xlabel('X ($\mu$m)')
    ax1.set_ylabel('Y ($\mu$m)')
    ax1.set_aspect('equal')
    plt.tight_layout()
    plt.show()
    return fig

def AtomicDen_Ring(X,Y,rad,thk,N):
    # Generates the 2D atomic density n2D of a ring 
    # Inputs:
    #   X,Y: X and Y as 2D arrays of plane of the atoms
    #   rad, thk: radius and thichkness of ring
    #   N: Number of atoms in the ring
    # Output:
    #   n2D: 2D density of atoms

    if (X.shape != Y.shape): 
        raise Exception('4Fimaging::plotProbe::X, Y should have same dimensions.')
    R = np.sqrt(X**2+Y**2)
    TH = np.arctan2(Y,X)
    n2D = np.ones(np.shape(R))
    n2D[(R >= (rad+thk/2)) | (R <= (rad-thk/2))] = 0
    n2D = N*n2D/(np.sum(n2D))
    return n2D
AtomicDen_RingVec = np.vectorize(AtomicDen_Ring)

# %%  Derived Values  #########################################################
k = 2*np.pi*10**(3)/lamda                       # Wave number [um-1]

# %%  Input Field: Point Object Field #########################################
# Define Array ================================================================
ROI_rad = 1.5*max(w0, Rad_airy)             # radius of ROI, defines small ROI [um]
x = np.linspace(-ROI_rad,ROI_rad,Nx)
y = np.linspace(-ROI_rad,ROI_rad,Ny)
[X1,Y1] = np.meshgrid(x,y)
R1 = np.sqrt(X1**2+Y1**2)
# Generating Data =============================================================
E1 = GB.PointSourceVec(1,R1[:,:],0,w0,k)
I1 = GB.BeamInt(E1,P)
plotInt(1,X1,Y1,I1,ROI_rad,1,1)
# Generating the PSF, ASF = sqrt(PSF) and generating the Image ================
ASF = OE.ASF(X1,Y1,Rad_airy,kind='airy')
def func(x):
    return np.sum(OE.ImageViaPSF(X1, Y1, E1, ASF, norm=x))-np.sum(I1)    
norm = optimize.fsolve(func, 1/np.sum(ASF))
plotInt(2,X1,Y1,norm*ASF,ROI_rad,1,1, title = 0)
I1_blur = OE.ImageViaPSF(X1, Y1, E1, ASF, norm=norm)
plotInt(3,X1,Y1,I1_blur,ROI_rad,1,1)

# %% Input Field: Atomic ring under constant intensity probe ##################
# Define Array ================================================================
w0 = 100;                              # radius of the flat beam beam [um]
ROI_rad = 1.5*max(w0, Rad_airy)         # radius of ROI, defines small ROI [um]
x = np.linspace(-ROI_rad,ROI_rad,Nx)
y = np.linspace(-ROI_rad,ROI_rad,Ny)
[X1,Y1] = np.meshgrid(x,y)
R1 = np.sqrt(X1**2+Y1**2)
#Generating Data ==============================================================
E1 = GB.PointSourceVec(1,R1[:,:],0,w0,k)
I1 = GB.BeamInt(E1,P)
n2D = AtomicDen_Ring(X1,Y1,Rad_ring,Thk_ring,N);
I1 = I1*np.exp(-sigma*n2D)
E1 = np.sqrt(I1)
plotInt(4,X1,Y1,I1,ROI_rad,1,1)
# Generating the PSF, ASF = sqrt(PSF) and generating the Image ================
ASF = OE.ASF(X1,Y1,Rad_airy,kind='airy')
def func(x):
    return np.sum(OE.ImageViaPSF(X1, Y1, E1, ASF, norm=x))-np.sum(I1)    
norm = optimize.fsolve(func, 1/np.sum(ASF))
plotInt(5,X1,Y1,norm*ASF,ROI_rad,1,1,title = 0)
I1_blur = OE.ImageViaPSF(X1, Y1, E1, ASF, norm=norm)
plotInt(6,X1,Y1,I1_blur,ROI_rad,1,1)

# %% Input Field: Atomic ring under gaussian intensity probe ##################
# Define Array ================================================================
w0 = 7000;                              # radius of gaussian beam [um]
ROI_rad = 1.5*max(Rad_ring, Rad_airy)   # radius of ROI, defines small ROI [um]
x = np.linspace(-ROI_rad,ROI_rad,Nx)
y = np.linspace(-ROI_rad,ROI_rad,Ny)
[X1,Y1] = np.meshgrid(x,y)
R1 = np.sqrt(X1**2+Y1**2)
#Generating Data ==============================================================
E1 = GB.PlaneWaveFieldVec(1,R1[:,:],0,w0,k)
I1 = GB.BeamInt(E1,P)
n2D = AtomicDen_Ring(X1,Y1,Rad_ring,Thk_ring,N);
I1 = I1*np.exp(-sigma*n2D)
E1 = np.sqrt(I1)
figure = plotInt(7,X1,Y1,I1,ROI_rad,1,1)
cwd = os.getcwd()
os.chdir(startUp.saveDir)
figure.set_size_inches(10, 10)
figure.savefig('ImageBlur_AtomPlane.png')
os.chdir(cwd)
#Generating the PSF, ASF = sqrt(PSF) and generating the Image =================
ASF = OE.ASF(X1,Y1,Rad_airy,kind='airy')
def func(x):
    return np.sum(OE.ImageViaPSF(X1, Y1, E1, ASF, norm=x))-np.sum(I1)    
norm = optimize.fsolve(func, 1/np.sum(ASF))
plotInt(8,X1,Y1,norm*ASF,ROI_rad,1,1,title = 0)
I1_blur = OE.ImageViaPSF(X1, Y1, E1, ASF, norm=norm)
figure  = plotInt(9,X1,Y1,I1_blur,ROI_rad,1,1)
cwd = os.getcwd()
os.chdir(startUp.saveDir)
figure.set_size_inches(10, 10)
figure.savefig('ImageBlur_ImagePlane.png')
os.chdir(cwd)


# Pixelate the image on Camera ================================================
[X1_cam,Y1_cam,I1_blur_cam,PixSize_cam] = OE.PixelizeImage(I1_blur,X1,Y1,0.3)
plotInt(10,X1_cam,Y1_cam,I1_blur_cam,ROI_rad,1,1)





