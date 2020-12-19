# -*- coding: utf-8 -*-
"""
This code simulates a 4F imaging system.
Created on Thu July 19 12:30:00 2020
@author: Swarnav Banik
sbanik1@umd.edu
"""

from __future__ import division
import startUp
import os as os
import GaussianBeam as GB
import OpticalElements as OE
import numpy as np
from scipy import optimize
from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt

# %% Add path and initialize ##################################################
startUp

# %% Input Parameters #########################################################
# Probe and Imaging Properties ================================================
lamda = 589                                 # Wavelength [nm]
w0 = 30000                                  # Waist of incoming collimated probe beam [um]
P = 1                                       # Total power in probe beam[mW]
f1 = +300                                   # Focal length of first lens [mm]
Ap_rad = 1                                  # radius of aperture [mm]
f2 = +300                                   # Focal length of second lens [mm]
sigma = 3*(589*10**-9)**2/(2*np.pi)*10**12  # Absorption cross section [um^2]
# Atomic Distribution Properties ==============================================
RingRad = 100                               # Radius of atomic ring [um]
RingThk = 20                                # Radius of atomic ring [um]
N = 100e+3                                  # Number of atoms
# Array Parameters ============================================================
Nx = 2**8                                   # No. of point in X
Ny = 2**8                                   # No. of point in Y
Nz = 5                                      # No. of point in Z
# %% Derived Values ###########################################################
k = 2*np.pi*10**(3)/lamda                   # Wave number [um-1]
zmin = -4*GB.RayleighLength(k,w0)           # Min Z value [um]
zmax = 4*GB.RayleighLength(k,w0)            # Max Z value [um]
# Larger Array ==============================================================
ROI_rad_L = 50                              # radius of ROI, defines small ROI [mm]
x = np.linspace(-10*w0,10*w0,Nx)
y = np.linspace(-10*w0,10*w0,Ny)
z = np.linspace(zmin,zmax,Nz)
[X_L,Y_L,Z_L] = np.meshgrid(x,y,z)
R_L = np.sqrt(X_L**2+Y_L**2)
THETA_L = np.arctan2(Y_L,X_L)
# Smaller Array ===============================================================
ROI_rad_S = 500                             # radius of ROI, defines small ROI [um]
x = np.linspace(-ROI_rad_S,ROI_rad_S,Nx)
y = np.linspace(-ROI_rad_S,ROI_rad_S,Ny)
z = np.linspace(zmin,zmax,Nz)
[X_S,Y_S,Z_S] = np.meshgrid(x,y,z)
R_S = np.sqrt(X_S**2+Y_S**2)
THETA_S = np.arctan2(Y_S,X_S)

# %% Common Functions ##############################################

def plotGaussianBeam(num,X,Y,I,ROI_radius,alpha_X,alpha_Y, **kwargs):
    # Plot a gaussian beam
    # Inputs:
    #   num: figure number
    #   X,Y: X and Y as 2D arrays of plane of the beam [mm]
    #   I: Intensity profile of the beam as 2D array 
    #   I_fit: Fitted gaussian beam
    #   params: Fitted gaussian beam params 
    #   alpha_X and alpha_Y: Dimensionless constants that help define an ROI
    
    
    plt.rcParams.update({'font.size': 8})
    for key, value in kwargs.items(): 
        if key == 'Fit':
            I_fit = np.asarray(value)  
            if (Y.shape != I_fit.shape):
                raise Exception('4Fimaging::plotGaussianBeam::I_fit should have same dimensions as X,Y and I.')
        if key == 'params':
            parameter = value
    if (I.shape != X.shape or X.shape != Y.shape): 
        raise Exception('4Fimaging::plotGaussianBeam::X, Y, I and I_fit should have same dimensions.')

    [Ny,Nx] = X.shape
    #Plot Data
    fig = plt.figure(num)
    gs=GridSpec(3,3)
    fig.clf()
    ax1=fig.add_subplot(gs[0:2,0:2]) 
    c = ax1.pcolor(X[0,:]*10**-3, Y[:,0]*10**-3,I, cmap='Blues', vmin=0, vmax=np.max(I))
    ax1.set_title('Total Power = {0} mW'.format(round(np.sum(I),2)))
    fig.colorbar(c, ax=ax1, label = 'Beam Intensity')
    ax1.axis([-alpha_X*ROI_radius, alpha_X*ROI_radius, -alpha_Y*ROI_radius, alpha_Y*ROI_radius])
    ax=fig.add_subplot(gs[2,0:2])
    if len(kwargs) == 0:
        ax.plot(X[np.int(Ny/2),:]*10**-3,I[np.int(Ny/2),:],'.',)
        ax.axis([-alpha_X*ROI_radius, alpha_X*ROI_radius, 0, 1.1*np.max(I)])
    else:
        ax.plot(X[np.int(Ny/2),:]*10**-3,I[np.int(Ny/2),:],'.',X[np.int(Ny/2),:]*10**-3,I_fit[np.int(Ny/2),:],'--')
        ax.axis([-alpha_X*ROI_radius, alpha_X*ROI_radius, 0, 1.1*np.max(I_fit)])
        ax.set_title('X waist = {0} $\mu$m'.format(round(parameter[3],2)))
    ax.set_xlabel('X (mm)')    
    ax=fig.add_subplot(gs[0:2,2])
    if len(kwargs) == 0:
        ax.plot(I[:,np.int(Nx/2)],Y[:,np.int(Nx/2)]*10**-3,'.')
        ax.axis([0, 1.1*np.max(I), -alpha_Y*ROI_radius, alpha_Y*ROI_radius])
    else:
        ax.plot(I[:,np.int(Nx/2)],Y[:,np.int(Nx/2)]*10**-3,'.',I_fit[:,np.int(Nx/2)],Y[:,np.int(Nx/2)]*10**-3,'--')
        ax.axis([0, 1.1*np.max(I_fit), -alpha_Y*ROI_radius, alpha_Y*ROI_radius])
        ax.set_title('Y waist = {0} $\mu$m'.format(round(parameter[4],2)))
    ax.set_ylabel('Y (mm)')        
    plt.tight_layout()
    plt.show()
    
def plotProbe(num,X,Y,Probe,Atoms,n2D,ROI_rad,alpha_X,alpha_Y):
    # Plot the probe beam and how its altered by the atomic distribution n2D
    # Inputs:
    #   num: figure number
    #   X,Y: X and Y as 2D arrays of plane of the beam [mm]
    #   Probe: Intensity profile of the probe beam w/o atoms as 2D array 
    #   Atoms: Intensity profile of the probe beam with atoms as 2D array 
    #   n2D: atomic distribution
    #   alpha_X and alpha_Y: Dimensionless constants that help define an ROI
    if (Probe.shape != X.shape or X.shape != Y.shape or Y.shape != Atoms.shape or Atoms.shape != n2D.shape): 
        raise Exception('4Fimaging::plotProbe::X, Y, Probe, Atoms and n2D should have same dimensions.')
    # if (num.size != 3): 
    #     raise Exception('4Fimaging::plotProbe::A list of 3 figure numbers needed')

    plt.rcParams.update({'font.size': 12})
    [Ny,Nx] = X.shape
    #Plot Data
    fig = plt.figure(num[0]) 
    fig.clf()
    ax1=fig.add_subplot(111)
    c = ax1.pcolor(X[0,:], Y[:,0],Probe, cmap='Blues', vmin=np.min(Probe), vmax=np.max(Probe))
    ax1.set_title('Probe: \n Total Power = {0} $\mu$W'.format(round(np.sum(Probe)*1e3,2)))
    fig.colorbar(c, ax=ax1, label = 'Beam Intensity')
    ax1.axis([-alpha_X*ROI_rad, alpha_X*ROI_rad, -alpha_Y*ROI_rad, alpha_Y*ROI_rad,])
    ax1.set_xlabel('X ($\mu$m)')
    ax1.set_ylabel('Y ($\mu$m)')
    fig = plt.figure(num[1]) 
    fig.clf()
    ax1=fig.add_subplot(111) 
    c = ax1.pcolor(X[0,:], Y[:,0],Atoms, cmap='Blues', vmin=np.min(Probe), vmax=np.max(Probe))
    ax1.set_title('Probe + Atoms: \n Total transmitted power = {0} $\mu$W'.format(round(np.sum(Atoms)*1e3,2)))
    fig.colorbar(c, ax=ax1, label = 'Beam Intensity')
    ax1.set_xlabel('X ($\mu$m)')
    ax1.set_ylabel('Y ($\mu$m)')
    ax1.axis([-alpha_X*ROI_rad, alpha_X*ROI_rad, -alpha_Y*ROI_rad, alpha_Y*ROI_rad,])
    if len(num) == 3 :
        fig = plt.figure(num[2]) 
        fig.clf()
        ax1=fig.add_subplot(111) 
        c = ax1.pcolor(X[0,:], Y[:,0],n2D, cmap='Blues', vmin=0, vmax=np.max(Probe))
        ax1.set_title('Atomic Distribution: \n Total Number = {0}k'.format(round(np.sum(n2D)*1e-3,2)))
        fig.colorbar(c, ax=ax1, label = 'Atomic Density')
        ax1.axis([-alpha_X*ROI_rad, alpha_X*ROI_rad, -alpha_Y*ROI_rad, alpha_Y*ROI_rad,])
        ax1.set_xlabel('X ($\mu$m)')
        ax1.set_ylabel('Y ($\mu$m)')
    plt.tight_layout()
    plt.show()

def AtomicDen_Ring(R,TH,rad,thk,N):
    n2D = np.ones(np.shape(R))
    n2D[(R >= (rad+thk/2)) | (R <= (rad-thk/2))] = 0
    n2D = N*n2D/(np.sum(n2D))
    return n2D
AtomicDen_RingVec = np.vectorize(AtomicDen_Ring)
    
def CircularAperture(E,X,Y,r0,x0,y0):
    # Evaluates the response of a circular aperture 
    # Inputs: E - 2D Field pattern before the the aperture
    #         X,Y - 2D grid representing co-ordinates
    #         [x0, y0] - Relative center of aperure [mm]
    #         r0 - radius of circular aperture [mm]
    # Outputs: E - 2D Field pattern after the the aperture
    if (E.shape != X.shape or X.shape != Y.shape ): 
        raise Exception('4Fimaging::CircularAperture::I, X and Y should have same dimensions.')
    r0 = r0*10**3
    x0 = x0*10**3
    y0 = y0*10**3
    X = X-x0
    Y = Y-y0
    R = np.sqrt(X**2+Y**2)
    E[R>r0] = 0    
    return E

# %% PROBE Field: Input Probe Beam ############################################
#Generating Data ==============================================================
X1 = X_L[:,:,0]
Y1 = Y_L[:,:,0]
E1 = GB.PlaneWaveFieldVec(1,R_L[:,:,0],0,w0,k)
I1 = GB.BeamInt(E1,P)
#Generating Fit ===============================================================
I1 = I1.astype(float)
p0 = [np.max(I1),0,0,w0,w0]
[params1,_] = optimize.curve_fit(GB.GaussBeamFit, [X1,Y1], I1.ravel(), p0=p0)
I1_fit = GB.GaussBeamFit([X1,Y1],*params1).reshape(Ny,Nx)
plotGaussianBeam(1,X1,Y1,I1,ROI_rad_L,1.5,1.5,Fit = I1_fit, params = params1)

# %% PROBE Field Evolution through the 4F SYSTEM ##############################
# Generating Data =============================================================
[E2, X2, Y2] = OE.SphLensAction(E1,X1,Y1,k,f1,FocussedAxis ='NONE')
I2 = GB.BeamInt(E2,P)
E2_ap = CircularAperture(E2,X2,Y2,Ap_rad,0,0)
I2_ap = GB.BeamInt(E2_ap,P)
[E3, X3, Y3] = OE.SphLensAction(E2_ap,X2,Y2,k,f2,FocussedAxis ='XY')
I3 = GB.BeamInt(E3,P)
#Generating Fits ==============================================================
I2 = I2.astype(float)
p0 = [np.max(I2),0,0,1,1]
[params2,_] = optimize.curve_fit(GB.GaussBeamFit, (X2,Y2), I2.ravel(),p0 = p0)
I2_fit = GB.GaussBeamFit((X2,Y2),*params2).reshape(Nx,Ny)
plotGaussianBeam(2,X2,Y2,I2,ROI_rad_L,0.0001,0.0001,Fit = I2_fit, params = params2)
I2_ap = I2_ap.astype(float)
p0 = [np.max(I2_ap),0,0,1,1]
[params2,_] = optimize.curve_fit(GB.GaussBeamFit, (X2,Y2), I2_ap.ravel(),p0 = p0)
I2_ap_fit = GB.GaussBeamFit((X2,Y2),*params2).reshape(Nx,Ny)
plotGaussianBeam(3,X2,Y2,I2_ap,ROI_rad_L,0.0001,0.0001,Fit = I2_ap_fit, params = params2)
I3 = I3.astype(float)
p0 = [np.max(I3),0,0,w0,w0]
[params3,_] = optimize.curve_fit(GB.GaussBeamFit, (X3,Y3), I3.ravel(),p0 = p0)
I3_fit = GB.GaussBeamFit((X3,Y3),*params3).reshape(Nx,Ny)
plotGaussianBeam(4,X3,Y3,I3,ROI_rad_L,1.5,1.5,Fit = I3_fit, params = params3)

# %% PROBE + ATOMS: at the plane of atoms #####################################
# Generating Data =============================================================
Xa = X_S[:,:,0]
Ya = Y_S[:,:,0]
Ea = GB.PlaneWaveFieldVec(1,R_S[:,:,0],0,w0,k)
Ia = GB.BeamInt(Ea,P*ROI_rad_S**2/w0**2)
n2D = AtomicDen_Ring(R_S[:,:,0],THETA_S[:,:,0],RingRad,RingThk,N);
AIa = Ia*np.exp(-sigma*n2D)
P_transmitted = np.sum(AIa)
AEa = np.sqrt(AIa)
plotProbe([5,6,7],Xa,Ya,Ia,AIa,n2D,ROI_rad_S,1,1)

# %% PROBE + ATOMS field evolution via 4F system: Spherical lens ##############
[Eb, Xb, Yb] = OE.SphLensAction(Ea,Xa,Ya,k,f1,FocussedAxis ='NONE')
[Ec, Xc, Yc] = OE.SphLensAction(Eb,Xb,Yb,k,f2,FocussedAxis ='XY')
Ib = GB.BeamInt(Eb,P*ROI_rad_S**2/w0**2)
Ic = GB.BeamInt(Ec,P*ROI_rad_S**2/w0**2)
[AEb, Xb, Yb] = OE.SphLensAction(AEa,Xa,Ya,k,f1,FocussedAxis ='NONE')
[AEc, Xc, Yc] = OE.SphLensAction(AEb,Xb,Yb,k,f2,FocussedAxis ='XY')
AIb = GB.BeamInt(AEb,P_transmitted)
AIc = GB.BeamInt(AEc,P_transmitted)
plotProbe([8,9],Xc,Yc,Ic,AIc,n2D,ROI_rad_S,1,1)


