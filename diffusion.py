#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  8 20:23:23 2021

@author: gabeseymour
"""

import numpy as np
import matplotlib.pyplot as plt

def FFT1(real):
    """
    

    Parameters
    ----------
    real : Real space 1D function

    Returns
    -------
    shifted FFT of real: fourier space representation of the passed function. 
        Shifted as to center the Fourier spectrum

    """
    return np.fft.fftshift(np.fft.fft(real))

def IFFT1(fourier):
    """
    

    Parameters
    ----------
    fourier : 1D Fourier spectrum centered in the middle of the array

    Returns
    -------
    inverse shifted real space function: inverse shifted and IFFT'd function

    """
    return np.fft.ifft(np.fft.ifftshift(fourier))

def find_nearest(array, value): 
    """
    Function from https://stackoverflow.com/questions/2566412/find-nearest-value-in-numpy-array
    because I'm too lazy to write this myself'
    """
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

    """
    Creates a 10*T*(MxNxO) real space numpy array and the corresponding Fourier 
    numpy array where T is the time period and M, N, and O are the number
    of time periods in x, y, and z dimensions
    """
    
    # Create 1D, 2D, or 3D array of specified size
    if Z==0 and Y!=0:
        realSpace = np.zeros([10*T*X, 10*T*Y]) # 2D
    elif Z==0 and Y==0:
        realSpace = np.zeros(10*T*X) #1D
    else:
        realSpace = np.zeros([10*T*X, 10*T*Y, 10*T*Z])
    
    fourierSpace = np.zeros_like(realSpace, complex)


    return realSpace, fourierSpace

def defineCoords(dimensions, steps):
    """
    Creates real space numpy array and the corresponding Fourier 
    numpy array based on physical (dimensions) X, Y, and Z dimensions, (steps) dx, dy, and dz spacial
    steps. 
    """
    X = dimensions[0];
    Y = dimensions[1];
    Z = dimensions[2];
    dx = steps[0];
    dy = steps[1];
    dz = steps[2];
    if dx != 0:
        Nx = int(X/dx)
    else:
        Nx = 1
        
    if dy != 0:
        Ny = int(Y/dy)
    else:
        Ny = 1
        
    if dz != 0:
        Nz = int(Z/dz)
    else:
        Nz = 1
    
    print("X = ", X)
    print("Y = ", Y)
    print("Z = ", Z)
    print("dx = ", dx)
    print("dy = ", dy)
    print("dz = ", dz)
    print("Nx = ", Nx)
    print("Ny = ", Ny)
    print("Nz = ", Nz)
    
    # Create 1D, 2D, or 3D array of specified size
    if Z==0 and Y!=0:
        print("2D")
        xCoords = np.arange(-X/2, X/2+dx, dx) # 2D
        print("xCoords len = ", len(xCoords))
        yCoords = np.arange(-Y/2, Y/2+dy, dy) # 2D
        fxCoords = np.arange(0, X/dx)
        fyCoords = np.arange(0, Y/dy)
        realCoords = [xCoords, yCoords]
        fourierCoords = [fxCoords, fyCoords]
        realSpace = np.zeros([Nx + 1, Ny + 1])
        xRealSpace = np.zeros(Nx)
        yRealSpace = np.zeros(Ny)
        
        print("realSpace shape = ", realSpace.shape)
        
    elif Z==0 and Y==0:
        print("1D")
        xCoords = np.arange(-X/2, X/2+dx, dx) # 1D
        fxCoords = np.arange(0, X/dx+dx)
        fxCoords = fxCoords - fxCoords[-1]/2 # Shift everything over
        realCoords = xCoords
        fourierCoords = fxCoords
        realSpace = np.zeros(Nx+1)
        xRealSpace = np.zeros(Nx)
        #realSpace = np.zeros(Nx) #1D
    else:
        print("3D")
        #realSpace = np.zeros([10*T*X, 10*T*Y, 10*T*Z])
    
    fourierSpace = np.zeros_like(realSpace, complex)


    return realCoords, fourierCoords, realSpace, fourierSpace

    
def testCoords1D(realSpace, fourierSpace, X, dx):
    """
    Tests a set of 1D fourier transforms to make sure our coordinates
    are defined correctly
    """
    
    def cosines(freqs, Nyquist, real, fourier, x):
        for i in range(0, len(freqs)):
            real = real + np.cos(2*np.pi * freqs[i] * x)
            
        fourier = np.fft.fftshift(np.fft.fft(real))
        return real, fourier
    
    def gauss(sigma, mu, real, fourier, x): # Still need to fix
        real = 1/(sigma * np.sqrt(2 * np.pi)) * np.exp( - (x - mu)**2 / (2 * sigma**2) )
        fourier = np.fft.fftshift(real)
        return real, fourier
    
    def window(w, real, fourier, x):
        real = real*0 + 1
        mask = np.abs(x)-w/2 < 0
        real = real*mask
        #fourier = np.fft.fftshift(np.fft.fft(real))
        fourier = FFT1(real)
        #f1 = np.fft.fft(real)
        #ifft = np.fft.ifft(np.fft.fftshift(fourier))
        #ifft = np.fft.ifft(np.fft.ifftshift(fourier))
        ifft = IFFT1(fourier)
        return real, fourier, ifft
    
    N = len(realSpace)
    xVals = np.arange(-X/2, X/2 + dx, dx)
    f = np.linspace(0, X, N)/X
    const = 0*realSpace + 1
    impulse = np.fft.fft(const)
    #plt.plot(f, impulse)
    cosines, impulses = cosines([1], X/2, const*0, f, xVals)
    #plt.plot(xVals, cosines)
   # plt.plot(f*X, np.abs(impulses))
    gauss, gaussFT = gauss(1, 0, const*0, f, xVals)
    window, sinc, iwindow = window(X/4, const, f, xVals)
    plt.plot(f - f[-1]/2, np.abs(sinc))
    #plt.plot(xVals, iwindow)
   # plt.plot(xVals, np.abs(gauss))
    #plt.plot(f, np.abs(gaussFT))
    
    
def diffusion1(species, diffusivity, t_D, f_coords, x_coords):
    """
    

    Parameters
    ----------
    species : Array of floats
        Each element's value corresponds with the number of that species in 
        that real-space coordinate
    diffusivity : Float
        D of the material. In µm^2/s. For my holographic material, the 
        diffusivity is 1 µm^2/s. This is fairly low. For testing, use 10x this
        which is arbitrarily picked
    t_D : float
        Diffusion time
    f_coords : Array of floats
        The frequency coordinate of each element
        
    Returns
    -------
    species : Array of floats
        The species distribution after diffusion

    """
    print("len(f_coords) = ", len(f_coords))
    print("len(species) = ", len(species))
    print("len(f_coords) = ", len(f_coords))
    f_c = 1/np.sqrt(diffusivity*t_D) # Characteristic spatial frequency
    print("f_c = ", f_c)
    H = np.exp(-np.square(np.abs(f_coords)/f_c)) # Transfer function
    #plt.figure(1)
    print("f_coords = ", np.fft.fftshift(f_coords))
    print("H = ", H)
    #plt.plot(f_coords, H)
    #plt.plot(x_coords, species)
    f_species = FFT1(species)
    f_species = f_species/f_species[int(len(f_species)/2-.5)] # Gotta normalize
        # So species at f = 0 is 1 
    #plt.figure(2)
    #plt.plot(f_coords, np.abs(f_species))
    f_species = np.multiply(species, H)
    plt.plot(f_coords, np.abs(f_species))
    species = IFFT1(f_species)
    #plt.plot(x_coords, species)
    
    return species
    

rlCoords, frCoords, rlSpace, frSpace = defineCoords([10, 0, 0], [.1, 0, 0])
#print("xcrds:\n", xcrds, "\n\n")
print("rlCoords:\n", rlCoords, "\n\n")
print("len(rlCoords = ", len(rlCoords))
rlSpace[int(len(rlCoords)/2-.5)] = 1 # Set impulse in center
print("rlCoords:\n", rlCoords, "\n\n")

print("rlSpace:\n", rlSpace, "\n\n")
##plt.plot(rlCoords, rlSpace)
rlSpace = diffusion1(rlSpace, 1, 5, frCoords, rlCoords)
##plt.plot(rlCoords, rlSpace)

##print("frSpace:\n", frSpace, "\n\n")
##testCoords1D(rlCoords[0], frCoords[0], 10, .1)
#rlspc, frspc = defineCoords(2, 100)
#testCoords1D(rlspc, frspc, 2, 100)
