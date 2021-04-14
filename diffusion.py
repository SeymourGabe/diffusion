#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  8 20:23:23 2021

@author: gabeseymour
"""

import numpy as np
import matplotlib.pyplot as plt

def find_nearest(array, value): 
    """
    Function from https://stackoverflow.com/questions/2566412/find-nearest-value-in-numpy-array
    because I'm too lazy to write this myself'
    """
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def defineCoords(T, X, Y = 0, Z = 0):
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

def defineCoords2(dimensions, steps):
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
        #realSpace = np.zeros(Nx) #1D
    else:
        print("3D")
        #realSpace = np.zeros([10*T*X, 10*T*Y, 10*T*Z])
    
    fourierSpace = np.zeros_like(realSpace, complex)


    return realCoords, fourierCoords, realSpace, fourierSpace


def testCoords1D01(realSpace, fourierSpace, T, M):
    """
    Function to verify coordinate system is playing nicely in 1D
    """
    # Test FT of const
    test1 = np.ones_like(realSpace)
    FT1 = np.fft.fft(test1)
    pwr1 = np.abs(FT1)
    N = len(test1)
    f = np.arange(0, N)/T
    #plt.plot(f, pwr1) # Should be an impulse at f = 0
    test2 = np.arange(0, len(realSpace))
    F = 10 # Arbitraty
    test2 = np.cos(F*test2)
    t = np.linspace(0, 100, M*T)
    #plt.plot(t[0:int(len(t)/10)], test2[0:int(len(t)/10)]) # Plotting M periods of cosine
    FT2 = np.fft.fft(test2)
    pwr2 = np.abs(FT2)
    N = len(test2)
    f = np.arange(0, N)/T
    Nyquist = N/(2*T)
    NyquistIndex =  int(f.flat[np.abs(f - Nyquist).argmin()])
    max1 = pwr2[0:NyquistIndex].argmax()
    max2 = pwr2[NyquistIndex:].argmax()
    print("N = ", N)
    print("Nyquist = ", Nyquist)
    print("NyquistIndex = ", NyquistIndex)
    print("f[NyquistIndex] = ", f[NyquistIndex])
    print("Max1 = ", max1)
    print("Max2 = ", max2)
    print("f[max1] = ", f[max1])
    print("f[max2] = ", f[max2])
    print("diffFreq1= ", np.abs(f[NyquistIndex] - f[max1]))
    print("diffFreq1= ", np.abs(f[NyquistIndex] - f[max2]))
    plt.plot(f[:], pwr2[:])
    
def testCoords1D0(realSpace, fourierSpace, X, dx):
    x = np.arange(-X/2, X/2 + dx, dx)
    realSpace1 = 0*realSpace + 1
    fourierSpace1 = np.fft.fft(realSpace1)
    plt1 = np.abs(fourierSpace1)
    N = len(realSpace)
    f = np.linspace(0, X, N)/X
    #f = np.arange(0, X+dx, dx)/X
    print("X/dx = ", X/dx)
    print("len(f) = ", len(f))
    print("len(fourierSpace1) = ", len(fourierSpace1))
    print("len(realSpace1) = ", len(realSpace1))
    plt.plot(f, plt1)
    
    realSpace2 = np.arange(0, len(realSpace))*dx
    realSpace2 = np.cos(2*np.pi* 2 * realSpace2) + np.cos(2*np.pi* 4 * realSpace2)
    print(realSpace2)
    fourierSpace2 = np.fft.fft(realSpace2)
    plt2 = np.abs(fourierSpace2)
    #plt.plot(x, realSpace2)
    plt.plot(f*X, plt2)
    
def testCoords1D(realSpace, fourierSpace, X, dx):
    """
    Tests a set of 1D fourier transforms to make sure our coordinates
    are defined correctly
    """
    
    def cosines(freqs, Nyquist, real, fourier, x):
        for i in range(0, len(freqs)):
            real = real + np.cos(2*np.pi * freqs[i] * x)
            
        fourier = np.fft.fft(real)
        return real, fourier
    
    def gauss(sigma, mu, real, fourier, x): # Still need to fix
        real = 1/(sigma * np.sqrt(2 * np.pi)) * np.exp( - (x - mu)**2 / (2 * sigma**2) )
        fourier = np.fft.fft(real)
        return real, fourier
    
    def window(w, real, fourier, x):
        real = real*0 + 1
        mask = np.abs(x)-w/2 < 0
        real = real*mask
        fourier = np.fft.fft(real)
        return real, fourier
    
    N = len(realSpace)
    xVals = np.arange(-X/2, X/2 + dx, dx)
    f = np.linspace(0, X, N)/X
    const = 0*realSpace + 1
    impulse = np.fft.fft(const)
    #plt.plot(f, impulse)
    cosines, impulses = cosines([0.1, 2, 4], X/2, const*0, f, xVals)
   # plt.plot(xVals, cosines)
    #plt.plot(f*X, np.abs(impulses))
    gauss, gaussFT = gauss(1, 0, const*0, f, xVals)
    window, sinc = window(X/2, const, f, xVals)
    plt.plot(f[0:int(len(f)/2)], np.abs(sinc)[0:int(len(f)/2)])
    #plt.plot(xVals, gauss)
   # plt.plot(f, gaussFT)
    
    
    
    
    

rlCoords, frCoords, rlSpace, frSpace = defineCoords2([10, 5, 0], [.1, .2, 0])
#print("xcrds:\n", xcrds, "\n\n")
print("rlCoords:\n", rlCoords, "\n\n")
print("rlSpace:\n", rlSpace, "\n\n")
print("frSpace:\n", frSpace, "\n\n")
testCoords1D(rlCoords[0], frCoords[0], 10, .1)
#rlspc, frspc = defineCoords(2, 100)
#testCoords1D(rlspc, frspc, 2, 100)